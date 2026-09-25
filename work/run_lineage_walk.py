#!/usr/bin/env python3
"""Strict Phase 06 lineage-walk benchmark under partial observation."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import time
from collections import deque
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from build_fitness_graph import parse_mutations
from run_oracle_benchmark import load_graph
from run_partial_observation import append_registry, connected_component, read_tasks, task_record

METHODS = (
    "RANDOM_WALK",
    "ADAPTIVE_GREEDY",
    "ADAPTIVE_LOOKAHEAD",
    "BUDGETED_WALK_ORACLE",
)


def lineage_calibration_nodes(start, component, adj, count):
    """Keep current one-step actions hidden while selecting shared training points."""
    if count <= 1:
        return [start]
    excluded = {start, *adj[start]}
    candidates = sorted(
        component - excluded,
        key=lambda node: (
            hashlib.sha256(f"{start}|{node}".encode()).hexdigest(), node,
        ),
    )
    return [start, *candidates[:count - 1]]


def mutation_sets(nodes):
    return {node: parse_mutations(row["mutations"]) for node, row in nodes.items()}


def feature_kernel(left, right, mutations):
    return np.asarray([
        [1.0 + len(mutations[a] & mutations[b]) for b in right]
        for a in left
    ], dtype=float)


def fit_ridge(observed_fitness, mutations, alpha):
    observed = sorted(observed_fitness)
    y = np.asarray([observed_fitness[node] for node in observed], dtype=float)
    kernel = feature_kernel(observed, observed, mutations)
    coefficients = np.linalg.solve(kernel + alpha * np.eye(len(observed)), y)
    return observed, coefficients


def predict_ridge(model, candidates, mutations):
    observed, coefficients = model
    values = feature_kernel(candidates, observed, mutations) @ coefficients
    return dict(zip(candidates, map(float, values)))


def planning_nodes(current, observed, adj, horizon):
    queue = deque([(current, 0)])
    seen = {current}
    nodes = set()
    while queue:
        node, depth = queue.popleft()
        if depth >= horizon:
            continue
        for nxt in sorted(adj[node]):
            if nxt in observed or nxt in seen:
                continue
            seen.add(nxt)
            nodes.add(nxt)
            queue.append((nxt, depth + 1))
    return sorted(nodes)


def legal_neighbors(current, observed, adj):
    return sorted(node for node in adj[current] if node not in observed)


def validate_lineage(path, adj):
    if len(path) != len(set(path)):
        return False, "repeated_node"
    for left, right in zip(path, path[1:]):
        if right not in adj[left]:
            return False, f"non_edge:{left}->{right}"
    return True, ""


def enumerate_predicted_paths(current, observed, adj, predictions, horizon):
    paths = []

    def visit(path):
        if len(path) - 1 >= horizon:
            paths.append(path)
            return
        candidates = [node for node in sorted(adj[path[-1]]) if node not in observed and node not in path]
        if not candidates:
            paths.append(path)
            return
        for node in candidates:
            visit(path + [node])

    visit([current])
    return [path for path in paths if len(path) > 1]


def plan_first_action(current, observed, adj, predictions, horizon, beam_width=16):
    """Maximize predicted terminal fitness and execute only the first action."""
    beam = [[current]]
    completed = []
    for _ in range(horizon):
        expanded = []
        for path in beam:
            candidates = [
                node for node in sorted(adj[path[-1]])
                if node not in observed and node not in path
            ]
            if not candidates:
                completed.append(path)
            else:
                expanded.extend(path + [node] for node in candidates)
        if not expanded:
            break
        expanded.sort(
            key=lambda path: (
                predictions[path[-1]], predictions[path[1]],
                tuple(-node for node in path[1:]),
            ),
            reverse=True,
        )
        beam = expanded[:beam_width]
    paths = completed + beam
    paths = [path for path in paths if len(path) > 1]
    if not paths:
        return None, []
    best_path = max(
        paths,
        key=lambda path: (
            predictions[path[-1]],
            predictions[path[1]],
            tuple(-node for node in path[1:]),
        ),
    )
    return best_path[1], best_path


def budgeted_walk_oracle(start, budget, adj, nodes, blocked=None):
    """Best measured node reachable in <= budget legal steps and a shortest path."""
    blocked = set(blocked or ()) - {start}
    queue = deque([start])
    parent = {start: None}
    distance = {start: 0}
    while queue:
        current = queue.popleft()
        if distance[current] >= budget:
            continue
        for node in sorted(adj[current]):
            if node not in blocked and node not in distance:
                distance[node] = distance[current] + 1
                parent[node] = current
                queue.append(node)
    target = max(
        distance,
        key=lambda node: (nodes[node]["fitness"], -distance[node], -node),
    )
    path = []
    current = target
    while current is not None:
        path.append(current)
        current = parent[current]
    path.reverse()
    return nodes[target]["fitness"], path


def threshold_metrics(curve, threshold, initially_reached):
    if initially_reached:
        return False, None, "threshold_reached_during_calibration"
    hits = [index for index, value in enumerate(curve[1:], 1) if value >= threshold]
    if hits:
        return True, min(hits), ""
    return False, None, "threshold_not_reached"


def lineage_valley_metrics(path, nodes, target_fitness, valley_required, epsilon):
    fitness_path = [nodes[node]["fitness"] for node in path]
    drops = [fitness_path[i] - fitness_path[i + 1] for i in range(len(path) - 1)]
    downhill = [i + 1 for i, drop in enumerate(drops) if drop > epsilon]
    first_downhill = downhill[0] if downhill else None
    recovery = None
    if first_downhill is not None:
        recovery = next(
            (index for index in range(first_downhill + 1, len(path))
             if fitness_path[index] >= target_fitness - epsilon),
            None,
        )
    success = bool(
        valley_required
        and nodes[path[0]]["fitness"] < target_fitness - epsilon
        and first_downhill is not None
        and recovery is not None
    )
    return {
        "executed_fitness_path": fitness_path,
        "number_downhill_steps": len(downhill),
        "max_fitness_drop": max(drops, default=0.0),
        "minimum_path_fitness": min(fitness_path),
        "first_downhill_query": first_downhill,
        "recovery_query": recovery,
        "strict_valley_crossing_success": success if valley_required else None,
        "strict_valley_crossing_success_na_reason": "" if valley_required else "task_not_valley_required",
    }


def prequery_diagnostics(events):
    deployable = [event for event in events if event["predicted_fitness"] is not None]
    if not deployable:
        return None, None, None
    predicted = np.asarray([event["predicted_fitness"] for event in deployable])
    actual = np.asarray([event["true_fitness"] for event in deployable])
    rmse = float(np.sqrt(np.mean((predicted - actual) ** 2)))
    rho = None
    if len(deployable) >= 3 and len(set(predicted)) > 1 and len(set(actual)) > 1:
        value = float(spearmanr(predicted, actual).statistic)
        rho = value if math.isfinite(value) else None
    top_accuracy = float(np.mean([
        event["surrogate_top_choice_true_rank"] == 1 for event in deployable
    ]))
    return rmse, rho, top_accuracy


def normalized_regret(oracle_fitness, best_fitness, start_fitness, epsilon):
    denominator = oracle_fitness - start_fitness
    if denominator <= epsilon:
        return None, "budgeted_oracle_does_not_improve_start"
    return (oracle_fitness - best_fitness) / denominator, ""


def make_record(
    method, assay, task, seed, requested_budget, calibration, executed_path,
    best_curve, events, nodes, component, adj, runtime, alpha, code_version,
    oracle_results=None,
):
    actual_budget = len(executed_path) - 1
    if oracle_results is None:
        oracle_fitness, oracle_path = budgeted_walk_oracle(
            task["start_node"], requested_budget, adj, nodes, blocked=calibration,
        )
    else:
        oracle_fitness, oracle_path = oracle_results[requested_budget]
    global_fitness = max(nodes[node]["fitness"] for node in component)
    start_fitness = nodes[task["start_node"]]["fitness"]
    initial_best = best_curve[0]
    best_fitness = max(best_curve)
    epsilon = max(1e-12, 1e-9 * max(abs(oracle_fitness - start_fitness), 1.0))
    fitnesses = np.asarray([nodes[node]["fitness"] for node in component], dtype=float)
    top1_threshold = float(np.quantile(fitnesses, 0.99))
    top5_threshold = float(np.quantile(fitnesses, 0.95))
    initial_top1 = initial_best >= top1_threshold
    initial_top5 = initial_best >= top5_threshold
    top1_reached, q_top1, q_top1_reason = threshold_metrics(best_curve, top1_threshold, initial_top1)
    top5_reached, q_top5, q_top5_reason = threshold_metrics(best_curve, top5_threshold, initial_top5)
    regret, regret_reason = normalized_regret(
        oracle_fitness, best_fitness, start_fitness, epsilon,
    )
    raw_auc = float(np.trapz(np.asarray(best_curve), dx=1.0)) if actual_budget else 0.0
    oracle_gain = oracle_fitness - start_fitness
    normalized_gain_auc = None
    normalized_gain_auc_reason = "budgeted_oracle_does_not_improve_start"
    if oracle_gain > epsilon:
        gain_curve = (np.asarray(best_curve) - start_fitness) / oracle_gain
        normalized_gain_auc = (
            float(np.trapz(gain_curve, dx=1.0) / actual_budget)
            if actual_budget else 0.0
        )
        normalized_gain_auc_reason = ""
    valley = lineage_valley_metrics(
        executed_path, nodes, task["best_reachable_fitness"],
        task["valley_required"], epsilon,
    )
    rmse, rho, top_accuracy = prequery_diagnostics(events)
    valid, path_error = validate_lineage(executed_path, adj)
    return {
        "phase": "06_lineage_walk",
        "protocol": "LINEAGE_WALK",
        "assay": assay,
        "task_id": task["task_id"],
        "start_state": task["start_node"],
        "seed": seed,
        "method": method,
        "planning_horizon": 1 if method == "ADAPTIVE_GREEDY" else 3 if method == "ADAPTIVE_LOOKAHEAD" else None,
        "lookahead_objective": "predicted_terminal_fitness" if method in {"ADAPTIVE_GREEDY", "ADAPTIVE_LOOKAHEAD"} else "not_applicable",
        "surrogate": "binary_mutation_features_kernel_ridge" if method not in {"RANDOM_WALK", "BUDGETED_WALK_ORACLE"} else "not_applicable",
        "initial_observation_count": len(calibration),
        "initial_observation_nodes": json.dumps(calibration, separators=(",", ":")),
        "initial_best_fitness": initial_best,
        "initial_top1_reached": initial_top1,
        "initial_top5_reached": initial_top5,
        "initially_solved": initial_top5,
        "requested_budget": requested_budget,
        "actual_budget": actual_budget,
        "candidate_pool_size": len(component),
        "start_fitness": start_fitness,
        "best_fitness_seen": best_fitness,
        "final_current_fitness": nodes[executed_path[-1]]["fitness"],
        "best_fitness_gain": best_fitness - start_fitness,
        "budgeted_oracle_fitness": oracle_fitness,
        "budgeted_oracle_path": json.dumps(oracle_path, separators=(",", ":")),
        "clairvoyant_global_ceiling": global_fitness,
        "normalized_regret": regret,
        "normalized_regret_na_reason": regret_reason,
        "auc_best_fitness_vs_query": raw_auc,
        "auc_best_fitness_per_query": raw_auc / actual_budget if actual_budget else 0.0,
        "normalized_best_gain_auc": normalized_gain_auc,
        "normalized_best_gain_auc_na_reason": normalized_gain_auc_reason,
        "top1_threshold": top1_threshold,
        "top5_threshold": top5_threshold,
        "top1_reached_after_calibration": top1_reached,
        "top5_reached_after_calibration": top5_reached,
        "queries_to_top1": q_top1,
        "queries_to_top1_na_reason": q_top1_reason,
        "queries_to_top5": q_top5,
        "queries_to_top5_na_reason": q_top5_reason,
        "valley_required": task["valley_required"],
        "target_fitness": task["best_reachable_fitness"],
        "executed_path": json.dumps(executed_path, separators=(",", ":")),
        "executed_fitness_path": json.dumps(valley.pop("executed_fitness_path"), separators=(",", ":")),
        "path_valid": valid,
        "path_validation_error": path_error,
        **valley,
        "prequery_RMSE": rmse,
        "prequery_Spearman": rho,
        "top_choice_ranking_accuracy": top_accuracy,
        "prequery_events": json.dumps(events, separators=(",", ":")),
        "ridge_alpha": alpha,
        "runtime": runtime,
        "code_version": code_version,
    }


def run_deployable(
    method, assay, task, nodes, adj, mutations, budgets, calibration_count,
    alpha, seed, code_version, oracle_results=None,
):
    start = task["start_node"]
    component = connected_component(start, adj)
    calibration = lineage_calibration_nodes(start, component, adj, min(calibration_count, len(component)))
    observed_fitness = {node: nodes[node]["fitness"] for node in calibration}
    observed = set(calibration)
    current = start
    path = [start]
    best_curve = [max(observed_fitness.values())]
    events = []
    rng_seed = int(hashlib.sha256(
        f"{assay}|{task['task_id']}|{seed}|{method}".encode()
    ).hexdigest()[:16], 16)
    rng = random.Random(rng_seed)
    snapshots = []
    started = time.perf_counter()
    for query_index in range(1, max(budgets) + 1):
        candidates = legal_neighbors(current, observed, adj)
        if not candidates:
            break
        predicted = None
        planned_path = []
        if method == "RANDOM_WALK":
            query = rng.choice(candidates)
        else:
            model = fit_ridge(observed_fitness, mutations, alpha)
            horizon = 1 if method == "ADAPTIVE_GREEDY" else 3
            prediction_nodes = planning_nodes(current, observed, adj, horizon)
            predictions = predict_ridge(model, prediction_nodes, mutations)
            query, planned_path = plan_first_action(
                current, observed, adj, predictions, horizon,
            )
            predicted = predictions[query]
        true_values = sorted(
            ((nodes[node]["fitness"], -node, node) for node in candidates),
            reverse=True,
        )
        true_rank = next(index for index, (_, _, node) in enumerate(true_values, 1) if node == query)
        predicted_rank = None
        surrogate_top_choice_true_rank = None
        if predicted is not None:
            ranked = sorted(candidates, key=lambda node: (predictions[node], -node), reverse=True)
            predicted_rank = ranked.index(query) + 1
            surrogate_top_choice = ranked[0]
            surrogate_top_choice_true_rank = next(
                index for index, (_, _, node) in enumerate(true_values, 1)
                if node == surrogate_top_choice
            )
        event = {
            "query": query_index,
            "current_before": current,
            "selected_node": query,
            "legal_candidate_count": len(candidates),
            "predicted_fitness": predicted,
            "true_fitness": nodes[query]["fitness"],
            "prediction_error": None if predicted is None else predicted - nodes[query]["fitness"],
            "prediction_rank": predicted_rank,
            "true_rank": true_rank,
            "surrogate_top_choice_true_rank": surrogate_top_choice_true_rank,
            "planned_path": planned_path,
            "training_observation_count": len(observed_fitness),
        }
        events.append(event)
        current = query
        path.append(query)
        observed.add(query)
        observed_fitness[query] = nodes[query]["fitness"]
        best_curve.append(max(best_curve[-1], nodes[query]["fitness"]))
        if query_index in budgets:
            snapshots.append(make_record(
                method, assay, task, seed, query_index, calibration, path.copy(),
                best_curve.copy(), events.copy(), nodes, component, adj,
                time.perf_counter() - started, alpha, code_version,
                oracle_results,
            ))
    for budget in budgets:
        if budget > len(path) - 1:
            snapshots.append(make_record(
                method, assay, task, seed, budget, calibration, path.copy(),
                best_curve.copy(), events.copy(), nodes, component, adj,
                time.perf_counter() - started, alpha, code_version,
                oracle_results,
            ))
    return sorted(snapshots, key=lambda row: row["requested_budget"])


def run_budgeted_oracle(
    assay, task, nodes, adj, budgets, calibration_count, seed, code_version,
    oracle_results=None,
):
    start = task["start_node"]
    component = connected_component(start, adj)
    calibration = lineage_calibration_nodes(start, component, adj, min(calibration_count, len(component)))
    rows = []
    started = time.perf_counter()
    for budget in budgets:
        if oracle_results is None:
            _, path = budgeted_walk_oracle(
                start, budget, adj, nodes, blocked=calibration,
            )
        else:
            _, path = oracle_results[budget]
        curve = [max(nodes[node]["fitness"] for node in calibration)]
        for node in path[1:]:
            curve.append(max(curve[-1], nodes[node]["fitness"]))
        rows.append(make_record(
            "BUDGETED_WALK_ORACLE", assay, task, seed, budget, calibration,
            path, curve, [], nodes, component, adj,
            time.perf_counter() - started, 0.0, code_version,
            oracle_results,
        ))
    return rows


def read_stratified_unique_tasks(path, assays, tasks_per_stratum):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected = []
    for assay in assays:
        assay_rows = [row for row in rows if row["assay"] == assay]
        unique = {}
        for row in sorted(assay_rows, key=lambda item: item["task_id"]):
            unique.setdefault((row["start_node"], row["valley_required"]), row)
        for label in ("True", "False"):
            stratum = sorted(
                (row for (_, value), row in unique.items() if value == label),
                key=lambda row: (
                    hashlib.sha256(
                        f"{row['assay']}|{row['start_node']}|{label}".encode()
                    ).hexdigest(),
                    row["task_id"],
                ),
            )
            selected.extend(task_record(row) for row in stratum[:tasks_per_stratum])
    return selected


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmarks", type=Path, required=True)
    parser.add_argument("--tasks", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--assay", action="append", required=True)
    parser.add_argument("--tasks-per-assay", type=int, default=2)
    parser.add_argument("--tasks-per-stratum", type=int)
    parser.add_argument("--task-selection", choices=("smoke_valley_first", "stratified_hash", "stratified_unique_all_seeds"), default="stratified_hash")
    parser.add_argument("--seed", action="append", type=int, default=[])
    parser.add_argument("--budgets", type=int, nargs="+", default=[10, 20, 50, 100])
    parser.add_argument("--calibration-count", type=int, default=1)
    parser.add_argument("--ridge-alpha", type=float, default=1.0)
    parser.add_argument("--code-version", default="phase06-lineage-v1")
    parser.add_argument("--registry", type=Path)
    args = parser.parse_args()

    budgets = sorted(set(args.budgets))
    seeds = args.seed or [1]
    if args.task_selection == "stratified_unique_all_seeds":
        if not args.tasks_per_stratum:
            parser.error("--tasks-per-stratum is required for stratified_unique_all_seeds")
        tasks = read_stratified_unique_tasks(
            args.tasks, args.assay, args.tasks_per_stratum,
        )
    else:
        tasks = read_tasks(
            args.tasks, args.assay, args.tasks_per_assay, args.task_selection,
        )
    records = []
    for task in tasks:
        assay_dir = args.benchmarks / task["assay"]
        nodes, _, adj = load_graph(assay_dir / "nodes.csv", assay_dir / "edges.csv")
        mutations = mutation_sets(nodes)
        start = task["start_node"]
        component = connected_component(start, adj)
        calibration = lineage_calibration_nodes(
            start, component, adj, min(args.calibration_count, len(component)),
        )
        oracle_results = {
            budget: budgeted_walk_oracle(
                start, budget, adj, nodes, blocked=calibration,
            )
            for budget in budgets
        }
        for seed in seeds:
            for method in METHODS:
                if method == "BUDGETED_WALK_ORACLE":
                    records.extend(run_budgeted_oracle(
                        task["assay"], task, nodes, adj, budgets,
                        args.calibration_count, seed, args.code_version,
                        oracle_results,
                    ))
                else:
                    records.extend(run_deployable(
                        method, task["assay"], task, nodes, adj, mutations,
                        budgets, args.calibration_count, args.ridge_alpha,
                        seed, args.code_version, oracle_results,
                    ))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_csv(args.out.with_suffix(".csv"), records)
    payload = {
        "setting": "closed-pool partial observation",
        "protocol": "LINEAGE_WALK",
        "status": "VALIDATION" if len(tasks) <= 4 else "PAIRED_BENCHMARK",
        "task_selection": args.task_selection,
        "methods": list(METHODS),
        "budgets": budgets,
        "records": records,
    }
    args.out.with_suffix(".json").write_text(json.dumps(payload, indent=2) + "\n")
    appended = append_registry(args.registry, records) if args.registry else 0
    print(json.dumps({
        "status": payload["status"], "tasks": len(tasks),
        "methods": len(METHODS), "records": len(records),
        "registry_rows_appended": appended, "out": str(args.out),
    }, indent=2))


if __name__ == "__main__":
    main()
