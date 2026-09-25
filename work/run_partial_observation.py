#!/usr/bin/env python3
"""Closed-pool partial-observation benchmark with a shared Ridge surrogate."""
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

METHODS = (
    "RANDOM",
    "SURROGATE_GREEDY",
    "SURROGATE_BEAM",
    "STATIC_SURROGATE_LOOKAHEAD",
    "CLAIRVOYANT_REFERENCE",
)


def connected_component(start, adj):
    stack, seen = [start], {start}
    while stack:
        cur = stack.pop()
        for nxt in adj[cur]:
            if nxt not in seen:
                seen.add(nxt); stack.append(nxt)
    return seen


def calibration_nodes(start, component, adj, count):
    selected, queue, seen = [start], deque([start]), {start}
    while queue and len(selected) < count:
        cur = queue.popleft()
        for nxt in sorted(adj[cur]):
            if nxt in component and nxt not in seen:
                seen.add(nxt); queue.append(nxt); selected.append(nxt)
                if len(selected) >= count:
                    break
    return selected


def mutation_sets(nodes):
    return {node: parse_mutations(row["mutations"]) for node, row in nodes.items()}


def feature_kernel(left, right, mutations):
    # Linear kernel for binary mutation-token features plus an intercept feature.
    return np.asarray([[1.0 + len(mutations[a] & mutations[b]) for b in right] for a in left], dtype=float)


def ridge_predict(observed_fitness, candidates, mutations, alpha):
    observed = sorted(observed_fitness)
    y = np.asarray([observed_fitness[node] for node in observed], dtype=float)
    kernel = feature_kernel(observed, observed, mutations)
    coefficients = np.linalg.solve(kernel + alpha * np.eye(len(observed)), y)
    predictions = feature_kernel(candidates, observed, mutations) @ coefficients
    fitted = kernel @ coefficients
    if len(y) > 1 and np.sum((y - y.mean()) ** 2) > 0:
        r2 = 1.0 - float(np.sum((y - fitted) ** 2) / np.sum((y - y.mean()) ** 2))
        rho = float(spearmanr(y, fitted).statistic)
        if not math.isfinite(rho): rho = None
    else:
        r2, rho = None, None
    return dict(zip(candidates, map(float, predictions))), r2, rho


def frontier(observed, component, adj):
    return sorted({nxt for node in observed for nxt in adj[node] if nxt in component and nxt not in observed})


def choose_random(candidates, rng):
    return rng.choice(candidates)


def choose_greedy(candidates, predictions):
    return max(candidates, key=lambda node: (predictions[node], -node))


def choose_beam(candidates, observed, component, adj, predictions, horizon, width):
    beam = [[node] for node in candidates]
    beam.sort(key=lambda path: (predictions[path[-1]], -path[-1]), reverse=True)
    beam = beam[:width]
    best = list(beam)
    for _ in range(1, horizon):
        expanded = list(beam)
        for path in beam:
            for nxt in sorted(adj[path[-1]]):
                if nxt in component and nxt not in observed and nxt not in path:
                    expanded.append(path + [nxt])
        expanded.sort(key=lambda path: (predictions.get(path[-1], -math.inf), -len(path), -path[-1]), reverse=True)
        beam = expanded[:width]
        best.extend(beam)
    path = max(best, key=lambda item: (predictions.get(item[-1], -math.inf), -len(item), -item[-1]))
    return path[0]


def choose_lookahead(candidates, observed, component, adj, predictions, horizon):
    best_for_first = []
    for first in candidates:
        queue = deque([(first, 1)]); seen_depth = {first: 1}; best_score = predictions[first]
        while queue:
            cur, depth = queue.popleft()
            best_score = max(best_score, predictions.get(cur, -math.inf))
            if depth >= horizon:
                continue
            for nxt in adj[cur]:
                if nxt in component and nxt not in observed and seen_depth.get(nxt, horizon + 1) > depth + 1:
                    seen_depth[nxt] = depth + 1; queue.append((nxt, depth + 1))
        best_for_first.append((best_score, predictions[first], -first, first))
    return max(best_for_first)[-1]


def shortest_next_step_to_best(candidates, observed, component, adj, nodes):
    target = max((node for node in component if node not in observed), key=lambda node: (nodes[node]["fitness"], -node))
    queue, parent = deque([target]), {target: None}
    candidate_set = set(candidates)
    reached = None
    while queue and reached is None:
        cur = queue.popleft()
        if cur in candidate_set:
            reached = cur; break
        for nxt in sorted(adj[cur]):
            if nxt in component and nxt not in observed and nxt not in parent:
                parent[nxt] = cur; queue.append(nxt)
    return reached if reached is not None else max(candidates, key=lambda node: (nodes[node]["fitness"], -node))


def select_deployable_query(method, candidates, observed, component, adj, predictions, rng, horizon, beam_width):
    """Select without access to hidden measured fitness."""
    if method == "RANDOM": return choose_random(candidates, rng)
    if method == "SURROGATE_GREEDY": return choose_greedy(candidates, predictions)
    if method == "SURROGATE_BEAM": return choose_beam(candidates, observed, component, adj, predictions, horizon, beam_width)
    if method == "STATIC_SURROGATE_LOOKAHEAD": return choose_lookahead(candidates, observed, component, adj, predictions, horizon)
    raise ValueError(method)


def select_oracle_query(candidates, observed, component, adj, nodes):
    """Upper-bound selector; this is the only selector allowed hidden fitness."""
    return shortest_next_step_to_best(candidates, observed, component, adj, nodes)


def auc_best_curve(curve):
    return float(np.trapz(np.asarray(curve, dtype=float), dx=1.0)) if len(curve) > 1 else 0.0


def normalized_auc_best_curve(curve, pool_best, epsilon):
    """Mean query-time improvement over calibration, scaled by remaining headroom."""
    if len(curve) <= 1:
        return 0.0
    denominator = pool_best - curve[0]
    if denominator <= epsilon:
        return None
    normalized = (np.asarray(curve, dtype=float) - curve[0]) / denominator
    return float(np.trapz(normalized, dx=1.0) / (len(curve) - 1))


def strict_downhill_path_exists(start, observed, adj, nodes, target, epsilon, max_steps):
    """Whether an observed simple path reaches the target with a strict drop."""
    queue = deque([(start, False, (start,))])
    while queue:
        cur, crossed_downhill, path = queue.popleft()
        if crossed_downhill and nodes[cur]["fitness"] >= target - epsilon:
            return True
        if len(path) - 1 >= max_steps:
            continue
        for nxt in sorted(adj[cur]):
            if nxt not in observed or nxt in path:
                continue
            next_crossed = crossed_downhill or (
                nodes[nxt]["fitness"] < nodes[cur]["fitness"] - epsilon
            )
            queue.append((nxt, next_crossed, path + (nxt,)))
    return False


def snapshot(method, assay, task, seed, requested_budget, calibration, queried, curve, nodes, component, adj, diagnostics, runtime, horizon, beam_width, alpha, code_version, phase):
    fitnesses = sorted(nodes[node]["fitness"] for node in component)
    pool_best, pool_start = fitnesses[-1], nodes[task["start_node"]]["fitness"]
    epsilon = max(1e-12, 1e-9 * max(abs(pool_best - pool_start), 1.0))
    denominator = pool_best - pool_start
    observed = calibration + queried
    observed_set = set(observed)
    calibration_best = curve[0]
    best = max(nodes[node]["fitness"] for node in observed)
    top1_threshold = float(np.quantile(fitnesses, 0.99)); top5_threshold = float(np.quantile(fitnesses, 0.95))
    calibration_top1_reached = calibration_best >= top1_threshold
    calibration_top5_reached = calibration_best >= top5_threshold
    top1_hits = [index for index, value in enumerate(curve[1:], 1) if value >= top1_threshold]
    top5_hits = [index for index, value in enumerate(curve[1:], 1) if value >= top5_threshold]
    queries_to_top1 = None if calibration_top1_reached else (min(top1_hits) if top1_hits else None)
    queries_to_top5 = None if calibration_top5_reached else (min(top5_hits) if top5_hits else None)
    target = task["best_reachable_fitness"]
    target_fitness_reached = best >= target - epsilon
    downhill_path_observed = strict_downhill_path_exists(
        task["start_node"], observed_set, adj, nodes, target, epsilon, horizon
    ) if task["valley_required"] else None
    return {
        "phase": phase,
        "assay": assay,
        "task_id": task["task_id"],
        "start_state": task["start_node"],
        "seed": seed,
        "method": method,
        "surrogate": "binary_mutation_features_kernel_ridge",
        "planning_semantics": (
            "static_surrogate_graph_search"
            if method in {"SURROGATE_BEAM", "STATIC_SURROGATE_LOOKAHEAD"}
            else "not_applicable"
        ),
        "reference_semantics": (
            "clairvoyant_path_to_global_pool_best_not_finite_budget_optimal"
            if method == "CLAIRVOYANT_REFERENCE"
            else "not_applicable"
        ),
        "initial_observation_count": len(calibration),
        "initial_observation_nodes": json.dumps(calibration, separators=(",", ":")),
        "calibration_best_fitness": calibration_best,
        "requested_budget": requested_budget,
        "actual_budget": len(queried),
        "candidate_pool_size": len(component),
        "budget_fraction": len(queried) / max(1, len(component) - len(calibration)),
        "fraction_of_pool_queried": len(queried) / len(component),
        "horizon": horizon,
        "beam_width": beam_width,
        "best_fitness_discovered": best,
        "terminal_fitness": nodes[queried[-1]]["fitness"] if queried else nodes[calibration[-1]]["fitness"],
        "normalized_regret": (pool_best - best) / denominator if denominator > epsilon else None,
        "normalized_regret_na_reason": "" if denominator > epsilon else "pool_best_does_not_improve_start",
        "top1_threshold": top1_threshold,
        "top5_threshold": top5_threshold,
        "calibration_top1_reached": calibration_top1_reached,
        "calibration_top5_reached": calibration_top5_reached,
        "policy_top1_reached": bool(top1_hits) and not calibration_top1_reached,
        "policy_top5_reached": bool(top5_hits) and not calibration_top5_reached,
        "queries_to_top1_after_calibration": queries_to_top1,
        "queries_to_top1_na_reason": (
            "threshold_reached_during_calibration" if calibration_top1_reached
            else "" if queries_to_top1 is not None else "top1_threshold_not_reached"
        ),
        "queries_to_top5_after_calibration": queries_to_top5,
        "queries_to_top5_na_reason": (
            "threshold_reached_during_calibration" if calibration_top5_reached
            else "" if queries_to_top5 is not None else "top5_threshold_not_reached"
        ),
        "auc_best_fitness_raw": auc_best_curve(curve[:len(queried) + 1]),
        "normalized_auc_best_fitness": normalized_auc_best_curve(
            curve[:len(queried) + 1], pool_best, epsilon
        ),
        "normalized_auc_na_reason": "" if pool_best - calibration_best > epsilon else "calibration_at_pool_best",
        "target_fitness_reached": target_fitness_reached,
        "target_fitness_reached_na_reason": "" if task["valley_required"] else "reported_for_all_tasks_not_valley_evidence",
        "strict_downhill_target_path_observed": downhill_path_observed,
        "strict_downhill_target_path_na_reason": "" if task["valley_required"] else "task_not_valley_required",
        "valley_required": task["valley_required"],
        "surrogate_in_sample_R2": diagnostics[-1][0] if diagnostics and method != "CLAIRVOYANT_REFERENCE" else None,
        "surrogate_in_sample_spearman": diagnostics[-1][1] if diagnostics and method != "CLAIRVOYANT_REFERENCE" else None,
        "surrogate_diagnostic_scope": "observed_training_set_in_sample" if method != "CLAIRVOYANT_REFERENCE" else "not_applicable",
        "runtime": runtime,
        "code_version": code_version,
        "ridge_alpha": alpha,
        "queried_nodes": json.dumps(queried, separators=(",", ":")),
    }


def run_policy(method, assay, task, nodes, adj, mutations, budgets, calibration_count, horizon, beam_width, alpha, seed, code_version="phase06-smoke-v3", phase="06_partial_observation_smoke"):
    start = task["start_node"]
    component = connected_component(start, adj)
    calibration = calibration_nodes(start, component, adj, min(calibration_count, len(component)))
    observed, queried = set(calibration), []
    observed_fitness = {node: nodes[node]["fitness"] for node in calibration}
    curve, diagnostics = [max(observed_fitness.values())], []
    rng_seed = int(hashlib.sha256(f"{assay}|{task['task_id']}|{seed}|{method}".encode()).hexdigest()[:16], 16)
    rng = random.Random(rng_seed)
    max_budget = min(max(budgets), max(0, len(component) - len(calibration)))
    snapshots, start_time = [], time.perf_counter()
    for query_index in range(1, max_budget + 1):
        candidates = frontier(observed, component, adj)
        if not candidates:
            break
        prediction_nodes = sorted(component - observed)
        if method == "CLAIRVOYANT_REFERENCE":
            predictions, r2, rho = {}, None, None
            query = select_oracle_query(candidates, observed, component, adj, nodes)
        else:
            predictions, r2, rho = ridge_predict(observed_fitness, prediction_nodes, mutations, alpha)
            query = select_deployable_query(method, candidates, observed, component, adj, predictions, rng, horizon, beam_width)
        diagnostics.append((r2, rho))
        # This is the measurement/reveal boundary. The selected node's measured
        # fitness becomes available only after the policy commits to the query.
        observed.add(query); queried.append(query); observed_fitness[query] = nodes[query]["fitness"]
        curve.append(max(curve[-1], nodes[query]["fitness"]))
        if query_index in budgets:
            snapshots.append(snapshot(method, assay, task, seed, query_index, calibration, queried, curve, nodes, component, adj, diagnostics, time.perf_counter() - start_time, horizon, beam_width, alpha, code_version, phase))
    for budget in budgets:
        if budget > len(queried):
            snapshots.append(snapshot(method, assay, task, seed, budget, calibration, queried, curve, nodes, component, adj, diagnostics, time.perf_counter() - start_time, horizon, beam_width, alpha, code_version, phase))
    snapshots.sort(key=lambda row: row["requested_budget"])
    return snapshots


def task_hash(row):
    key = f"{row['assay']}|{row['start_node']}|{row['valley_required']}"
    return hashlib.sha256(key.encode()).hexdigest()


def task_record(row):
    return {
        "assay": row["assay"], "task_id": row["task_id"],
        "start_node": int(row["start_node"]),
        "best_reachable_fitness": float(row["best_reachable_fitness"]),
        "valley_required": row["valley_required"] == "True",
    }


def read_tasks(path, assays, tasks_per_assay, selection):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected = []
    for assay in assays:
        candidates = [row for row in rows if row["assay"] == assay and int(row["seed"]) == 1]
        if selection == "smoke_valley_first":
            candidates.sort(key=lambda row: (row["valley_required"] != "True", row["task_id"]))
            chosen = candidates[:tasks_per_assay]
        else:
            if tasks_per_assay % 2:
                raise ValueError("stratified_hash requires an even tasks_per_assay")
            per_stratum = tasks_per_assay // 2
            chosen = []
            for label in ("True", "False"):
                stratum = sorted(
                    (row for row in candidates if row["valley_required"] == label),
                    key=lambda row: (task_hash(row), row["task_id"]),
                )
                if len(stratum) < per_stratum:
                    raise ValueError(f"{assay} has only {len(stratum)} tasks in stratum {label}")
                chosen.extend(stratum[:per_stratum])
            chosen.sort(key=lambda row: (row["valley_required"] != "True", task_hash(row)))
        selected.extend(task_record(row) for row in chosen)
    return selected


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def append_registry(path, records):
    """Append Phase 06 snapshots without replacing any earlier-phase row."""
    if path.exists():
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            old_fields = reader.fieldnames or []
            existing = list(reader)
    else:
        old_fields, existing = [], []
    phase06_fields = list(records[0]) if records else []
    fields = old_fields + [field for field in phase06_fields if field not in old_fields]
    key_fields = ("phase", "assay", "task_id", "seed", "method", "requested_budget", "code_version")
    keys = {tuple(row.get(field, "") for field in key_fields) for row in existing if row.get("phase")}
    appended = 0
    for record in records:
        row = {field: record.get(field, "") for field in fields}
        if "dataset" in fields:
            row["dataset"] = row.get("dataset") or "ProteinGym"
        if "task_type" in fields:
            row["task_type"] = row.get("task_type") or "closed_pool_partial_observation"
        key = tuple(str(row.get(field, "")) for field in key_fields)
        if key not in keys:
            existing.append(row); keys.add(key); appended += 1
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(existing)
    return appended


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmarks", type=Path, required=True)
    parser.add_argument("--tasks", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--assay", action="append", required=True)
    parser.add_argument("--tasks-per-assay", type=int, default=1)
    parser.add_argument("--task-selection", choices=("smoke_valley_first", "stratified_hash"), default="smoke_valley_first")
    parser.add_argument("--phase", default="06_partial_observation_smoke")
    parser.add_argument("--seed", action="append", type=int, default=[])
    parser.add_argument("--budgets", type=int, nargs="+", default=[10, 20, 50, 100])
    parser.add_argument("--calibration-count", type=int, default=5)
    parser.add_argument("--horizon", type=int, default=3)
    parser.add_argument("--beam-width", type=int, default=4)
    parser.add_argument("--ridge-alpha", type=float, default=1.0)
    parser.add_argument("--code-version", default="phase06-smoke-v3")
    parser.add_argument("--registry", type=Path)
    args = parser.parse_args()
    seeds = args.seed or [1]
    tasks = read_tasks(args.tasks, args.assay, args.tasks_per_assay, args.task_selection)
    records = []
    for task in tasks:
        assay_dir = args.benchmarks / task["assay"]
        nodes, _, adj = load_graph(assay_dir / "nodes.csv", assay_dir / "edges.csv")
        mutations = mutation_sets(nodes)
        for seed in seeds:
            for method in METHODS:
                records.extend(run_policy(method, task["assay"], task, nodes, adj, mutations, sorted(args.budgets), args.calibration_count, args.horizon, args.beam_width, args.ridge_alpha, seed, args.code_version, args.phase))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_csv(args.out.with_suffix(".csv"), records)
    run_status = "SMOKE_TEST" if args.phase.endswith("smoke") else "FULL_BENCHMARK"
    payload = {
        "setting": "closed-pool partial observation",
        "status": run_status,
        "hidden_fitness_rule": "deployable policies receive only observed fitness and surrogate predictions; CLAIRVOYANT_REFERENCE is isolated and is not labeled finite-budget optimal",
        "calibration": "start plus deterministic breadth-first neighbors, up to calibration_count; excluded from query budget and identical across policies",
        "task_selection": args.task_selection,
        "budgets": sorted(args.budgets), "methods": list(METHODS), "records": records,
    }
    args.out.with_suffix(".json").write_text(json.dumps(payload, indent=2) + "\n")
    appended = append_registry(args.registry, records) if args.registry else 0
    print(json.dumps({"status": run_status, "tasks": len(tasks), "methods": len(METHODS), "records": len(records), "registry_rows_appended": appended, "out": str(args.out)}, indent=2))


if __name__ == "__main__":
    main()
