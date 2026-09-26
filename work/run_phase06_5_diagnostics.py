#!/usr/bin/env python3
"""Phase 06.5 diagnostics on the frozen strict LINEAGE_WALK benchmark."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from collections import defaultdict, deque
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from run_lineage_walk import (
    fit_ridge,
    legal_neighbors,
    mutation_sets,
    normalized_regret,
    planning_nodes,
    predict_ridge,
    validate_lineage,
)
from run_oracle_benchmark import load_graph

TRUE_METHOD = "TRUE_FITNESS_POLICY_DIAGNOSTIC"
SURROGATE_METHOD = "SURROGATE_HORIZON_DIAGNOSTIC"
DEFAULT_HORIZONS = (1, 2, 3, 4, 5)
EPSILON_FLOOR = 1e-12


def parse_bool(value):
    if value in (True, "True", "true", "1", 1):
        return True
    if value in (False, "False", "false", "0", 0):
        return False
    raise ValueError(f"invalid boolean: {value!r}")


def optional_float(value):
    return None if value in (None, "") else float(value)


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_frozen_design(path):
    """Read immutable task/policy records and enforce one provenance tuple per task."""
    rows = read_csv(path)
    by_task = defaultdict(list)
    for row in rows:
        by_task[(row["assay"], row["task_id"])].append(row)
    tasks = []
    errors = []
    for key, members in sorted(by_task.items()):
        provenance_fields = (
            "start_state", "initial_observation_nodes", "valley_required",
            "target_fitness", "start_fitness", "candidate_pool_size",
        )
        for field in provenance_fields:
            if len({row[field] for row in members}) != 1:
                errors.append(f"{key}: non-unique {field}")
        methods = {row["method"] for row in members}
        seeds = sorted({int(row["seed"]) for row in members})
        budgets = sorted({int(row["requested_budget"]) for row in members})
        if methods != {"RANDOM_WALK", "ADAPTIVE_GREEDY", "ADAPTIVE_LOOKAHEAD", "BUDGETED_WALK_ORACLE"}:
            errors.append(f"{key}: method set mismatch")
        if seeds != [1, 2, 3] or budgets != [10, 20, 50, 100]:
            errors.append(f"{key}: seed/budget mismatch")
        first = members[0]
        tasks.append({
            "assay": key[0],
            "task_id": key[1],
            "start": int(first["start_state"]),
            "start_fitness": float(first["start_fitness"]),
            "calibration": json.loads(first["initial_observation_nodes"]),
            "valley_required": parse_bool(first["valley_required"]),
            "target_fitness": float(first["target_fitness"]),
            "candidate_pool_size": int(first["candidate_pool_size"]),
            "seeds": seeds,
            "budgets": budgets,
        })
    if errors:
        raise ValueError("frozen design errors:\n" + "\n".join(errors[:20]))
    if len(tasks) != 154:
        raise ValueError(f"expected 154 frozen tasks, found {len(tasks)}")
    return tasks, rows


def epsilon_for(start_fitness, oracle_fitness):
    return max(EPSILON_FLOOR, 1e-9 * max(abs(oracle_fitness - start_fitness), 1.0))


def path_key(path, values):
    """Phase 06 terminal objective and deterministic first-action/node tie break."""
    return (values[path[-1]], values[path[1]], tuple(-node for node in path[1:]))


def exact_first_action(current, observed, adj, values, horizon, max_expansions):
    """Search legal paths exactly, using an admissible terminal-fitness bound."""
    available = set(values) - set(observed)
    bound_cache = {}

    def terminal_upper_bound(node, remaining):
        key = (node, remaining)
        if key in bound_cache:
            return bound_cache[key]
        if remaining == 0:
            result = values[node]
        else:
            candidates = [nxt for nxt in adj[node] if nxt in available]
            result = max(
                [values[node]] + [
                    terminal_upper_bound(nxt, remaining - 1)
                    for nxt in candidates
                ]
            )
        bound_cache[key] = result
        return result

    best_path = None
    expansions = 0
    capped = False

    def visit(path):
        nonlocal best_path, expansions, capped
        if capped:
            return
        remaining = horizon - (len(path) - 1)
        candidates = [
            node for node in adj[path[-1]]
            if node in available and node not in path
        ]
        if remaining == 0 or not candidates:
            if len(path) > 1 and (best_path is None or path_key(path, values) > path_key(best_path, values)):
                best_path = path.copy()
            return
        ranked = sorted(
            candidates,
            key=lambda node: (
                terminal_upper_bound(node, remaining - 1),
                values[path[1]] if len(path) > 1 else values[node],
                -node,
            ),
            reverse=True,
        )
        for node in ranked:
            upper = terminal_upper_bound(node, remaining - 1)
            if best_path is not None and upper < values[best_path[-1]]:
                continue
            expansions += 1
            if max_expansions and expansions > max_expansions:
                capped = True
                return
            visit(path + [node])

    visit([current])
    return (
        None if best_path is None else best_path[1],
        [] if best_path is None else best_path,
        expansions,
        not capped,
    )


def beam_first_action(current, observed, adj, predictions, horizon, beam_width=16):
    """Frozen Phase 06 beam semantics, generalized to diagnostic horizons."""
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
        expanded.sort(key=lambda path: path_key(path, predictions), reverse=True)
        beam = expanded[:beam_width]
    paths = [path for path in completed + beam if len(path) > 1]
    if not paths:
        return None, []
    best = max(paths, key=lambda path: path_key(path, predictions))
    return best[1], best


def rank_with_tie_break(nodes, values):
    return sorted(nodes, key=lambda node: (values[node], -node), reverse=True)


def local_decision_event(current, candidates, chosen, predictions, nodes, planned_path, query):
    true_values = {node: nodes[node]["fitness"] for node in candidates}
    true_ranking = rank_with_tie_break(candidates, true_values)
    predicted_ranking = rank_with_tie_break(candidates, predictions)
    true_best = true_ranking[0]
    predicted_best = predicted_ranking[0]
    rho = None
    if len(candidates) >= 3:
        x = [predictions[node] for node in candidates]
        y = [true_values[node] for node in candidates]
        if len(set(x)) > 1 and len(set(y)) > 1:
            value = float(spearmanr(x, y).statistic)
            rho = value if math.isfinite(value) else None
    return {
        "query": query,
        "current_before": current,
        "selected_node": chosen,
        "number_legal_neighbors": len(candidates),
        "true_best_neighbor": true_best,
        "predicted_best_neighbor": predicted_best,
        "rank_of_true_best_neighbor": predicted_ranking.index(true_best) + 1,
        "predicted_rank_of_chosen_neighbor": predicted_ranking.index(chosen) + 1,
        "true_rank_of_chosen_neighbor": true_ranking.index(chosen) + 1,
        "true_fitness_of_best_neighbor": true_values[true_best],
        "true_fitness_of_chosen_neighbor": true_values[chosen],
        "predicted_fitness_of_chosen_neighbor": predictions[chosen],
        "top1_action_accuracy": predicted_best == true_best,
        "top3_action_recall": true_best in predicted_ranking[:3],
        "local_action_regret": true_values[true_best] - true_values[chosen],
        "local_spearman": rho,
        "planned_path": planned_path,
    }


def trajectory_snapshots(path, best_curve, budgets):
    rows = {}
    max_steps = len(path) - 1
    for budget in budgets:
        step = min(budget, max_steps)
        rows[budget] = {
            "executed_path": path[:step + 1],
            "best_curve": best_curve[:step + 1],
            "actual_budget": step,
            "best_fitness_seen": max(best_curve[:step + 1]),
        }
    return rows


def recovery_headroom(chosen, blocked, threshold, nodes, adj):
    """Shortest legal continuation from chosen that recovers a fitness threshold."""
    blocked = set(blocked) - {chosen}
    queue = deque([(chosen, 0)])
    seen = {chosen}
    while queue:
        current, distance = queue.popleft()
        if distance > 0 and nodes[current]["fitness"] >= threshold:
            return True, distance
        for nxt in sorted(adj[current]):
            if nxt not in blocked and nxt not in seen:
                seen.add(nxt)
                queue.append((nxt, distance + 1))
    return False, None


def true_branch_terminal(current, chosen, observed, nodes, adj, horizon, max_expansions):
    """Best true terminal value conditional on executing the chosen first action."""
    if horizon <= 1:
        return nodes[chosen]["fitness"], [current, chosen], True
    values = {node: row["fitness"] for node, row in nodes.items()}
    action, suffix, _, exact = exact_first_action(
        chosen, set(observed) | {chosen}, adj, values, horizon - 1, max_expansions,
    )
    if action is None:
        return nodes[chosen]["fitness"], [current, chosen], exact
    return nodes[suffix[-1]]["fitness"], [current, *suffix], exact


def annotate_recovery(events, path, nodes, horizon):
    """Attach downhill/recovery evidence after the complete trajectory is known."""
    fitness_path = [nodes[node]["fitness"] for node in path]
    for event in events:
        query = event["query"]
        before = fitness_path[query - 1]
        chosen = fitness_path[query]
        true_downhill = chosen < before
        recovery = next(
            (index for index in range(query + 1, len(fitness_path))
             if fitness_path[index] >= before),
            None,
        )
        event.update({
            "current_true_fitness": before,
            "was_downhill": true_downhill,
            "true_downhill": true_downhill,
            "eventually_recovered": recovery is not None if true_downhill else None,
            "queries_until_recovery": recovery - query if true_downhill and recovery is not None else None,
            "recovery_query": recovery if true_downhill else None,
            "recovery_within_horizon": (
                recovery is not None and recovery - query <= horizon
                if true_downhill else None
            ),
        })
    return events


def run_true_policy(task, nodes, adj, horizon, max_expansions):
    observed = set(task["calibration"])
    current = task["start"]
    path = [current]
    best_curve = [max(nodes[node]["fitness"] for node in task["calibration"])]
    events = []
    all_exact = True
    started = time.perf_counter()
    for query in range(1, max(task["budgets"]) + 1):
        candidates = legal_neighbors(current, observed, adj)
        if not candidates:
            break
        values = {node: nodes[node]["fitness"] for node in nodes}
        chosen, planned, expansions, exact = exact_first_action(
            current, observed, adj, values, horizon, max_expansions,
        )
        all_exact = all_exact and exact
        if not exact or chosen is None:
            break
        event = local_decision_event(
            current, candidates, chosen,
            {node: nodes[node]["fitness"] for node in candidates},
            nodes, planned, query,
        )
        event.update({
            "search_expansions": expansions,
            "search_exact": exact,
            "predicted_downhill": nodes[chosen]["fitness"] < nodes[current]["fitness"],
            "future_true_gain_within_H": nodes[planned[-1]]["fitness"] - nodes[current]["fitness"],
            "future_predicted_gain_within_H": nodes[planned[-1]]["fitness"] - nodes[current]["fitness"],
            "true_branch_path": planned,
            "true_branch_exact": exact,
        })
        events.append(event)
        current = chosen
        observed.add(chosen)
        path.append(chosen)
        best_curve.append(max(best_curve[-1], nodes[chosen]["fitness"]))
    events = annotate_recovery(events, path, nodes, horizon)
    return trajectory_snapshots(path, best_curve, task["budgets"]), events, all_exact, time.perf_counter() - started


def run_surrogate_policy(
    task, nodes, adj, mutations, horizon, alpha, beam_width,
    true_max_expansions,
):
    observed = set(task["calibration"])
    observed_fitness = {node: nodes[node]["fitness"] for node in task["calibration"]}
    current = task["start"]
    path = [current]
    best_curve = [max(observed_fitness.values())]
    events = []
    started = time.perf_counter()
    for query in range(1, max(task["budgets"]) + 1):
        candidates = legal_neighbors(current, observed, adj)
        if not candidates:
            break
        model = fit_ridge(observed_fitness, mutations, alpha)
        prediction_nodes = planning_nodes(current, observed, adj, horizon)
        predictions = predict_ridge(model, prediction_nodes, mutations)
        predicted_current = predict_ridge(model, [current], mutations)[current]
        chosen, planned = beam_first_action(
            current, observed, adj, predictions, horizon, beam_width,
        )
        if chosen is None:
            break
        event = local_decision_event(
            current, candidates, chosen, predictions, nodes, planned, query,
        )
        event.update({
            "training_observation_count": len(observed_fitness),
            "predicted_downhill": predictions[chosen] < predicted_current,
            "future_predicted_gain_within_H": predictions[planned[-1]] - predicted_current,
            "future_true_gain_within_H": nodes[planned[-1]]["fitness"] - nodes[current]["fitness"],
        })
        if nodes[chosen]["fitness"] < nodes[current]["fitness"]:
            branch_fitness, branch_path, branch_exact = true_branch_terminal(
                current, chosen, observed, nodes, adj, horizon,
                true_max_expansions,
            )
            recovery_exists, recovery_steps = recovery_headroom(
                chosen, observed | {current}, nodes[current]["fitness"], nodes, adj,
            )
            event.update({
                "future_true_gain_within_H": branch_fitness - nodes[current]["fitness"],
                "true_branch_path": branch_path,
                "true_branch_exact": branch_exact,
                "legal_recovery_path_exists": recovery_exists,
                "shortest_legal_recovery_steps": recovery_steps,
            })
        events.append(event)
        current = chosen
        observed.add(chosen)
        observed_fitness[chosen] = nodes[chosen]["fitness"]
        path.append(chosen)
        best_curve.append(max(best_curve[-1], nodes[chosen]["fitness"]))
    events = annotate_recovery(events, path, nodes, horizon)
    return trajectory_snapshots(path, best_curve, task["budgets"]), events, time.perf_counter() - started


def resolve_original_target(task, benchmark_root, nodes):
    """Recover the Phase 04 target node encoded by the frozen task identity."""
    parts = task["task_id"].rsplit(":", 2)
    if len(parts) != 3 or not parts[1].startswith("s") or not parts[2].startswith("t"):
        raise ValueError(f"cannot parse frozen task id: {task['task_id']}")
    seed = int(parts[1][1:])
    index = int(parts[2][1:])
    rows = read_csv(benchmark_root / task["assay"] / f"oracle_compare_seed{seed}.csv")
    row = rows[index]
    if int(row["start_state"]) != task["start"]:
        raise ValueError(f"{task['task_id']}: Phase 04 start mismatch")
    target = int(row["oracle_target"])
    if not math.isclose(nodes[target]["fitness"], task["target_fitness"], rel_tol=0.0, abs_tol=1e-10):
        raise ValueError(f"{task['task_id']}: Phase 04 target fitness mismatch")
    return target


def path_to_target_metrics(task, nodes, adj, target_node=None):
    """Lexicographic graph searches for Phase 06 actionable-valley evidence."""
    start = task["start"]
    target_fitness = task["target_fitness"]
    eps = epsilon_for(task["start_fitness"], target_fitness)
    targets = {target_node} if target_node is not None else {
        node for node, row in nodes.items()
        if row["fitness"] >= target_fitness - eps
    }

    # BFS shortest target path.
    queue = deque([start])
    parent = {start: None}
    target = None
    while queue:
        cur = queue.popleft()
        if cur in targets and cur != start:
            target = cur
            break
        for nxt in sorted(adj[cur]):
            if nxt not in parent:
                parent[nxt] = cur
                queue.append(nxt)
    shortest = []
    if target is not None:
        cur = target
        while cur is not None:
            shortest.append(cur)
            cur = parent[cur]
        shortest.reverse()

    # 0-1 BFS: minimum downhill-edge count to any target.
    dist = {start: 0}
    dq = deque([start])
    while dq:
        cur = dq.popleft()
        for nxt in adj[cur]:
            weight = int(nodes[nxt]["fitness"] < nodes[cur]["fitness"] - eps)
            candidate = dist[cur] + weight
            if candidate < dist.get(nxt, math.inf):
                dist[nxt] = candidate
                (dq.append if weight else dq.appendleft)(nxt)
    target_downhill = [dist[node] for node in targets if node != start and node in dist]
    min_downhill = min(target_downhill) if target_downhill else None

    # Minimax Dijkstra: minimum possible largest single-step drop.
    import heapq
    drop_cost = {start: 0.0}
    heap = [(0.0, start)]
    while heap:
        cost, cur = heapq.heappop(heap)
        if cost != drop_cost[cur]:
            continue
        for nxt in adj[cur]:
            drop = max(0.0, nodes[cur]["fitness"] - nodes[nxt]["fitness"])
            candidate = max(cost, drop)
            if candidate < drop_cost.get(nxt, math.inf):
                drop_cost[nxt] = candidate
                heapq.heappush(heap, (candidate, nxt))
    target_drops = [drop_cost[node] for node in targets if node != start and node in drop_cost]
    min_drop = min(target_drops) if target_drops else None
    reachable = bool(shortest)
    no_monotonic_target = min_downhill is not None and min_downhill > 0
    return {
        "minimum_true_downhill_depth": min_drop,
        "minimum_number_downhill_steps": min_downhill,
        "minimum_path_length_to_superior_target": len(shortest) - 1 if shortest else None,
        "required_horizon": len(shortest) - 1 if shortest and no_monotonic_target else None,
        "reachable_under_strict_lineage": reachable,
        "requires_downhill": no_monotonic_target,
        "shortest_target_path": shortest,
    }


def count_simple_maximal_paths(start, adj, cap):
    count = 0
    stack = [(start, frozenset([start]))]
    while stack:
        cur, visited = stack.pop()
        candidates = [node for node in adj[cur] if node not in visited]
        if not candidates:
            count += 1
            if cap and count >= cap:
                return count, False
        else:
            stack.extend((node, visited | {node}) for node in candidates)
    return count, True


def make_trajectory_rows(task, policy_type, horizon, snapshots, events, runtime, exact, nodes, adj):
    rows = []
    start = task["start_fitness"]
    for budget, snapshot in sorted(snapshots.items()):
        oracle_fitness = max(
            nodes[node]["fitness"]
            for node in shortest_reachable_nodes(task["start"], budget, adj, blocked=set(task["calibration"]) - {task["start"]})
        )
        regret, reason = normalized_regret(
            oracle_fitness, snapshot["best_fitness_seen"], start,
            epsilon_for(start, oracle_fitness),
        )
        performance = None if regret is None else 1.0 - regret
        budget_events = annotate_recovery(
            [dict(event) for event in events[:snapshot["actual_budget"]]],
            snapshot["executed_path"], nodes, horizon,
        )
        rows.append({
            "assay": task["assay"],
            "task_id": task["task_id"],
            "policy_type": policy_type,
            "horizon": horizon,
            "requested_budget": budget,
            "actual_budget": snapshot["actual_budget"],
            "start_state": task["start"],
            "start_fitness": start,
            "initial_observation_nodes": json.dumps(task["calibration"], separators=(",", ":")),
            "valley_required": task["valley_required"],
            "target_fitness": task["target_fitness"],
            "best_fitness_seen": snapshot["best_fitness_seen"],
            "budgeted_oracle_fitness": oracle_fitness,
            "normalized_regret": regret,
            "normalized_regret_na_reason": reason,
            "normalized_performance": performance,
            "executed_path": json.dumps(snapshot["executed_path"], separators=(",", ":")),
            "path_valid": validate_lineage(snapshot["executed_path"], adj)[0],
            "diagnostic_exact": exact,
            "events": json.dumps(budget_events, separators=(",", ":")),
            "runtime_full_trajectory": runtime,
        })
    return rows


def shortest_reachable_nodes(start, budget, adj, blocked=None):
    blocked = set(blocked or ()) - {start}
    distance = {start: 0}
    queue = deque([start])
    while queue:
        cur = queue.popleft()
        if distance[cur] >= budget:
            continue
        for nxt in adj[cur]:
            if nxt not in blocked and nxt not in distance:
                distance[nxt] = distance[cur] + 1
                queue.append(nxt)
    return distance


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase06", type=Path, required=True)
    parser.add_argument("--benchmarks", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--horizons", type=int, nargs="+", default=list(DEFAULT_HORIZONS))
    parser.add_argument("--task-limit-per-assay", type=int)
    parser.add_argument("--true-max-expansions", type=int, default=5_000_000)
    parser.add_argument("--simple-path-cap", type=int, default=5_000_000)
    parser.add_argument("--ridge-alpha", type=float, default=1.0)
    parser.add_argument("--beam-width", type=int, default=16)
    parser.add_argument("--skip-surrogate", action="store_true")
    args = parser.parse_args()

    tasks, _ = load_frozen_design(args.phase06)
    if args.task_limit_per_assay:
        selected = []
        counts = defaultdict(int)
        for task in tasks:
            if counts[task["assay"]] < args.task_limit_per_assay:
                selected.append(task)
                counts[task["assay"]] += 1
        tasks = selected

    trajectories = []
    task_metadata = []
    graph_cache = {}
    for task in tasks:
        assay = task["assay"]
        if assay not in graph_cache:
            graph_cache[assay] = load_graph(
                args.benchmarks / assay / "nodes.csv",
                args.benchmarks / assay / "edges.csv",
            )
        nodes, _, adj = graph_cache[assay]
        mutations = mutation_sets(nodes)
        target_node = resolve_original_target(task, args.benchmarks, nodes)
        metrics = path_to_target_metrics(task, nodes, adj, target_node=target_node)
        component = shortest_reachable_nodes(task["start"], len(nodes), adj)
        path_count, path_count_exact = (None, None)
        if assay.startswith("GCN4"):
            path_count, path_count_exact = count_simple_maximal_paths(
                task["start"], adj, args.simple_path_cap,
            )
        task_metadata.append({
            "assay": assay,
            "task_id": task["task_id"],
            "start_state": task["start"],
            "valley_required": task["valley_required"],
            "candidate_pool_size": task["candidate_pool_size"],
            "reachable_component_size": len(component),
            **{key: (json.dumps(value, separators=(",", ":")) if isinstance(value, list) else value) for key, value in metrics.items()},
            "unique_legal_lineage_paths": path_count,
            "unique_legal_lineage_paths_exact": path_count_exact,
        })
        for horizon in sorted(set(args.horizons)):
            snapshots, events, exact, runtime = run_true_policy(
                task, nodes, adj, horizon, args.true_max_expansions,
            )
            trajectories.extend(make_trajectory_rows(
                task, TRUE_METHOD, horizon, snapshots, events, runtime, exact, nodes, adj,
            ))
            if not args.skip_surrogate:
                snapshots, events, runtime = run_surrogate_policy(
                    task, nodes, adj, mutations, horizon,
                    args.ridge_alpha, args.beam_width,
                    args.true_max_expansions,
                )
                trajectories.extend(make_trajectory_rows(
                    task, SURROGATE_METHOD, horizon, snapshots, events,
                    runtime, True, nodes, adj,
                ))
    write_csv(args.out.with_name(args.out.name + "_trajectories.csv"), trajectories)
    write_csv(args.out.with_name(args.out.name + "_task_metadata.csv"), task_metadata)
    metadata = {
        "protocol": "LINEAGE_WALK",
        "phase": "06.5",
        "frozen_phase06_source": str(args.phase06),
        "tasks": len(tasks),
        "horizons": sorted(set(args.horizons)),
        "true_policy": TRUE_METHOD,
        "true_objective": "true_terminal_fitness_receding_horizon_first_action",
        "surrogate_policy": None if args.skip_surrogate else SURROGATE_METHOD,
        "true_max_expansions_per_decision": args.true_max_expansions,
        "records": len(trajectories),
        "all_true_rows_exact": all(
            parse_bool(row["diagnostic_exact"])
            for row in trajectories if row["policy_type"] == TRUE_METHOD
        ),
    }
    args.out.with_name(args.out.name + "_run.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
