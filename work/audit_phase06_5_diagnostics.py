#!/usr/bin/env python3
"""Independent audit for Phase 06.5 diagnostic outputs."""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

from run_oracle_benchmark import load_graph

TRUE = "TRUE_FITNESS_POLICY_DIAGNOSTIC"
SURROGATE = "SURROGATE_HORIZON_DIAGNOSTIC"


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def parse_bool(value):
    return value in (True, "True", "true", "1", 1)


def rank(nodes, values):
    return sorted(nodes, key=lambda node: (values[node], -node), reverse=True)


def audit(trajectories, task_meta, frozen, benchmark_root, expected_tasks):
    errors = []
    graphs = {}
    meta = {(row["assay"], row["task_id"]): row for row in task_meta}
    if len(meta) != expected_tasks:
        errors.append(f"task metadata count {len(meta)} != {expected_tasks}")
    expected_trajectory_rows = expected_tasks * 2 * 5 * 4
    if len(trajectories) != expected_trajectory_rows:
        errors.append(f"trajectory rows {len(trajectories)} != {expected_trajectory_rows}")
    keys = set()
    for row in trajectories:
        label = f"{row['assay']} {row['task_id']} {row['policy_type']} H{row['horizon']} B{row['requested_budget']}"
        key = (row["assay"], row["task_id"], row["policy_type"], int(row["horizon"]), int(row["requested_budget"]))
        if key in keys:
            errors.append(f"{label}: duplicate")
        keys.add(key)
        assay = row["assay"]
        if assay not in graphs:
            graphs[assay] = load_graph(
                benchmark_root / assay / "nodes.csv",
                benchmark_root / assay / "edges.csv",
            )
        nodes, _, adj = graphs[assay]
        task = meta.get((assay, row["task_id"]))
        if task is None:
            errors.append(f"{label}: missing task metadata")
            continue
        path = json.loads(row["executed_path"])
        calibration = set(json.loads(row["initial_observation_nodes"]))
        if path[0] != int(row["start_state"]):
            errors.append(f"{label}: wrong start")
        if len(path) != len(set(path)):
            errors.append(f"{label}: repeated node")
        if any(right not in adj[left] for left, right in zip(path, path[1:])):
            errors.append(f"{label}: non-edge")
        if len(path) - 1 != int(row["actual_budget"]) or int(row["actual_budget"]) > int(row["requested_budget"]):
            errors.append(f"{label}: budget mismatch")
        expected_best = max(nodes[node]["fitness"] for node in path)
        if not math.isclose(expected_best, float(row["best_fitness_seen"]), rel_tol=0.0, abs_tol=1e-10):
            errors.append(f"{label}: best fitness mismatch")
        if row["policy_type"] == TRUE and not parse_bool(row["diagnostic_exact"]):
            errors.append(f"{label}: non-exact true diagnostic")
        observed = set(calibration)
        events = json.loads(row["events"])
        for index, event in enumerate(events, 1):
            current = path[index - 1]
            chosen = path[index]
            legal = sorted(node for node in adj[current] if node not in observed)
            if event["current_before"] != current or event["selected_node"] != chosen:
                errors.append(f"{label} q{index}: event/path mismatch")
                break
            if event["number_legal_neighbors"] != len(legal) or chosen not in legal:
                errors.append(f"{label} q{index}: legal set mismatch")
                break
            true_values = {node: nodes[node]["fitness"] for node in legal}
            true_rank = rank(legal, true_values)
            if event["true_best_neighbor"] != true_rank[0]:
                errors.append(f"{label} q{index}: true best mismatch")
                break
            expected_regret = true_values[true_rank[0]] - true_values[chosen]
            if not math.isclose(expected_regret, event["local_action_regret"], rel_tol=0.0, abs_tol=1e-10):
                errors.append(f"{label} q{index}: local regret mismatch")
                break
            if row["policy_type"] == TRUE:
                if event["predicted_best_neighbor"] != true_rank[0] or event["rank_of_true_best_neighbor"] != 1:
                    errors.append(f"{label} q{index}: true policy rank mismatch")
                    break
            observed.add(chosen)
    # Frozen surrogate H1/H3 must match all three Phase 06 seeds exactly.
    diag = {
        (row["assay"], row["task_id"], int(row["horizon"]), int(row["requested_budget"])): row
        for row in trajectories if row["policy_type"] == SURROGATE and int(row["horizon"]) in (1, 3)
    }
    frozen_pairs = 0
    for row in frozen:
        if (row["assay"], row["task_id"]) not in meta or row["method"] not in {"ADAPTIVE_GREEDY", "ADAPTIVE_LOOKAHEAD"}:
            continue
        horizon = 1 if row["method"] == "ADAPTIVE_GREEDY" else 3
        key = (row["assay"], row["task_id"], horizon, int(row["requested_budget"]))
        candidate = diag.get(key)
        frozen_pairs += 1
        if candidate is None:
            errors.append(f"frozen {key}: missing diagnostic")
        elif json.loads(candidate["executed_path"]) != json.loads(row["executed_path"]):
            errors.append(f"frozen {key} seed {row['seed']}: path mismatch")
        elif not math.isclose(float(candidate["best_fitness_seen"]), float(row["best_fitness_seen"]), rel_tol=0.0, abs_tol=1e-12):
            errors.append(f"frozen {key} seed {row['seed']}: fitness mismatch")
    for key, row in meta.items():
        path = json.loads(row["shortest_target_path"])
        assay = row["assay"]
        nodes, _, adj = graphs[assay]
        if path:
            if path[0] != int(row["start_state"]):
                errors.append(f"{key}: actionable witness wrong start")
            if any(right not in adj[left] for left, right in zip(path, path[1:])):
                errors.append(f"{key}: actionable witness non-edge")
            downhill = sum(nodes[right]["fitness"] < nodes[left]["fitness"] - 1e-12 for left, right in zip(path, path[1:]))
            if parse_bool(row["requires_downhill"]) and downhill == 0:
                errors.append(f"{key}: actionable witness lacks downhill")
    return {
        "status": "PASS" if not errors else "FAIL",
        "tasks": len(meta),
        "trajectory_rows": len(trajectories),
        "frozen_h1_h3_seed_comparisons": frozen_pairs,
        "all_true_rows_exact": all(parse_bool(row["diagnostic_exact"]) for row in trajectories if row["policy_type"] == TRUE),
        "errors": errors,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trajectories", type=Path, required=True)
    parser.add_argument("--task-metadata", type=Path, required=True)
    parser.add_argument("--phase06", type=Path, required=True)
    parser.add_argument("--benchmarks", type=Path, required=True)
    parser.add_argument("--expected-tasks", type=int, default=154)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = audit(
        read_csv(args.trajectories), read_csv(args.task_metadata),
        read_csv(args.phase06), args.benchmarks, args.expected_tasks,
    )
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
