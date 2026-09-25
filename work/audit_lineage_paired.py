#!/usr/bin/env python3
"""Audit the completed Phase 06 paired LINEAGE_WALK benchmark."""
from __future__ import annotations

import argparse
import csv
import io
import json
import math
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from run_oracle_benchmark import load_graph

METHODS = {
    "RANDOM_WALK",
    "ADAPTIVE_GREEDY",
    "ADAPTIVE_LOOKAHEAD",
    "BUDGETED_WALK_ORACLE",
}


def optional_float(value):
    return None if value in (None, "") else float(value)


def optional_int(value):
    return None if value in (None, "") else int(value)


def parse_bool(value):
    if value in (True, "True", "true", "1", 1):
        return True
    if value in (False, "False", "false", "0", 0):
        return False
    return None


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames or [], list(reader)


def audit_rows(rows, benchmark_root):
    errors = []
    graphs = {}
    groups = defaultdict(list)
    tasks = set()
    strata = defaultdict(Counter)

    for index, row in enumerate(rows, 2):
        label = f"line {index} {row.get('assay')} {row.get('task_id')} {row.get('method')} B={row.get('requested_budget')}"
        try:
            assay = row["assay"]
            if assay not in graphs:
                graph_dir = benchmark_root / assay
                graphs[assay] = load_graph(graph_dir / "nodes.csv", graph_dir / "edges.csv")
            nodes, _, adj = graphs[assay]
            method = row["method"]
            task_key = (assay, row["task_id"])
            pair_key = (assay, row["task_id"], int(row["seed"]), int(row["requested_budget"]))
            groups[pair_key].append(row)
            tasks.add(task_key)
            strata[assay]["VALLEY_REQUIRED" if parse_bool(row["valley_required"]) else "NO_VALLEY_REQUIRED"] += 1

            path = json.loads(row["executed_path"])
            fitness_path = json.loads(row["executed_fitness_path"])
            requested = int(row["requested_budget"])
            actual = int(row["actual_budget"])
            start = int(row["start_state"])
            if path[0] != start:
                errors.append(f"{label}: wrong_start")
            if len(path) != len(set(path)):
                errors.append(f"{label}: repeated_node")
            if any(right not in adj[left] for left, right in zip(path, path[1:])):
                errors.append(f"{label}: non_edge_transition")
            if len(path) - 1 != actual or actual > requested:
                errors.append(f"{label}: budget_or_path_length_mismatch")
            expected_fitness = [nodes[node]["fitness"] for node in path]
            if fitness_path != expected_fitness:
                errors.append(f"{label}: fitness_path_mismatch")
            if not parse_bool(row["path_valid"]):
                errors.append(f"{label}: stored_path_invalid")
            if not math.isclose(float(row["final_current_fitness"]), expected_fitness[-1], rel_tol=0.0, abs_tol=1e-12):
                errors.append(f"{label}: final_fitness_mismatch")
            initial_best = float(row["initial_best_fitness"])
            expected_best = max([initial_best, *expected_fitness[1:]])
            if not math.isclose(float(row["best_fitness_seen"]), expected_best, rel_tol=0.0, abs_tol=1e-12):
                errors.append(f"{label}: best_fitness_mismatch")

            oracle = float(row["budgeted_oracle_fitness"])
            denominator = oracle - float(row["start_fitness"])
            epsilon = max(1e-12, 1e-9 * max(abs(denominator), 1.0))
            stored_regret = optional_float(row["normalized_regret"])
            expected_regret = None if denominator <= epsilon else (oracle - expected_best) / denominator
            if stored_regret is None and expected_regret is not None:
                errors.append(f"{label}: missing_normalized_regret")
            elif stored_regret is not None and expected_regret is None:
                errors.append(f"{label}: unexpected_normalized_regret")
            elif stored_regret is not None and not math.isclose(stored_regret, expected_regret, rel_tol=0.0, abs_tol=1e-10):
                errors.append(f"{label}: normalized_regret_mismatch")

            drops = [fitness_path[i] - fitness_path[i + 1] for i in range(len(fitness_path) - 1)]
            downhill = [i + 1 for i, drop in enumerate(drops) if drop > epsilon]
            if int(row["number_downhill_steps"]) != len(downhill):
                errors.append(f"{label}: downhill_count_mismatch")
            first = downhill[0] if downhill else None
            target = float(row["target_fitness"])
            recovery = None
            if first is not None:
                recovery = next((i for i in range(first + 1, len(fitness_path)) if fitness_path[i] >= target - epsilon), None)
            expected_valley = bool(
                parse_bool(row["valley_required"])
                and fitness_path[0] < target - epsilon
                and first is not None
                and recovery is not None
            )
            stored_valley = parse_bool(row["strict_valley_crossing_success"])
            if parse_bool(row["valley_required"]) and stored_valley != expected_valley:
                errors.append(f"{label}: strict_valley_mismatch")
            if optional_int(row["first_downhill_query"]) != first:
                errors.append(f"{label}: first_downhill_mismatch")
            if optional_int(row["recovery_query"]) != recovery:
                errors.append(f"{label}: recovery_mismatch")

            events = json.loads(row["prequery_events"])
            if method in {"ADAPTIVE_GREEDY", "ADAPTIVE_LOOKAHEAD"}:
                expected_counts = list(range(1, len(events) + 1))
                if [event["training_observation_count"] for event in events] != expected_counts:
                    errors.append(f"{label}: retraining_count_mismatch")
                horizon = 1 if method == "ADAPTIVE_GREEDY" else 3
                if int(row["planning_horizon"]) != horizon:
                    errors.append(f"{label}: horizon_mismatch")
                for event in events:
                    planned = event["planned_path"]
                    if not planned or event["selected_node"] != planned[1] or len(planned) - 1 > horizon:
                        errors.append(f"{label}: first_action_or_horizon_violation")
                        break
            elif events and method == "BUDGETED_WALK_ORACLE":
                errors.append(f"{label}: oracle_has_prequery_events")

            if method == "BUDGETED_WALK_ORACLE":
                oracle_path = json.loads(row["budgeted_oracle_path"])
                if path != oracle_path:
                    errors.append(f"{label}: oracle_path_mismatch")
                if not math.isclose(expected_best, oracle, rel_tol=0.0, abs_tol=1e-12):
                    errors.append(f"{label}: oracle_best_mismatch")
        except Exception as exc:
            errors.append(f"{label}: exception:{type(exc).__name__}:{exc}")

    for key, members in groups.items():
        if {row["method"] for row in members} != METHODS:
            errors.append(f"group {key}: method_set_mismatch")
        if len(members) != len(METHODS):
            errors.append(f"group {key}: member_count_mismatch")
        if len({row["initial_observation_nodes"] for row in members}) != 1:
            errors.append(f"group {key}: calibration_mismatch")
        if len({row["budgeted_oracle_fitness"] for row in members}) != 1:
            errors.append(f"group {key}: oracle_fitness_mismatch")

    return {
        "records": len(rows),
        "tasks": len(tasks),
        "paired_groups": len(groups),
        "methods": sorted({row["method"] for row in rows}),
        "seeds": sorted({int(row["seed"]) for row in rows}),
        "budgets": sorted({int(row["requested_budget"]) for row in rows}),
        "strata_record_counts": {assay: dict(counts) for assay, counts in sorted(strata.items())},
        "initially_solved_records": sum(parse_bool(row["initially_solved"]) is True for row in rows),
        "strict_valley_success_records": sum(parse_bool(row["strict_valley_crossing_success"]) is True for row in rows),
        "clipped_records": sum(int(row["actual_budget"]) < int(row["requested_budget"]) for row in rows),
        "errors": errors,
    }


def audit_registry(repo, registry, benchmark_rows):
    old_bytes = subprocess.check_output([
        "git", "-C", str(repo), "show", "HEAD:work/results/experiment_registry.csv",
    ])
    old_reader = csv.DictReader(io.StringIO(old_bytes.decode()))
    old_fields = old_reader.fieldnames or []
    old_rows = list(old_reader)
    _, current = read_csv(registry)
    prefix = current[:len(old_rows)]
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=old_fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows([{field: row.get(field, "") for field in old_fields} for row in prefix])
    phase_rows = [row for row in current if row.get("phase") == "06_lineage_walk"]
    return {
        "prior_rows": len(old_rows),
        "current_rows": len(current),
        "lineage_rows": len(phase_rows),
        "prior_prefix_byte_identical": output.getvalue().encode() == old_bytes,
        "lineage_rows_match_output": len(phase_rows) == len(benchmark_rows),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--benchmarks", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    _, rows = read_csv(args.input)
    result = audit_rows(rows, args.benchmarks)
    result["registry"] = audit_registry(args.repo, args.registry, rows)
    result["status"] = "PASS" if (
        not result["errors"]
        and result["records"] == 7392
        and result["tasks"] == 154
        and result["paired_groups"] == 1848
        and result["registry"]["prior_prefix_byte_identical"]
        and result["registry"]["lineage_rows_match_output"]
    ) else "FAIL"
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
