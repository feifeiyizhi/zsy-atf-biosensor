#!/usr/bin/env python3
"""Analyze Phase 06.5 planning-opportunity and information-bottleneck diagnostics."""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

TRUE = "TRUE_FITNESS_POLICY_DIAGNOSTIC"
SURROGATE = "SURROGATE_HORIZON_DIAGNOSTIC"
HORIZONS = (1, 2, 3, 4, 5)
BUDGETS = (10, 20, 50, 100)


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def optional_float(value):
    return None if value in (None, "") else float(value)


def parse_bool(value):
    if value in (True, "True", "true", "1", 1):
        return True
    if value in (False, "False", "false", "0", 0):
        return False
    return None


def mean_defined(values):
    values = [value for value in values if value is not None and math.isfinite(value)]
    return float(np.mean(values)) if values else None


def percentile(values, q):
    values = sorted(values)
    if not values:
        return None
    position = (len(values) - 1) * q
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return float(values[low])
    return float(values[low] * (high - position) + values[high] * (position - low))


def bootstrap_ci(task_values, seed, iterations):
    values = [value for value in task_values if value is not None and math.isfinite(value)]
    if not values:
        return None, None
    rng = np.random.default_rng(seed)
    estimates = np.mean(rng.choice(values, size=(iterations, len(values)), replace=True), axis=1)
    return percentile(estimates, 0.025), percentile(estimates, 0.975)


def macro_bootstrap_ci(subset, field, seed, iterations):
    by_assay = defaultdict(list)
    for row in subset:
        value = row[field]
        if value is not None and math.isfinite(value):
            by_assay[row["assay"]].append(value)
    if not by_assay:
        return None, None
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(iterations):
        assay_estimates = []
        for values in by_assay.values():
            sample = rng.choice(values, size=len(values), replace=True)
            assay_estimates.append(float(np.mean(sample)))
        estimates.append(float(np.mean(assay_estimates)))
    return percentile(estimates, 0.025), percentile(estimates, 0.975)


def local_summary(events):
    if not events:
        return {key: None for key in (
            "top1_action_accuracy", "top3_action_recall", "mean_rank_of_true_best",
            "mean_local_action_regret", "mean_local_spearman", "mean_legal_neighbors",
        )}
    return {
        "top1_action_accuracy": mean_defined([float(event["top1_action_accuracy"]) for event in events]),
        "top3_action_recall": mean_defined([float(event["top3_action_recall"]) for event in events]),
        "mean_rank_of_true_best": mean_defined([event["rank_of_true_best_neighbor"] for event in events]),
        "mean_local_action_regret": mean_defined([event["local_action_regret"] for event in events]),
        "mean_local_spearman": mean_defined([event["local_spearman"] for event in events]),
        "mean_legal_neighbors": mean_defined([event["number_legal_neighbors"] for event in events]),
    }


def opportunity_label(value, epsilon, exact):
    if not exact or value is None:
        return "AMBIGUOUS"
    if value > epsilon:
        return "POSITIVE_PLANNING_OPPORTUNITY"
    if value < -epsilon:
        return "NEGATIVE_DEEPER_PLANNING"
    return "NO_CLEAR_PLANNING_OPPORTUNITY"


def downhill_failure(event, horizon):
    if not event.get("true_downhill"):
        return None
    true_gain = event.get("future_true_gain_within_H")
    predicted_gain = event.get("future_predicted_gain_within_H")
    recovered = event.get("eventually_recovered")
    within = event.get("recovery_within_horizon")
    recovery_exists = event.get("legal_recovery_path_exists")
    recovery_steps = event.get("shortest_legal_recovery_steps")
    if true_gain is not None and true_gain >= 0 and within:
        return "USEFUL_VALLEY_RECOVERED"
    if predicted_gain is not None and predicted_gain > 0 and (true_gain is None or true_gain <= 0) and not recovered:
        return "FALSE_VALLEY_PREDICTION"
    if true_gain is not None and true_gain > 0 and not recovered:
        return "SURROGATE_RECOVERY_FAILURE"
    if recovered and not within:
        return "HORIZON_TOO_SHORT"
    if recovery_exists and recovery_steps is not None and recovery_steps > horizon:
        return "HORIZON_TOO_SHORT"
    if not recovered and recovery_exists is False:
        return "GRAPH_BLOCKED"
    return "UNRESOLVED"


def build_task_rows(trajectories, task_meta, frozen):
    trajectory = {
        (row["assay"], row["task_id"], row["policy_type"], int(row["horizon"]), int(row["requested_budget"])): row
        for row in trajectories
    }
    metadata = {(row["assay"], row["task_id"]): row for row in task_meta}
    analysis_tasks = set(metadata)
    frozen_map = defaultdict(dict)
    for row in frozen:
        if (
            (row["assay"], row["task_id"]) in analysis_tasks
            and row["method"] in {"ADAPTIVE_GREEDY", "ADAPTIVE_LOOKAHEAD"}
        ):
            frozen_map[(row["assay"], row["task_id"], int(row["seed"]), int(row["requested_budget"]))][row["method"]] = row
    output = []
    for pair_key, methods in sorted(frozen_map.items()):
        assay, task_id, seed, budget = pair_key
        if set(methods) != {"ADAPTIVE_GREEDY", "ADAPTIVE_LOOKAHEAD"}:
            raise ValueError(f"incomplete frozen pair {pair_key}")
        meta = metadata[(assay, task_id)]
        true_rows = {h: trajectory[(assay, task_id, TRUE, h, budget)] for h in HORIZONS}
        surrogate_rows = {h: trajectory[(assay, task_id, SURROGATE, h, budget)] for h in HORIZONS}
        true_h1, true_h3 = true_rows[1], true_rows[3]
        frozen_h1, frozen_h3 = methods["ADAPTIVE_GREEDY"], methods["ADAPTIVE_LOOKAHEAD"]
        absolute_opportunity = float(true_h3["best_fitness_seen"]) - float(true_h1["best_fitness_seen"])
        performance_h1 = optional_float(true_h1["normalized_performance"])
        performance_h3 = optional_float(true_h3["normalized_performance"])
        normalized_opportunity = None if performance_h1 is None or performance_h3 is None else performance_h3 - performance_h1
        oracle_gain = float(true_h1["budgeted_oracle_fitness"]) - float(true_h1["start_fitness"])
        epsilon = max(1e-12, 1e-9 * max(abs(oracle_gain), 1.0))
        surrogate_abs_gain = float(frozen_h3["best_fitness_seen"]) - float(frozen_h1["best_fitness_seen"])
        frozen_r1 = optional_float(frozen_h1["normalized_regret"])
        frozen_r3 = optional_float(frozen_h3["normalized_regret"])
        realized_gain = None if frozen_r1 is None or frozen_r3 is None else frozen_r1 - frozen_r3
        capture = None
        capture_reason = "planning_opportunity_approximately_zero_or_na"
        if normalized_opportunity is not None and abs(normalized_opportunity) > epsilon and realized_gain is not None:
            capture = realized_gain / normalized_opportunity
            capture_reason = ""
        actionable = bool(
            parse_bool(frozen_h1["valley_required"])
            and parse_bool(meta["requires_downhill"])
            and parse_bool(meta["reachable_under_strict_lineage"])
            and int(meta["minimum_path_length_to_superior_target"]) <= budget
        )
        valley_status = (
            "PHASE06_ACTIONABLE_VALLEY" if actionable
            else "NON_VALLEY" if not parse_bool(frozen_h1["valley_required"])
            else "OTHER_UNRESOLVED"
        )
        h1_events = json.loads(surrogate_rows[1]["events"])
        h3_events = json.loads(surrogate_rows[3]["events"])
        local_h1 = local_summary(h1_events)
        local_h3 = local_summary(h3_events)
        failure_counts = Counter(
            failure for failure in (downhill_failure(event, 3) for event in h3_events)
            if failure is not None
        )
        row = {
            "assay": assay,
            "task_id": task_id,
            "seed": seed,
            "budget": budget,
            "valley_required": parse_bool(frozen_h1["valley_required"]),
            "valley_status": valley_status,
            "phase06_actionable_valley": actionable,
            "start_fitness": float(frozen_h1["start_fitness"]),
            "TRUE_H1": float(true_h1["best_fitness_seen"]),
            "TRUE_H3": float(true_h3["best_fitness_seen"]),
            "planning_opportunity_H3_absolute": absolute_opportunity,
            "planning_opportunity_H3": normalized_opportunity,
            "planning_opportunity_label": opportunity_label(
                normalized_opportunity if normalized_opportunity is not None else absolute_opportunity,
                epsilon,
                all(parse_bool(true_rows[h]["diagnostic_exact"]) for h in HORIZONS),
            ),
            "SURROGATE_H1": float(frozen_h1["best_fitness_seen"]),
            "SURROGATE_H3": float(frozen_h3["best_fitness_seen"]),
            "realized_gain_H3_absolute": surrogate_abs_gain,
            "realized_gain_H3": realized_gain,
            "opportunity_capture_ratio": capture,
            "opportunity_capture_ratio_na_reason": capture_reason,
            "normalized_regret_H1": frozen_r1,
            "normalized_regret_H3": frozen_r3,
            "budgeted_oracle_fitness": float(true_h1["budgeted_oracle_fitness"]),
            "minimum_true_downhill_depth": optional_float(meta["minimum_true_downhill_depth"]),
            "minimum_number_downhill_steps": int(meta["minimum_number_downhill_steps"]) if meta["minimum_number_downhill_steps"] else None,
            "minimum_path_length_to_superior_target": int(meta["minimum_path_length_to_superior_target"]) if meta["minimum_path_length_to_superior_target"] else None,
            "required_horizon": int(meta["required_horizon"]) if meta["required_horizon"] else None,
            "budget_feasible": int(meta["minimum_path_length_to_superior_target"]) <= budget if meta["minimum_path_length_to_superior_target"] else False,
            "reachable_under_strict_lineage": parse_bool(meta["reachable_under_strict_lineage"]),
            "candidate_pool_size": int(meta["candidate_pool_size"]),
            "reachable_component_size": int(meta["reachable_component_size"]),
            "unique_legal_lineage_paths": int(meta["unique_legal_lineage_paths"]) if meta["unique_legal_lineage_paths"] else None,
            "unique_legal_lineage_paths_exact": parse_bool(meta["unique_legal_lineage_paths_exact"]),
            **{f"H1_{key}": value for key, value in local_h1.items()},
            **{f"H3_{key}": value for key, value in local_h3.items()},
            "adaptive_h3_downhill_actions": sum(failure_counts.values()),
            "downhill_failure_mode_counts": json.dumps(dict(sorted(failure_counts.items())), separators=(",", ":")),
            "all_true_horizons_exact": all(parse_bool(true_rows[h]["diagnostic_exact"]) for h in HORIZONS),
        }
        for h in HORIZONS:
            row[f"TRUE_H{h}_performance"] = optional_float(true_rows[h]["normalized_performance"])
            row[f"SURROGATE_H{h}_performance"] = optional_float(surrogate_rows[h]["normalized_performance"])
            row[f"TRUE_H{h}_best"] = float(true_rows[h]["best_fitness_seen"])
            row[f"SURROGATE_H{h}_best"] = float(surrogate_rows[h]["best_fitness_seen"])
            row[f"TRUE_H{h}_actual_budget"] = int(true_rows[h]["actual_budget"])
            row[f"SURROGATE_H{h}_actual_budget"] = int(surrogate_rows[h]["actual_budget"])
        output.append(row)
    return output


def collapse_tasks(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["assay"], row["task_id"], row["budget"])].append(row)
    output = []
    for key, members in grouped.items():
        base = dict(members[0])
        for field in ("realized_gain_H3", "realized_gain_H3_absolute", "normalized_regret_H1", "normalized_regret_H3"):
            base[field] = mean_defined([member[field] for member in members])
        output.append(base)
    return output


def build_assay_summary(task_rows, iterations):
    tasks = collapse_tasks(task_rows)
    output = []
    strata = ("ALL_TASKS", "PHASE06_ACTIONABLE_VALLEY", "NON_VALLEY", "OTHER_UNRESOLVED")
    scopes = ["ALL_ASSAYS", *sorted({row["assay"] for row in tasks})]
    policies = ((TRUE, "TRUE"), (SURROGATE, "SURROGATE"))
    seed_counter = 81000
    for scope in scopes:
        for budget in BUDGETS:
            for stratum in strata:
                subset = [
                    row for row in tasks
                    if row["budget"] == budget
                    and (scope == "ALL_ASSAYS" or row["assay"] == scope)
                    and (stratum == "ALL_TASKS" or row["valley_status"] == stratum)
                ]
                if not subset:
                    continue
                for policy_name, prefix in policies:
                    for horizon in HORIZONS:
                        values = [row[f"{prefix}_H{horizon}_performance"] for row in subset]
                        ci_low, ci_high = bootstrap_ci(values, seed_counter, iterations)
                        seed_counter += 1
                        assay_means = [
                            mean_defined([row[f"{prefix}_H{horizon}_performance"] for row in subset if row["assay"] == assay])
                            for assay in sorted({row["assay"] for row in subset})
                        ]
                        delta_field = f"{prefix}_H{horizon}_vs_H1_delta"
                        for row in subset:
                            h_value = row[f"{prefix}_H{horizon}_performance"]
                            h1_value = row[f"{prefix}_H1_performance"]
                            row[delta_field] = None if h_value is None or h1_value is None else h_value - h1_value
                        delta_values = [row[delta_field] for row in subset]
                        delta_ci_low, delta_ci_high = bootstrap_ci(
                            delta_values, seed_counter + 100000, iterations,
                        )
                        macro_ci_low, macro_ci_high = macro_bootstrap_ci(
                            subset, delta_field, seed_counter + 200000, iterations,
                        )
                        assay_delta_means = [
                            mean_defined([row[delta_field] for row in subset if row["assay"] == assay])
                            for assay in sorted({row["assay"] for row in subset})
                        ]
                        output.append({
                            "scope": scope,
                            "budget": budget,
                            "valley_status": stratum,
                            "policy_type": policy_name,
                            "horizon": horizon,
                            "unique_tasks": len(subset),
                            "assays": len({row["assay"] for row in subset}),
                            "micro_task_weighted_mean_performance": mean_defined(values),
                            "macro_equal_assay_mean_performance": mean_defined(assay_means),
                            "within_scope_task_bootstrap_ci_low": ci_low,
                            "within_scope_task_bootstrap_ci_high": ci_high,
                            "micro_H_vs_H1_mean_delta": mean_defined(delta_values),
                            "micro_H_vs_H1_paired_ci_low": delta_ci_low,
                            "micro_H_vs_H1_paired_ci_high": delta_ci_high,
                            "macro_H_vs_H1_mean_delta": mean_defined(assay_delta_means),
                            "macro_H_vs_H1_bootstrap_ci_low": macro_ci_low,
                            "macro_H_vs_H1_bootstrap_ci_high": macro_ci_high,
                            "mean_actual_budget": mean_defined([row[f"{prefix}_H{horizon}_actual_budget"] for row in subset]),
                            "top1_action_accuracy": mean_defined([row[f"H{horizon}_top1_action_accuracy"] for row in subset]) if prefix == "SURROGATE" and horizon in (1, 3) else None,
                            "top3_action_recall": mean_defined([row[f"H{horizon}_top3_action_recall"] for row in subset]) if prefix == "SURROGATE" and horizon in (1, 3) else None,
                            "mean_rank_true_best": mean_defined([row[f"H{horizon}_mean_rank_of_true_best"] for row in subset]) if prefix == "SURROGATE" and horizon in (1, 3) else None,
                            "mean_local_action_regret": mean_defined([row[f"H{horizon}_mean_local_action_regret"] for row in subset]) if prefix == "SURROGATE" and horizon in (1, 3) else None,
                            "mean_local_spearman": mean_defined([row[f"H{horizon}_mean_local_spearman"] for row in subset]) if prefix == "SURROGATE" and horizon in (1, 3) else None,
                            "cross_assay_inference": "EXPLORATORY_N4" if scope == "ALL_ASSAYS" else "WITHIN_ASSAY",
                        })
    return output


def gcn4_saturation(task_rows):
    tasks = collapse_tasks(task_rows)
    gcn = [row for row in tasks if row["assay"].startswith("GCN4")]
    ceiling = {
        (row["assay"], row["task_id"]): row["budgeted_oracle_fitness"]
        for row in gcn if row["budget"] == 100
    }
    earliest = {}
    for row in sorted(gcn, key=lambda item: item["budget"]):
        key = (row["assay"], row["task_id"])
        target = ceiling.get(key)
        if target is not None and math.isclose(row["budgeted_oracle_fitness"], target, rel_tol=0.0, abs_tol=1e-10):
            earliest.setdefault(key, row["budget"])
    rows = []
    for budget in BUDGETS:
        subset = [row for row in gcn if row["budget"] == budget]
        saturated = [earliest.get((row["assay"], row["task_id"]), math.inf) <= budget for row in subset]
        rows.append({
            "budget": budget,
            "tasks": len(subset),
            "fraction_reachable_ceiling_saturated": mean_defined([float(value) for value in saturated]),
            "median_earliest_saturation_budget": float(np.median([
                earliest[(row["assay"], row["task_id"])] for row in subset
                if (row["assay"], row["task_id"]) in earliest
            ])) if subset else None,
            "median_component_size": float(np.median([row["reachable_component_size"] for row in subset])) if subset else None,
            "median_unique_legal_lineage_paths": float(np.median([row["unique_legal_lineage_paths"] for row in subset if row["unique_legal_lineage_paths"] is not None])) if subset else None,
            "path_count_exact_fraction": mean_defined([
                float(row["unique_legal_lineage_paths_exact"])
                for row in subset if row["unique_legal_lineage_paths_exact"] is not None
            ]),
        })
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trajectories", type=Path, required=True)
    parser.add_argument("--task-metadata", type=Path, required=True)
    parser.add_argument("--phase06", type=Path, required=True)
    parser.add_argument("--task-out", type=Path, required=True)
    parser.add_argument("--assay-out", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument("--bootstrap-iterations", type=int, default=10000)
    args = parser.parse_args()
    trajectories = read_csv(args.trajectories)
    metadata = read_csv(args.task_metadata)
    frozen = read_csv(args.phase06)
    task_rows = build_task_rows(trajectories, metadata, frozen)
    assay_rows = build_assay_summary(task_rows, args.bootstrap_iterations)
    write_csv(args.task_out, task_rows)
    write_csv(args.assay_out, assay_rows)
    payload = {
        "phase": "06.5",
        "comparison": "H3 - H1; positive performance delta favors H3",
        "task_rows": len(task_rows),
        "unique_tasks": len({(row['assay'], row['task_id']) for row in task_rows}),
        "paired_groups": len({(row['assay'], row['task_id'], row['seed'], row['budget']) for row in task_rows}),
        "bootstrap_iterations": args.bootstrap_iterations,
        "bootstrap_unit": "unique task after averaging policy seeds",
        "micro": "task-weighted",
        "macro": "equal assay weight",
        "cross_assay_inference": "EXPLORATORY_N4",
        "opportunity_labels": Counter(row["planning_opportunity_label"] for row in task_rows),
        "valley_status": Counter(row["valley_status"] for row in task_rows),
        "gcn4_saturation": gcn4_saturation(task_rows),
    }
    payload["opportunity_labels"] = dict(payload["opportunity_labels"])
    payload["valley_status"] = dict(payload["valley_status"])
    args.json_out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
