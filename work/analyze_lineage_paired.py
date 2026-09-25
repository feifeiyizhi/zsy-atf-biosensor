#!/usr/bin/env python3
"""Paired task-clustered analysis for the Phase 06 LINEAGE_WALK benchmark."""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

GREEDY = "ADAPTIVE_GREEDY"
LOOKAHEAD = "ADAPTIVE_LOOKAHEAD"
STRATA = ("ALL_TASKS", "VALLEY_REQUIRED", "NO_VALLEY_REQUIRED")


def parse_optional_float(value):
    return None if value in (None, "") else float(value)


def parse_bool(value):
    if value in (True, "True", "true", "1", 1):
        return True
    if value in (False, "False", "false", "0", 0):
        return False
    return None


def read_pairs(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    by_key = defaultdict(dict)
    for row in rows:
        if row["method"] in {GREEDY, LOOKAHEAD}:
            key = (
                row["assay"], row["task_id"], int(row["seed"]),
                int(row["requested_budget"]),
            )
            if row["method"] in by_key[key]:
                raise ValueError(f"duplicate pair member: {key} {row['method']}")
            by_key[key][row["method"]] = row
    incomplete = [key for key, methods in by_key.items() if set(methods) != {GREEDY, LOOKAHEAD}]
    if incomplete:
        raise ValueError(f"incomplete pairs: {incomplete[:5]}")
    return [(key, methods[GREEDY], methods[LOOKAHEAD]) for key, methods in sorted(by_key.items())]


def in_stratum(greedy, stratum):
    valley = parse_bool(greedy["valley_required"])
    return (
        stratum == "ALL_TASKS"
        or (stratum == "VALLEY_REQUIRED" and valley)
        or (stratum == "NO_VALLEY_REQUIRED" and not valley)
    )


def pair_value(greedy, lookahead, metric):
    if metric == "normalized_regret_delta":
        left = parse_optional_float(greedy["normalized_regret"])
        right = parse_optional_float(lookahead["normalized_regret"])
    elif metric == "normalized_best_gain_auc_delta":
        left = parse_optional_float(greedy["normalized_best_gain_auc"])
        right = parse_optional_float(lookahead["normalized_best_gain_auc"])
    elif metric == "strict_valley_success_delta":
        left = parse_bool(greedy["strict_valley_crossing_success"])
        right = parse_bool(lookahead["strict_valley_crossing_success"])
        left = None if left is None else float(left)
        right = None if right is None else float(right)
    else:
        raise ValueError(metric)
    return None if left is None or right is None else right - left


def clustered_values(pairs, metric):
    by_task = defaultdict(list)
    raw = []
    for key, greedy, lookahead in pairs:
        value = pair_value(greedy, lookahead, metric)
        if value is not None:
            by_task[(key[0], key[1])].append(value)
            raw.append(value)
    cluster_means = np.asarray([np.mean(values) for values in by_task.values()], dtype=float)
    return np.asarray(raw, dtype=float), cluster_means


def bootstrap_mean(cluster_means, seed, iterations):
    if len(cluster_means) == 0:
        return None, None
    if len(cluster_means) == 1:
        value = float(cluster_means[0])
        return value, value
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(cluster_means), size=(iterations, len(cluster_means)))
    means = cluster_means[indices].mean(axis=1)
    low, high = np.quantile(means, [0.025, 0.975])
    return float(low), float(high)


def summarize_metric(pairs, metric, seed, iterations):
    raw, clusters = clustered_values(pairs, metric)
    low, high = bootstrap_mean(clusters, seed, iterations)
    return {
        "paired_seed_records": int(len(raw)),
        "unique_task_clusters": int(len(clusters)),
        "mean_delta": float(raw.mean()) if len(raw) else None,
        "median_delta": float(np.median(raw)) if len(raw) else None,
        "cluster_mean_delta": float(clusters.mean()) if len(clusters) else None,
        "bootstrap_ci_low": low,
        "bootstrap_ci_high": high,
        "lookahead_wins": int(np.sum(raw < 0)) if metric == "normalized_regret_delta" else int(np.sum(raw > 0)),
        "ties": int(np.sum(np.isclose(raw, 0.0, atol=1e-12))),
        "lookahead_losses": int(np.sum(raw > 0)) if metric == "normalized_regret_delta" else int(np.sum(raw < 0)),
    }


def build_summaries(pairs, iterations):
    assays = sorted({key[0] for key, _, _ in pairs})
    scopes = ["ALL_ASSAYS", *assays]
    metrics = (
        "normalized_regret_delta",
        "normalized_best_gain_auc_delta",
        "strict_valley_success_delta",
    )
    output = []
    for scope_index, scope in enumerate(scopes):
        for budget in sorted({key[3] for key, _, _ in pairs}):
            for stratum_index, stratum in enumerate(STRATA):
                subset = [
                    item for item in pairs
                    if item[0][3] == budget
                    and (scope == "ALL_ASSAYS" or item[0][0] == scope)
                    and in_stratum(item[1], stratum)
                ]
                for metric_index, metric in enumerate(metrics):
                    if metric == "strict_valley_success_delta" and stratum == "NO_VALLEY_REQUIRED":
                        continue
                    summary = summarize_metric(
                        subset, metric,
                        61000 + scope_index * 1000 + budget * 10 + stratum_index * 3 + metric_index,
                        iterations,
                    )
                    output.append({
                        "scope": scope,
                        "requested_budget": budget,
                        "stratum": stratum,
                        "metric": metric,
                        "delta_definition": "ADAPTIVE_LOOKAHEAD - ADAPTIVE_GREEDY",
                        "favorable_direction": "negative" if metric == "normalized_regret_delta" else "positive",
                        **summary,
                    })
    return output


def diagnostics(pairs):
    groups = defaultdict(list)
    for key, greedy, lookahead in pairs:
        for method, row in ((GREEDY, greedy), (LOOKAHEAD, lookahead)):
            rmse = parse_optional_float(row["prequery_RMSE"])
            regret = parse_optional_float(row["normalized_regret"])
            groups[(key[0], key[3], method)].append((rmse, regret, int(row["actual_budget"])))
    rows = []
    for (assay, budget, method), values in sorted(groups.items()):
        rmse = np.asarray([value[0] for value in values if value[0] is not None], dtype=float)
        paired = [(value[0], value[1]) for value in values if value[0] is not None and value[1] is not None]
        correlation = None
        if len(paired) >= 3:
            x = np.asarray([value[0] for value in paired])
            y = np.asarray([value[1] for value in paired])
            if len(set(x)) > 1 and len(set(y)) > 1:
                correlation = float(np.corrcoef(x, y)[0, 1])
        rows.append({
            "assay": assay,
            "requested_budget": budget,
            "method": method,
            "records": len(values),
            "mean_prequery_RMSE": float(rmse.mean()) if len(rmse) else None,
            "rmse_regret_pearson_exploratory": correlation,
            "mean_actual_budget": float(np.mean([value[2] for value in values])),
        })
    return rows


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--diagnostics", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--iterations", type=int, default=10000)
    args = parser.parse_args()
    pairs = read_pairs(args.input)
    summaries = build_summaries(pairs, args.iterations)
    diagnostic_rows = diagnostics(pairs)
    write_csv(args.summary, summaries)
    write_csv(args.diagnostics, diagnostic_rows)
    args.json.write_text(json.dumps({
        "comparison": "ADAPTIVE_LOOKAHEAD - ADAPTIVE_GREEDY",
        "pairing_key": ["assay", "task_id", "seed", "requested_budget"],
        "bootstrap_unit": "unique assay-task cluster after averaging policy seeds",
        "bootstrap_iterations": args.iterations,
        "paired_records": len(pairs),
        "summaries": summaries,
        "diagnostics": diagnostic_rows,
    }, indent=2) + "\n")
    print(json.dumps({
        "paired_records": len(pairs), "summary_rows": len(summaries),
        "diagnostic_rows": len(diagnostic_rows),
    }, indent=2))


if __name__ == "__main__":
    main()
