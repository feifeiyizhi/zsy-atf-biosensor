#!/usr/bin/env python3
"""Aggregate paired oracle runs and join graph ruggedness summaries."""
from __future__ import annotations
import argparse, csv, json, random
from pathlib import Path
from statistics import mean


def read_rows(paths):
    rows = []
    for path in sorted(paths):
        with path.open(newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                row["source_file"] = str(path)
                rows.append(row)
    return rows


def f(row, key):
    return float(row[key])


def bootstrap_ci(values, seed=17, n=10000):
    if not values:
        return [None, None]
    rng = random.Random(seed)
    samples = []
    for _ in range(n):
        samples.append(mean(values[rng.randrange(len(values))] for _ in values))
    samples.sort()
    return [samples[int(0.025 * (n - 1))], samples[int(0.975 * (n - 1))]]


def summarize(rows, landscape_root):
    grouped = {}
    for row in rows:
        grouped.setdefault(row["source_file"].split("/oracle_compare_seed")[0], []).append(row)
    output = []
    for key, items in sorted(grouped.items()):
        assay = items[0]["source_file"].split("/benchmarks/")[-1].split("/oracle_compare")[0]
        diffs = [f(r, "lookahead_regret") - f(r, "greedy_regret") for r in items]
        valley = [r for r in items if r["valley_required"].lower() == "true"]
        valley_diffs = [f(r, "lookahead_regret") - f(r, "greedy_regret") for r in valley]
        summary_path = landscape_root / assay / "landscape_summary.json"
        graph = json.loads(summary_path.read_text()) if summary_path.exists() else {}
        output.append({
            "assay": assay,
            "runs": len({r["source_file"] for r in items}),
            "tasks": len(items),
            "valley_tasks": len(valley),
            "mean_regret_delta_lookahead_minus_greedy": mean(diffs),
            "bootstrap_ci95": bootstrap_ci(diffs),
            "paired_lookahead_better_fraction": sum(x < 0 for x in diffs) / len(diffs),
            "valley_mean_regret_delta": mean(valley_diffs) if valley_diffs else None,
            "valley_bootstrap_ci95": bootstrap_ci(valley_diffs, seed=31) if valley_diffs else [None, None],
            "graph_nodes": graph.get("N_nodes"),
            "graph_edges": graph.get("N_edges"),
            "largest_component_fraction": graph.get("largest_connected_component_fraction"),
            "max_mutation_depth": graph.get("max_mutation_depth"),
            "local_maxima": graph.get("local_maxima"),
            "improving_neighbor_fraction": graph.get("fraction_nodes_with_improving_neighbors"),
        })
    return output


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmarks", type=Path, required=True)
    ap.add_argument("--landscapes", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    rows = read_rows(args.benchmarks.glob("*/oracle_compare_seed*.csv"))
    result = summarize(rows, args.landscapes)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"method": "paired_task_bootstrap", "bootstrap_replicates": 10000, "assays": result}, indent=2) + "\n")
    with args.out.with_suffix(".csv").open("w", newline="", encoding="utf-8") as f:
        fields = list(result[0]) if result else ["assay"]
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n"); w.writeheader(); w.writerows(result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
