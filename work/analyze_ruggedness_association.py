#!/usr/bin/env python3
"""Compute explicitly exploratory N=4 landscape/planning associations."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def rank(values):
    order = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    index = 0
    while index < len(order):
        end = index
        while end + 1 < len(order) and values[order[end + 1]] == values[order[index]]:
            end += 1
        average_rank = (index + end + 2) / 2.0
        for position in range(index, end + 1):
            ranks[order[position]] = average_rank
        index = end + 1
    return ranks


def spearman(x, y):
    if len(x) < 3:
        return None
    rx, ry = rank(x), rank(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    numerator = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    denominator = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return numerator / denominator if denominator else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--diagnostics", type=Path, required=True)
    parser.add_argument("--phase04", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    with args.diagnostics.open(newline="", encoding="utf-8") as handle:
        diagnostics = {row["assay"]: row for row in csv.DictReader(handle)}
    phase04 = json.loads(args.phase04.read_text())
    summaries = {row["assay"]: row for row in phase04["summary"] if row["group"] == "ALL_TASKS"}
    assays = sorted(set(diagnostics) & set(summaries))
    fields = (
        "largest_connected_component_fraction",
        "greedy_failure_rate",
        "local_maxima_density_all_nodes",
        "valley_required_task_fraction",
        "fraction_isolated_nodes",
        "neighbor_fitness_correlation",
    )
    values = []
    for assay in assays:
        row = diagnostics[assay]
        values.append({
            "assay": assay,
            "benchmark_quality_class": row["benchmark_quality_class"],
            "normalized_planning_gain": float(summaries[assay]["mean_normalized_planning_gain"]),
            **{field: float(row[field]) for field in fields},
        })
    associations = {
        field: spearman([row[field] for row in values], [row["normalized_planning_gain"] for row in values])
        for field in fields
    }
    output = {
        "label": "EXPLORATORY — N=4 ASSAYS",
        "n_assays": len(values),
        "method": "descriptive assay-level Spearman rank correlation",
        "planning_gain": "mean task-level normalized planning gain from corrected Phase 04",
        "associations": associations,
        "assay_values": values,
        "interpretation": "No formal generalization is supported. These four points define hypotheses and benchmark suitability only.",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
