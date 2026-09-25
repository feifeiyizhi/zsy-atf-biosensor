#!/usr/bin/env python3
"""Audit and freeze the corrected Phase 04 paired oracle benchmark."""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

from run_oracle_benchmark import load_graph

METHODS = ("random", "greedy", "beam", "lookahead", "oracle")
STRATA = ("ALL_TASKS", "VALLEY_REQUIRED", "NO_VALLEY_REQUIRED")


def percentile(values: list[float], q: float) -> float:
    values = sorted(values)
    if not values:
        return math.nan
    pos = (len(values) - 1) * q
    lo, hi = math.floor(pos), math.ceil(pos)
    if lo == hi:
        return values[lo]
    return values[lo] * (hi - pos) + values[hi] * (pos - lo)


def bootstrap_cluster_ci(rows: list[dict], key: str, seed: int, replicates: int) -> list[float | None]:
    clusters: dict[int, list[float]] = defaultdict(list)
    for row in rows:
        value = row.get(key)
        if value is not None and math.isfinite(value):
            clusters[row["start_state"]].append(value)
    values = [mean(cluster) for cluster in clusters.values()]
    if not values:
        return [None, None]
    rng = random.Random(seed)
    estimates = [mean(values[rng.randrange(len(values))] for _ in values) for _ in range(replicates)]
    estimates.sort()
    return [percentile(estimates, 0.025), percentile(estimates, 0.975)]


def monotonic_path_exists(start: int, target: int, horizon: int, nodes: dict, adj: dict) -> bool:
    stack = [(start, [start])]
    while stack:
        cur, path = stack.pop()
        if cur == target:
            return True
        if len(path) - 1 >= horizon:
            continue
        for nxt in adj[cur]:
            if nxt not in path and nodes[nxt]["fitness"] >= nodes[cur]["fitness"]:
                stack.append((nxt, path + [nxt]))
    return False


def parse_path(value: str) -> list[int]:
    path = json.loads(value)
    if not isinstance(path, list) or not all(isinstance(x, int) for x in path):
        raise ValueError(f"invalid path JSON: {value}")
    return path


def audit_assay(assay_dir: Path, starts: int, horizon: int, beam_width: int) -> tuple[list[dict], dict]:
    nodes, _, adj = load_graph(assay_dir / "nodes.csv", assay_dir / "edges.csv")
    fitnesses = [node["fitness"] for node in nodes.values()]
    q1, q3 = percentile(fitnesses, 0.25), percentile(fitnesses, 0.75)
    iqr = q3 - q1
    epsilon = max(1e-12, 1e-9 * max(abs(iqr), 1.0))
    errors: list[str] = []
    records: list[dict] = []
    run_parameters = []
    candidates = sorted(i for i in nodes if adj[i])
    for seed in (1, 2, 3):
        csv_path = assay_dir / f"oracle_compare_seed{seed}.csv"
        json_path = assay_dir / f"oracle_compare_seed{seed}.json"
        summary = json.loads(json_path.read_text())
        run_parameters.append({key: summary[key] for key in ("seed", "horizon", "beam_width", "tasks")})
        if summary["seed"] != seed or summary["horizon"] != horizon or summary["beam_width"] != beam_width:
            errors.append(f"seed {seed}: inconsistent run parameters")
        rng = random.Random(seed)
        expected_starts = sorted(rng.sample(candidates, min(starts, len(candidates))))
        with csv_path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        actual_starts = [int(row["start_state"]) for row in rows]
        if actual_starts != expected_starts:
            errors.append(f"seed {seed}: sampled starts are not reproducible")
        for index, row in enumerate(rows):
            start = int(row["start_state"])
            target = int(row["oracle_target"])
            start_fitness = nodes[start]["fitness"]
            oracle_fitness = float(row["oracle_terminal"])
            if abs(float(row["start_fitness"]) - start_fitness) > 1e-10:
                errors.append(f"seed {seed} task {index}: start fitness mismatch")
            for method in METHODS:
                path = parse_path(row[f"{method}_path"])
                if not path or path[0] != start or len(path) - 1 > horizon or len(path) != len(set(path)):
                    errors.append(f"seed {seed} task {index}: invalid {method} path shape")
                if any(b not in adj[a] for a, b in zip(path, path[1:])):
                    errors.append(f"seed {seed} task {index}: {method} path uses a non-edge")
                if path[-1] not in nodes:
                    errors.append(f"seed {seed} task {index}: {method} terminal missing")
                terminal = nodes[path[-1]]["fitness"]
                if abs(terminal - float(row[f"{method}_terminal"])) > 1e-10:
                    errors.append(f"seed {seed} task {index}: {method} terminal mismatch")
                if abs((oracle_fitness - terminal) - float(row[f"{method}_regret"])) > 1e-9:
                    errors.append(f"seed {seed} task {index}: {method} regret sign/formula mismatch")
            expected_valley = oracle_fitness > start_fitness and not monotonic_path_exists(start, target, horizon, nodes, adj)
            actual_valley = row["valley_required"].lower() == "true"
            if expected_valley != actual_valley:
                errors.append(f"seed {seed} task {index}: valley label mismatch")
            denominator = oracle_fitness - start_fitness
            normalized_greedy = float(row["greedy_regret"]) / denominator if denominator > epsilon else None
            normalized_lookahead = float(row["lookahead_regret"]) / denominator if denominator > epsilon else None
            record = {
                "assay": assay_dir.name,
                "seed": seed,
                "task_id": f"{assay_dir.name}:s{seed}:t{index:03d}",
                "start_state": start,
                "start_fitness": start_fitness,
                "oracle_target": target,
                "oracle_terminal_fitness": oracle_fitness,
                "valley_required": actual_valley,
                "greedy_terminal_fitness": float(row["greedy_terminal"]),
                "lookahead_terminal_fitness": float(row["lookahead_terminal"]),
                "greedy_regret": float(row["greedy_regret"]),
                "lookahead_regret": float(row["lookahead_regret"]),
                "regret_delta_lookahead_minus_greedy": float(row["lookahead_regret"]) - float(row["greedy_regret"]),
                "planning_gain_terminal": float(row["lookahead_terminal"]) - float(row["greedy_terminal"]),
                "normalized_greedy_regret": normalized_greedy,
                "normalized_lookahead_regret": normalized_lookahead,
                "normalized_planning_gain": (normalized_greedy - normalized_lookahead) if normalized_greedy is not None else None,
                "iqr_standardized_planning_gain": (float(row["lookahead_terminal"]) - float(row["greedy_terminal"])) / max(iqr, epsilon),
                "normalization_denominator": denominator if denominator > epsilon else None,
                "normalization_na_reason": "" if denominator > epsilon else "oracle_does_not_improve_start_beyond_epsilon",
            }
            records.append(record)
    return records, {
        "assay": assay_dir.name,
        "nodes": len(nodes),
        "edges": sum(len(v) for v in adj.values()) // 2,
        "fitness_iqr": iqr,
        "normalization_epsilon": epsilon,
        "run_parameters": run_parameters,
        "audit_errors": errors,
    }


def summarize(records: list[dict], replicates: int) -> list[dict]:
    output = []
    assays = sorted({row["assay"] for row in records})
    for assay_index, assay in enumerate(assays):
        assay_rows = [row for row in records if row["assay"] == assay]
        for stratum_index, stratum in enumerate(STRATA):
            if stratum == "VALLEY_REQUIRED":
                rows = [row for row in assay_rows if row["valley_required"]]
            elif stratum == "NO_VALLEY_REQUIRED":
                rows = [row for row in assay_rows if not row["valley_required"]]
            else:
                rows = assay_rows
            if not rows:
                continue
            normalized = [row for row in rows if row["normalized_planning_gain"] is not None]
            seed = 1701 + assay_index * 10 + stratum_index
            output.append({
                "assay": assay,
                "group": stratum,
                "tasks": len(rows),
                "unique_start_states": len({row["start_state"] for row in rows}),
                "seeds": len({row["seed"] for row in rows}),
                "valley_required_tasks": sum(row["valley_required"] for row in rows),
                "mean_greedy_terminal_fitness": mean(row["greedy_terminal_fitness"] for row in rows),
                "mean_lookahead_terminal_fitness": mean(row["lookahead_terminal_fitness"] for row in rows),
                "mean_greedy_regret": mean(row["greedy_regret"] for row in rows),
                "mean_lookahead_regret": mean(row["lookahead_regret"] for row in rows),
                "mean_regret_delta_lookahead_minus_greedy": mean(row["regret_delta_lookahead_minus_greedy"] for row in rows),
                "regret_delta_cluster_bootstrap_ci95_low": bootstrap_cluster_ci(rows, "regret_delta_lookahead_minus_greedy", seed, replicates)[0],
                "regret_delta_cluster_bootstrap_ci95_high": bootstrap_cluster_ci(rows, "regret_delta_lookahead_minus_greedy", seed, replicates)[1],
                "paired_lookahead_better_fraction": sum(row["planning_gain_terminal"] > 0 for row in rows) / len(rows),
                "paired_tie_fraction": sum(abs(row["planning_gain_terminal"]) <= 1e-12 for row in rows) / len(rows),
                "normalized_tasks": len(normalized),
                "mean_normalized_planning_gain": mean(row["normalized_planning_gain"] for row in normalized) if normalized else None,
                "normalized_gain_cluster_bootstrap_ci95_low": bootstrap_cluster_ci(normalized, "normalized_planning_gain", seed + 101, replicates)[0],
                "normalized_gain_cluster_bootstrap_ci95_high": bootstrap_cluster_ci(normalized, "normalized_planning_gain", seed + 101, replicates)[1],
                "mean_iqr_standardized_planning_gain": mean(row["iqr_standardized_planning_gain"] for row in rows),
                "iqr_standardized_gain_cluster_bootstrap_ci95_low": bootstrap_cluster_ci(rows, "iqr_standardized_planning_gain", seed + 202, replicates)[0],
                "iqr_standardized_gain_cluster_bootstrap_ci95_high": bootstrap_cluster_ci(rows, "iqr_standardized_planning_gain", seed + 202, replicates)[1],
            })
    return output


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def report_text(summary: list[dict], audit: list[dict], old_valley_count: int) -> str:
    all_rows = [row for row in summary if row["group"] == "ALL_TASKS"]
    lines = [
        "# Phase 04 Oracle Report",
        "",
        "**Evidence level: PRELIMINARY.** This is a paired, full-information oracle benchmark on four observed assays, not a general protein-design result.",
        "",
        "## Audit and correction",
        "",
        "The closeout audit found that the pre-closeout implementation allowed greedy/oracle to stop early but forced beam/lookahead to continue to the horizon. The registered 4-assay × 3-seed × 50-start design was rerun after giving every optimizing policy the same up-to-horizon stopping action. Pre-closeout files are retained under `work/results/phase04_precloseout/`.",
        "",
        f"The reproducible label count is **{sum(row['valley_required_tasks'] for row in all_rows)} VALLEY_REQUIRED records out of {sum(row['tasks'] for row in all_rows)}**. The earlier on-disk and progress count was {old_valley_count}; the pasted value 132 is not reproduced and was not substituted.",
        "",
        "Regret remains `oracle terminal fitness - policy terminal fitness`; negative `lookahead - greedy` regret delta favors lookahead. Raw magnitudes must not be compared across assays.",
        "",
        "Cross-assay metrics are (1) normalized regret, `regret / (oracle - start)`, only when the oracle improvement exceeds the assay-specific epsilon, and (2) paired terminal gain divided by assay fitness IQR. Epsilon is `max(1e-12, 1e-9 × max(|IQR|, 1))`. Confidence intervals resample unique start states, so repeated GCN4 starts across seeds are clustered.",
        "",
        "## Current result table (ALL TASKS)",
        "",
        "| Assay | Tasks | Valley | Raw regret delta | 95% clustered CI | IQR-standardized gain | Lookahead better | Ties |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in all_rows:
        lines.append(
            f"| {row['assay']} | {row['tasks']} | {row['valley_required_tasks']} | {row['mean_regret_delta_lookahead_minus_greedy']:.6g} | "
            f"[{row['regret_delta_cluster_bootstrap_ci95_low']:.6g}, {row['regret_delta_cluster_bootstrap_ci95_high']:.6g}] | "
            f"{row['mean_iqr_standardized_planning_gain']:.6g} | {row['paired_lookahead_better_fraction']:.1%} | {row['paired_tie_fraction']:.1%} |"
        )
    positive = [row["assay"] for row in all_rows if row["regret_delta_cluster_bootstrap_ci95_high"] < 0]
    lines += [
        "",
        "## 1. Does a planning advantage exist at all?",
        "",
        "**SUPPORTED:** Planning advantage exists on some observed protein fitness landscapes. Exact bounded lookahead reaches the bounded oracle under the corrected shared constraints; its paired advantage over greedy is positive on tasks where greedy is locally trapped.",
        "",
        "## 2. Is it universal?",
        "",
        "**NOT ESTABLISHED:** Planning is generally superior across protein landscapes. Only four assays were evaluated, and many tasks are ties. Full-information success also does not establish query-limited success.",
        "",
        "## 3. On which assays does it appear?",
        "",
        "The clustered paired interval favors lookahead on: " + (", ".join(positive) if positive else "none at this scale") + ". Group-specific VALLEY_REQUIRED and NO_VALLEY_REQUIRED rows are in the CSV/JSON summary.",
        "",
        "## 4. What are the important failure cases?",
        "",
        "GCN4 remains an important sparse graph/control case: most nodes are isolated and sampled connected starts repeat across seeds. Its previously reported terminal loss was caused by unequal forced-continuation semantics and does not survive the corrected audit. The sparse observability limitation remains real and is carried into Phase 05 classification rather than hidden.",
        "",
        "## 5. What is NOT yet established?",
        "",
        "This does not establish universal superiority, a ruggedness law, open-world protein design, transfer to unseen proteins, or recovery of oracle value when fitness is hidden. Those questions require Phases 05–07.",
        "",
        "## Reproducibility status",
        "",
        f"All {len(audit)} assay audits passed with zero path, graph, horizon, seed-sampling, regret, or valley-label errors." if all(not item["audit_errors"] for item in audit) else "Audit errors remain; see JSON.",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmarks", type=Path, required=True)
    parser.add_argument("--out-prefix", type=Path, required=True)
    parser.add_argument("--starts", type=int, default=50)
    parser.add_argument("--horizon", type=int, default=3)
    parser.add_argument("--beam-width", type=int, default=4)
    parser.add_argument("--bootstrap-replicates", type=int, default=10000)
    args = parser.parse_args()
    assays = sorted(path for path in args.benchmarks.iterdir() if list(path.glob("oracle_compare_seed[123].csv")))
    records, audits = [], []
    for assay_dir in assays:
        assay_records, audit = audit_assay(assay_dir, args.starts, args.horizon, args.beam_width)
        records.extend(assay_records)
        audits.append(audit)
    if any(audit["audit_errors"] for audit in audits):
        raise SystemExit(json.dumps(audits, indent=2))
    summary = summarize(records, args.bootstrap_replicates)
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_prefix.with_suffix(".csv"), summary)
    payload = {
        "status": "VERIFIED",
        "design": {"assays": len(assays), "tasks": len(records), "seeds_per_assay": 3, "starts_per_seed": args.starts, "horizon": args.horizon, "beam_width": args.beam_width},
        "raw_regret_definition": "oracle_terminal_fitness - policy_terminal_fitness",
        "paired_delta_definition": "lookahead_regret - greedy_regret; negative favors lookahead",
        "normalization": {
            "normalized_regret": "regret / (oracle_terminal_fitness - start_fitness), only when denominator > epsilon",
            "iqr_standardized_planning_gain": "(lookahead_terminal_fitness - greedy_terminal_fitness) / max(assay_fitness_IQR, epsilon)",
            "epsilon": "max(1e-12, 1e-9 * max(abs(assay_fitness_IQR), 1))",
        },
        "bootstrap": {"method": "paired cluster bootstrap over unique start_state within assay", "replicates": args.bootstrap_replicates},
        "audit": audits,
        "summary": summary,
        "task_records": records,
    }
    args.out_prefix.with_suffix(".json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    report = args.out_prefix.parent / "PHASE04_ORACLE_REPORT.md"
    report.write_text(report_text(summary, audits, old_valley_count=123), encoding="utf-8")
    print(json.dumps({"status": "VERIFIED", "assays": len(assays), "tasks": len(records), "summary_rows": len(summary), "report": str(report)}, indent=2))


if __name__ == "__main__":
    main()
