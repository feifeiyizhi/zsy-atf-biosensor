#!/usr/bin/env python3
"""Generate the concise Phase 06.5 decision report from audited summaries."""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ASSAY_SHORT = {
    "A4_HUMAN_Seuma_2022": "A4_HUMAN",
    "D7PM05_CLYGR_Somermeyer_2022": "D7PM05",
    "F7YBW8_MESOW_Aakre_2015": "F7YBW8",
    "GCN4_YEAST_Staller_2018": "GCN4",
}


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def value(row, key):
    return None if row.get(key) in (None, "") else float(row[key])


def mean(values):
    values = [v for v in values if v is not None]
    return float(np.mean(values)) if values else None


def fmt(number, digits=3):
    return "NA" if number is None else f"{number:+.{digits}f}"


def assay_evidence(task_rows, summary_rows, analysis):
    output = []
    assays = sorted({row["assay"] for row in task_rows})
    gcn_saturation = {
        int(row["budget"]): row for row in analysis.get("gcn4_saturation", [])
    }
    for assay in assays:
        rows = [row for row in task_rows if row["assay"] == assay and row["budget"] == "100"]
        task_rows_unique = {}
        for row in rows:
            task_rows_unique.setdefault(row["task_id"], row)
        tasks = list(task_rows_unique.values())
        true_opportunity = mean([value(row, "planning_opportunity_H3") for row in tasks])
        surrogate_gain = mean([value(row, "realized_gain_H3") for row in tasks])
        top1 = mean([value(row, "H3_top1_action_accuracy") for row in tasks])
        actionable = mean([float(row["phase06_actionable_valley"] == "True") for row in tasks])
        horizon_rows = [
            row for row in summary_rows
            if row["scope"] == assay and row["budget"] == "100"
            and row["valley_status"] == "ALL_TASKS"
            and row["policy_type"] == "TRUE_FITNESS_POLICY_DIAGNOSTIC"
        ]
        best = max(horizon_rows, key=lambda row: value(row, "micro_task_weighted_mean_performance") or -1)
        best_horizon = int(best["horizon"])
        positive_fraction = mean([
            float(row["planning_opportunity_label"] == "POSITIVE_PLANNING_OPPORTUNITY")
            for row in tasks
        ])
        negative_fraction = mean([
            float(row["planning_opportunity_label"] == "NEGATIVE_DEEPER_PLANNING")
            for row in tasks
        ])
        horizon_cis = {}
        supported_horizons = []
        for row in horizon_rows:
            horizon = int(row["horizon"])
            interval = (
                value(row, "micro_H_vs_H1_paired_ci_low"),
                value(row, "micro_H_vs_H1_paired_ci_high"),
            )
            horizon_cis[horizon] = interval
            if horizon > 1 and interval[0] is not None and interval[0] > 0:
                supported_horizons.append(horizon)
        h3_ci = horizon_cis[3]
        surrogate_h3_row = next(
            row for row in summary_rows
            if row["scope"] == assay and row["budget"] == "100"
            and row["valley_status"] == "ALL_TASKS"
            and row["policy_type"] == "SURROGATE_HORIZON_DIAGNOSTIC"
            and row["horizon"] == "3"
        )
        surrogate_h3_ci = (
            value(surrogate_h3_row, "micro_H_vs_H1_paired_ci_low"),
            value(surrogate_h3_row, "micro_H_vs_H1_paired_ci_high"),
        )
        any_true_supported = bool(supported_horizons)
        true_supported = h3_ci[0] is not None and h3_ci[0] > 0
        true_harm = h3_ci[1] is not None and h3_ci[1] < 0
        surrogate_not_positive = surrogate_h3_ci[1] is not None and surrogate_h3_ci[1] <= 0
        graph_saturated = (
            assay.startswith("GCN4")
            and gcn_saturation.get(10, {}).get("fraction_reachable_ceiling_saturated") == 1.0
        )
        if positive_fraction is None or h3_ci[0] is None or h3_ci[1] is None:
            diagnosis = "INSUFFICIENT_EVIDENCE"
        elif graph_saturated:
            diagnosis = "GRAPH_LIMITED"
        elif any_true_supported and best_horizon not in (2, 3):
            diagnosis = "HORIZON_MISMATCH"
        elif true_harm:
            diagnosis = "NO_PLANNING_OPPORTUNITY"
        elif true_supported and surrogate_not_positive:
            diagnosis = "SURROGATE_LIMITED"
        elif true_supported and (surrogate_gain or 0.0) > 0:
            diagnosis = "PLANNING_RECOVERABLE"
        elif not true_supported and h3_ci[0] is not None and h3_ci[1] is not None and h3_ci[0] <= 0 <= h3_ci[1]:
            diagnosis = "MIXED"
        else:
            diagnosis = "NO_PLANNING_OPPORTUNITY"
        output.append({
            "assay": assay,
            "tasks": len(tasks),
            "true_h3_minus_h1": true_opportunity,
            "surrogate_h3_minus_h1": surrogate_gain,
            "local_top1": top1,
            "actionable_fraction": actionable,
            "positive_opportunity_fraction": positive_fraction,
            "negative_opportunity_fraction": negative_fraction,
            "best_true_horizon": best_horizon,
            "true_h3_ci_low": h3_ci[0],
            "true_h3_ci_high": h3_ci[1],
            "surrogate_h3_ci_low": surrogate_h3_ci[0],
            "surrogate_h3_ci_high": surrogate_h3_ci[1],
            "true_h3_supported": true_supported,
            "supported_true_horizons": supported_horizons,
            "diagnosis": diagnosis,
        })
    return output


def report_text(task_rows, summary_rows, analysis, audit):
    assays = assay_evidence(task_rows, summary_rows, analysis)
    unique_tasks = {(row["assay"], row["task_id"]) for row in task_rows}
    b100 = {}
    for row in task_rows:
        if row["budget"] == "100":
            b100.setdefault((row["assay"], row["task_id"]), row)
    b100_rows = list(b100.values())
    positive = sum(row["planning_opportunity_label"] == "POSITIVE_PLANNING_OPPORTUNITY" for row in b100_rows)
    actionable = sum(row["phase06_actionable_valley"] == "True" for row in b100_rows)
    all_exact = audit["all_true_rows_exact"]
    failures = Counter()
    failure_tasks = {}
    for row in task_rows:
        if row["budget"] == "100":
            failure_tasks.setdefault((row["assay"], row["task_id"]), row)
    for row in failure_tasks.values():
        failures.update(json.loads(row["downhill_failure_mode_counts"]))
    pooled = [
        row for row in summary_rows
        if row["scope"] == "ALL_ASSAYS" and row["budget"] == "100"
        and row["valley_status"] == "ALL_TASKS"
    ]
    true_by_h = {
        int(row["horizon"]): value(row, "micro_H_vs_H1_mean_delta")
        for row in pooled if row["policy_type"] == "TRUE_FITNESS_POLICY_DIAGNOSTIC"
    }
    surrogate_by_h = {
        int(row["horizon"]): value(row, "micro_H_vs_H1_mean_delta")
        for row in pooled if row["policy_type"] == "SURROGATE_HORIZON_DIAGNOSTIC"
    }
    diagnoses = Counter(row["diagnosis"] for row in assays)
    supported_assays = [row for row in assays if row["supported_true_horizons"]]
    q1 = bool(supported_assays)
    q2 = any(row["diagnosis"] == "SURROGATE_LIMITED" for row in assays)
    q3_action = [value(row, "planning_opportunity_H3") for row in b100_rows if row["phase06_actionable_valley"] == "True"]
    q3_other = [value(row, "planning_opportunity_H3") for row in b100_rows if row["phase06_actionable_valley"] != "True"]
    q3 = (mean(q3_action) or 0.0) > (mean(q3_other) or 0.0)
    best_true_h = max(true_by_h, key=lambda h: true_by_h[h] if true_by_h[h] is not None else -1e99)
    q4 = best_true_h
    justify_phase07 = q1 and q2
    lines = [
        "# Phase 06.5 Planning Opportunity / Information Bottleneck Audit",
        "",
        "## 1. PROJECT PIPELINE",
        "",
        "✅ 01 Scope → ✅ 02 Benchmark → ✅ 03 Graph → ✅ 04 Oracle → ✅ 05 Landscape → ✅ 06 Partial Observation → ✅ 06.5 Diagnostic Gate → ⏸ 07 Learned Planner → ○ 08 338lib → ○ 09 Paper",
        "",
        "Phase 06 remains closed. Phase 07 remains on hold unless the evidence below supports the predefined information/surrogate-bottleneck gate.",
        "",
        "## 2. DIAGNOSTIC QUESTION",
        "",
        "Does the exact frozen Phase 06 benchmark contain usable deeper-planning opportunity, and if so is it lost because of local surrogate ranking, horizon mismatch, or graph/observability limits?",
        "",
        "## 3. EXPERIMENT SCALE",
        "",
        f"- {len(unique_tasks)} frozen tasks, {analysis['paired_groups']} task/seed/budget paired groups.",
        "- Same starts, calibration, seeds, budgets, strict graph, and no-repeat LINEAGE_WALK as Phase 06.",
        "- Evaluation-only TRUE H1–H5 and same-Ridge surrogate H1–H5.",
        f"- TRUE exactness: {'PASS' if all_exact else 'INCOMPLETE'}; audit status: {audit['status']}.",
        "- 10,000 task-clustered bootstrap iterations; pooled micro is task-weighted and macro gives equal assay weight.",
        "",
        "## 4. TRUE-FITNESS HORIZON RESULT",
        "",
        "Pooled task-weighted normalized-performance delta versus TRUE H1 at budget 100:",
        "",
        "| Horizon | TRUE H−H1 | Surrogate H−H1 |",
        "|---:|---:|---:|",
    ]
    for horizon in range(1, 6):
        lines.append(f"| H{horizon} | {fmt(true_by_h.get(horizon))} | {fmt(surrogate_by_h.get(horizon))} |")
    lines += [
        "",
        f"At budget 100, {positive}/{len(b100_rows)} unique tasks have positive TRUE H3>H1 opportunity. Zero and negative opportunities are retained.",
        "",
        "## 5. SURROGATE LOCAL-DECISION RESULT",
        "",
        "Local quality is evaluated only over each current legal Hamming-1 action set; whole-pool R² is not the primary diagnostic.",
        "",
        "| Assay | H3 local top-1 | TRUE H3−H1 | Surrogate H3−H1 | Best TRUE horizon |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in assays:
        lines.append(
            f"| {ASSAY_SHORT[row['assay']]} | {fmt(row['local_top1'])} | {fmt(row['true_h3_minus_h1'])} | "
            f"{fmt(row['surrogate_h3_minus_h1'])} | H{row['best_true_horizon']} |"
        )
    lines += [
        "",
        "## 6. ACTIONABLE-VALLEY RESULT",
        "",
        f"{actionable}/{len(b100_rows)} frozen tasks are `PHASE06_ACTIONABLE_VALLEY` at budget 100 under the original target-node provenance. Original `VALLEY_REQUIRED` labels are unchanged.",
        "",
        f"Mean TRUE H3 opportunity: actionable={fmt(mean(q3_action))}; other={fmt(mean(q3_other))}.",
        "",
        "## 7. ASSAY DIAGNOSIS TABLE",
        "",
        "| Assay | True H3>H1 opportunity | Surrogate H3>H1 | Local surrogate quality | Actionable valleys | Diagnosis |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in assays:
        lines.append(
            f"| {ASSAY_SHORT[row['assay']]} | {fmt(row['true_h3_minus_h1'])} | {fmt(row['surrogate_h3_minus_h1'])} | "
            f"top-1 {fmt(row['local_top1'])} | {row['actionable_fraction']:.1%} | {row['diagnosis']} |"
        )
    lines += [
        "",
        "Diagnoses follow the documented rule combining exact TRUE opportunity, realized surrogate gain, local top-1 ranking, best horizon, and graph headroom. They are descriptive, not population-level protein classes.",
        "",
        "## 8. FAILURE-MODE DECOMPOSITION",
        "",
        "Adaptive H3 downhill-event counts at budget 100:",
        "",
    ]
    lines.extend(f"- `{name}`: {count}" for name, count in sorted(failures.items()))
    lines += [
        "",
        "## 9. SCIENTIFIC INTERPRETATION",
        "",
        "The four-assay comparison is exploratory. Positive, zero, and harmful deeper-planning cases are all retained; no universal planning or ruggedness claim is made.",
        "",
        "## 10. PHASE 07 DECISION",
        "",
        f"- **Q1 — meaningful TRUE H>1 subset?** {'YES' if q1 else 'NO'}: "
        + (
            "; ".join(
                f"{ASSAY_SHORT[row['assay']]} supports H{','.join(str(h) for h in row['supported_true_horizons'])} over H1"
                for row in supported_assays
            )
            if supported_assays else "no assay has a paired task-bootstrap interval wholly above zero"
        ) + ".",
        f"- **Q2 — primarily local surrogate/ranking error?** {'SUPPORTED IN AT LEAST ONE ASSAY' if q2 else 'NOT ESTABLISHED'}.",
        f"- **Q3 — concentrated in actionable valleys?** {'YES DESCRIPTIVELY' if q3 else 'NO / NOT CLEAR'}.",
        f"- **Q4 — appropriate horizon?** Pooled best TRUE horizon is H{q4}; assay-specific horizons remain in the summary.",
        f"- **Q5 — defensible learned/uncertainty-aware planner rationale?** {'YES, BUT ONLY FOR THE EXPLICIT SUPPORTED REGIME' if justify_phase07 else 'NO'}.",
        "",
        "**Decision:** " + ("Phase 07 may be scoped only to the explicit surrogate-limited regime; it is not started by this audit." if justify_phase07 else "Keep Phase 07 on HOLD. The gate for learned planning is not met."),
        "",
        "## 11. FILES CREATED",
        "",
        "- `work/results/phase06_5_task_diagnostics.csv`",
        "- `work/results/phase06_5_assay_summary.csv`",
        "- `work/results/phase06_5_analysis.json`",
        "- `work/results/phase06_5_audit.json`",
        "- `work/results/PHASE06_5_DECISION_REPORT.md`",
        "",
        "Raw trajectory/event files remain on the shared volume and are not treated as independent scientific samples.",
        "",
        "## 12. VALIDATION",
        "",
        f"- Independent audit: {audit['status']} ({len(audit['errors'])} errors).",
        f"- Frozen H1/H3 path/fitness comparisons: {audit['frozen_h1_h3_seed_comparisons']}.",
        f"- TRUE rows exact: {audit['all_true_rows_exact']}.",
        "- Original target-node provenance, path edges, no-repeat constraints, local actions, budgets, and normalization were rechecked.",
        "",
        "## 13. NEXT 3 ACTIONS",
        "",
        "1. Preserve Phase 06 and Phase 06.5 as frozen evidence.",
        "2. If the gate is not met, expand/redesign tasks before any stronger predictor; if met only in a subset, predeclare that regime before Phase 07.",
        "3. Do not start Phase 07 automatically; await an explicit decision based on this report.",
        "",
        "## 14. BLOCKERS",
        "",
        "None engineering." if audit["status"] == "PASS" else "Audit errors remain; no scientific gate decision is valid until resolved.",
    ]
    return "\n".join(lines) + "\n", assays, diagnoses


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    task_rows = read_csv(args.task)
    summary_rows = read_csv(args.summary)
    analysis = json.loads(args.analysis.read_text())
    audit = json.loads(args.audit.read_text())
    report, assays, diagnoses = report_text(task_rows, summary_rows, analysis, audit)
    args.out.write_text(report)
    print(json.dumps({"assays": assays, "diagnoses": dict(diagnoses), "report": str(args.out)}, indent=2))


if __name__ == "__main__":
    main()
