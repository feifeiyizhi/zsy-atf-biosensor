#!/usr/bin/env python3
"""Build Phase 05 assay- and task-level landscape diagnostics."""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict, deque
from pathlib import Path
from statistics import mean, median

import numpy as np

from run_oracle_benchmark import load_graph


def percentile(values, q):
    return float(np.quantile(np.asarray(values, dtype=float), q)) if values else None


def component_map(nodes, adj):
    component_by_node, sizes = {}, []
    for start in nodes:
        if start in component_by_node:
            continue
        component_id = len(sizes)
        stack, members = [start], []
        component_by_node[start] = component_id
        while stack:
            cur = stack.pop(); members.append(cur)
            for nxt in adj[cur]:
                if nxt not in component_by_node:
                    component_by_node[nxt] = component_id; stack.append(nxt)
        sizes.append(len(members))
    return component_by_node, sizes


def shortest_path_length(start, target, adj, horizon):
    queue = deque([(start, 0)]); seen = {start}
    while queue:
        cur, distance = queue.popleft()
        if cur == target:
            return distance
        if distance >= horizon:
            continue
        for nxt in adj[cur]:
            if nxt not in seen:
                seen.add(nxt); queue.append((nxt, distance + 1))
    return None


def target_paths(start, target, adj, horizon):
    paths = []
    stack = [(start, [start])]
    while stack:
        cur, path = stack.pop()
        if cur == target:
            paths.append(path); continue
        if len(path) - 1 >= horizon:
            continue
        for nxt in sorted(adj[cur], reverse=True):
            if nxt not in path:
                stack.append((nxt, path + [nxt]))
    return paths


def valley_geometry(paths, nodes):
    candidates = []
    for path in paths:
        deltas = [nodes[b]["fitness"] - nodes[a]["fitness"] for a, b in zip(path, path[1:])]
        candidates.append((max([max(0.0, -delta) for delta in deltas], default=0.0), sum(delta < 0 for delta in deltas)))
    if not candidates:
        return None, None, "oracle_target_not_reachable_within_horizon"
    depth = min(item[0] for item in candidates)
    width = min(item[1] for item in candidates if abs(item[0] - depth) <= 1e-12)
    return depth, width, ""


def improving_access_fraction(nodes, adj, horizon):
    eligible = [node for node in nodes if adj[node]]
    successes = 0
    for start in eligible:
        start_fitness = nodes[start]["fitness"]
        queue = deque([(start, 0)]); seen = {start}; found = False
        while queue and not found:
            cur, distance = queue.popleft()
            if nodes[cur]["fitness"] > start_fitness:
                found = True; break
            if distance >= horizon:
                continue
            for nxt in adj[cur]:
                if nxt not in seen and nodes[nxt]["fitness"] >= nodes[cur]["fitness"]:
                    seen.add(nxt); queue.append((nxt, distance + 1))
        successes += found
    return successes / len(eligible) if eligible else None


def quality_class(nodes_count, edges_count, lcc_fraction, isolated_fraction, eligible_starts):
    if edges_count == 0 or eligible_starts < 20:
        return "UNSUITABLE_FOR_LOCAL_WALK", "no edges or fewer than 20 eligible connected starts"
    if lcc_fraction >= 0.90 and isolated_fraction <= 0.10:
        return "WELL_CONNECTED", "LCC fraction >=0.90 and isolated fraction <=0.10"
    if lcc_fraction < 0.10 or isolated_fraction >= 0.90:
        return "SPARSE_NEGATIVE_CONTROL", "has edges but LCC fraction <0.10 or isolated fraction >=0.90"
    if lcc_fraction >= 0.25 and isolated_fraction <= 0.60:
        return "MODERATELY_CONNECTED", "LCC fraction >=0.25 and isolated fraction <=0.60"
    return "UNSUITABLE_FOR_LOCAL_WALK", "falls between supported connectivity regimes"


def task_records(assay_dir, nodes, adj, component_by_node, component_sizes, horizon):
    records = []
    for seed in (1, 2, 3):
        path = assay_dir / f"oracle_compare_seed{seed}.csv"
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        for index, row in enumerate(rows):
            start, target = int(row["start_state"]), int(row["oracle_target"])
            paths = target_paths(start, target, adj, horizon)
            depth, width, geometry_reason = valley_geometry(paths, nodes)
            valley = row["valley_required"].lower() == "true"
            if not valley:
                width = None
                width_reason = "not_valley_required"
                depth = 0.0 if depth is not None else None
            else:
                width_reason = geometry_reason
            oracle = float(row["oracle_terminal"]); start_fitness = nodes[start]["fitness"]
            epsilon = max(1e-12, 1e-9 * max(abs(oracle - start_fitness), 1.0))
            denominator = oracle - start_fitness
            normalized_reason = "" if denominator > epsilon else "oracle_does_not_improve_start_beyond_epsilon"
            greedy_regret, lookahead_regret = float(row["greedy_regret"]), float(row["lookahead_regret"])
            records.append({
                "assay": assay_dir.name,
                "seed": seed,
                "task_id": f"{assay_dir.name}:s{seed}:t{index:03d}",
                "start_node": start,
                "start_fitness": start_fitness,
                "best_reachable_fitness": oracle,
                "valley_required": valley,
                "minimum_required_drop": depth if valley else 0.0,
                "minimum_required_drop_na_reason": geometry_reason,
                "valley_width": width,
                "valley_width_na_reason": width_reason,
                "shortest_oracle_path_length": shortest_path_length(start, target, adj, horizon),
                "greedy_terminal_fitness": float(row["greedy_terminal"]),
                "lookahead_terminal_fitness": float(row["lookahead_terminal"]),
                "oracle_terminal_fitness": oracle,
                "greedy_regret": greedy_regret,
                "lookahead_regret": lookahead_regret,
                "normalized_greedy_regret": greedy_regret / denominator if denominator > epsilon else None,
                "normalized_lookahead_regret": lookahead_regret / denominator if denominator > epsilon else None,
                "normalized_regret_na_reason": normalized_reason,
                "planning_gain": float(row["lookahead_terminal"]) - float(row["greedy_terminal"]),
                "component_size": component_sizes[component_by_node[start]],
                "start_degree": len(adj[start]),
            })
    return records


def assay_diagnostic(assay_dir, nodes, adj, component_sizes, tasks, horizon):
    node_count = len(nodes); edge_count = sum(len(value) for value in adj.values()) // 2
    degrees = [len(adj[node]) for node in nodes]
    nonisolated = [node for node in nodes if adj[node]]
    local_maxima = sum(nodes[node]["fitness"] >= max(nodes[n]["fitness"] for n in adj[node]) for node in nonisolated)
    local_minima = sum(nodes[node]["fitness"] <= min(nodes[n]["fitness"] for n in adj[node]) for node in nonisolated)
    improving = sum(any(nodes[n]["fitness"] > nodes[node]["fitness"] for n in adj[node]) for node in nodes)
    fitness = [nodes[node]["fitness"] for node in nodes]
    edge_x, edge_y = [], []
    for a in nodes:
        for b in adj[a]:
            if a < b:
                edge_x.extend([nodes[a]["fitness"], nodes[b]["fitness"]])
                edge_y.extend([nodes[b]["fitness"], nodes[a]["fitness"]])
    corr_reason = ""
    if len(edge_x) < 2 or np.std(edge_x) == 0 or np.std(edge_y) == 0:
        neighbor_corr = None; corr_reason = "fewer than two variable endpoint pairs"
    else:
        neighbor_corr = float(np.corrcoef(edge_x, edge_y)[0, 1])
    lcc = max(component_sizes, default=0); isolated = sum(degree == 0 for degree in degrees)
    graph_class, class_reason = quality_class(node_count, edge_count, lcc / node_count if node_count else 0, isolated / node_count if node_count else 0, len(nonisolated))
    valley_tasks = [task for task in tasks if task["valley_required"]]
    depths = [task["minimum_required_drop"] for task in valley_tasks if task["minimum_required_drop"] is not None]
    widths = [task["valley_width"] for task in valley_tasks if task["valley_width"] is not None]
    return {
        "assay": assay_dir.name,
        "nodes": node_count,
        "edges": edge_count,
        "edge_density": 2 * edge_count / (node_count * (node_count - 1)) if node_count > 1 else None,
        "edge_density_na_reason": "" if node_count > 1 else "fewer_than_two_nodes",
        "largest_connected_component_size": lcc,
        "largest_connected_component_fraction": lcc / node_count if node_count else None,
        "number_of_connected_components": len(component_sizes),
        "maximum_observed_mutation_depth": max(node["mutation_count"] for node in nodes.values()) if nodes else None,
        "median_degree": median(degrees) if degrees else None,
        "mean_degree": mean(degrees) if degrees else None,
        "degree_q25": percentile(degrees, 0.25),
        "degree_q75": percentile(degrees, 0.75),
        "maximum_degree": max(degrees, default=None),
        "fraction_isolated_nodes": isolated / node_count if node_count else None,
        "fraction_nodes_with_improving_neighbors": improving / node_count if node_count else None,
        "local_maxima_count_nonisolated": local_maxima,
        "local_maxima_density_all_nodes": local_maxima / node_count if node_count else None,
        "local_minima_count_nonisolated": local_minima,
        "local_minima_density_all_nodes": local_minima / node_count if node_count else None,
        "neighbor_fitness_correlation": neighbor_corr,
        "neighbor_fitness_correlation_na_reason": corr_reason,
        "greedy_failure_rate": sum(task["greedy_regret"] > 1e-12 for task in tasks) / len(tasks) if tasks else None,
        "valley_required_task_fraction": len(valley_tasks) / len(tasks) if tasks else None,
        "mean_valley_depth": mean(depths) if depths else None,
        "median_valley_depth": median(depths) if depths else None,
        "valley_depth_na_reason": "" if depths else "no_valley_required_tasks_with_supported_paths",
        "mean_valley_width": mean(widths) if widths else None,
        "median_valley_width": median(widths) if widths else None,
        "valley_width_na_reason": "" if widths else "no_valley_required_tasks_with_supported_paths",
        "accessible_improving_path_fraction_h3": improving_access_fraction(nodes, adj, horizon),
        "accessible_improving_path_fraction_na_reason": "" if nonisolated else "no_nonisolated_nodes",
        "fitness_variance_population": float(np.var(fitness)) if fitness else None,
        "fitness_iqr": percentile(fitness, 0.75) - percentile(fitness, 0.25) if fitness else None,
        "benchmark_quality_class": graph_class,
        "benchmark_quality_reason": class_reason,
        "phase06_suitable": graph_class in {"WELL_CONNECTED", "MODERATELY_CONNECTED"},
        "phase06_suitability_reason": "connected enough for local closed-pool walks" if graph_class in {"WELL_CONNECTED", "MODERATELY_CONNECTED"} else "retain as sparse control; not a primary Phase 06 benchmark",
    }


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmarks", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--horizon", type=int, default=3)
    args = parser.parse_args()
    all_tasks, diagnostics = [], []
    assay_dirs = sorted(path for path in args.benchmarks.iterdir() if list(path.glob("oracle_compare_seed[123].csv")))
    for assay_dir in assay_dirs:
        nodes, _, adj = load_graph(assay_dir / "nodes.csv", assay_dir / "edges.csv")
        component_by_node, component_sizes = component_map(nodes, adj)
        tasks = task_records(assay_dir, nodes, adj, component_by_node, component_sizes, args.horizon)
        all_tasks.extend(tasks)
        diagnostics.append(assay_diagnostic(assay_dir, nodes, adj, component_sizes, tasks, args.horizon))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "task_difficulty.csv", all_tasks)
    write_csv(args.out_dir / "assay_landscape_diagnostics.csv", diagnostics)
    metadata = {
        "status": "VERIFIED",
        "scope": "descriptive landscape/task characterization; assay-level relationships are EXPLORATORY — N=4 ASSAYS",
        "definitions": {
            "edge_density": "2E/(N(N-1)) for a simple undirected graph",
            "neighbor_fitness_correlation": "Pearson correlation over both orientations of each undirected edge",
            "greedy_failure_rate": "fraction of paired tasks with greedy_regret > 1e-12",
            "minimum_required_drop": "minimum across bounded paths to oracle target of the largest single-edge fitness decrease",
            "valley_width": "minimum downhill-edge count among paths attaining minimum_required_drop; NA for non-valley tasks",
            "accessible_improving_path_fraction_h3": "fraction of non-isolated nodes reaching a strictly fitter node within 3 non-decreasing edges",
        },
        "classification_rules": {
            "WELL_CONNECTED": "LCC fraction >=0.90 and isolated fraction <=0.10",
            "MODERATELY_CONNECTED": "LCC fraction >=0.25 and isolated fraction <=0.60",
            "SPARSE_NEGATIVE_CONTROL": "has edges and (LCC fraction <0.10 or isolated fraction >=0.90)",
            "UNSUITABLE_FOR_LOCAL_WALK": "no edges, fewer than 20 eligible starts, or outside supported connectivity regimes",
        },
        "assays": diagnostics,
    }
    (args.out_dir / "phase05_diagnostics_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"status": "VERIFIED", "assays": len(diagnostics), "tasks": len(all_tasks), "classes": Counter(row["benchmark_quality_class"] for row in diagnostics)}, default=dict, indent=2))


if __name__ == "__main__":
    main()
