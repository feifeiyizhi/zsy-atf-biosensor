#!/usr/bin/env python3
"""Deterministic method-gate tests for Phase 06.5 diagnostics."""
from __future__ import annotations

import csv
import json
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

from run_phase06_5_diagnostics import (
    beam_first_action,
    count_simple_maximal_paths,
    exact_first_action,
    load_frozen_design,
    local_decision_event,
    path_to_target_metrics,
)


def graph(edges):
    adj = defaultdict(set)
    for left, right in edges:
        adj[left].add(right)
        adj[right].add(left)
    return adj


def nodes(values):
    return {node: {"fitness": fitness} for node, fitness in values.items()}


class Phase065Tests(unittest.TestCase):
    def test_true_h1_selects_best_legal_neighbor(self):
        adj = graph([(0, 1), (0, 2), (0, 3)])
        values = {0: 0.0, 1: 1.0, 2: 3.0, 3: 2.0}
        action, path, _, exact = exact_first_action(0, {0}, adj, values, 1, 100)
        self.assertTrue(exact)
        self.assertEqual(action, 2)
        self.assertEqual(path, [0, 2])

    def test_true_h3_uses_terminal_value_and_first_action_only(self):
        adj = graph([(0, 1), (0, 2), (1, 3), (3, 4), (2, 5), (5, 6)])
        values = {0: 0.0, 1: -1.0, 2: 1.0, 3: -2.0, 4: 10.0, 5: 2.0, 6: 3.0}
        action, planned, _, exact = exact_first_action(0, {0}, adj, values, 3, 100)
        self.assertTrue(exact)
        self.assertEqual(action, 1)
        self.assertEqual(planned, [0, 1, 3, 4])
        self.assertEqual(len([action]), 1)

    def test_true_policy_never_uses_observed_or_repeats(self):
        adj = graph([(0, 1), (1, 2), (2, 0), (2, 3)])
        values = {0: 0.0, 1: 1.0, 2: 2.0, 3: 5.0}
        action, planned, _, exact = exact_first_action(0, {0, 1}, adj, values, 3, 100)
        self.assertTrue(exact)
        self.assertEqual(action, 2)
        self.assertNotIn(1, planned)
        self.assertEqual(len(planned), len(set(planned)))

    def test_exact_search_reports_cap(self):
        adj = graph([(0, 1), (0, 2), (1, 3), (2, 4)])
        values = {node: float(node) for node in range(5)}
        _, _, expansions, exact = exact_first_action(0, {0}, adj, values, 2, 1)
        self.assertFalse(exact)
        self.assertGreater(expansions, 1)

    def test_beam_h1_matches_frozen_tie_break(self):
        adj = graph([(0, 1), (0, 2)])
        predictions = {1: 1.0, 2: 1.0}
        action, path = beam_first_action(0, {0}, adj, predictions, 1, 16)
        self.assertEqual(action, 1)
        self.assertEqual(path, [0, 1])

    def test_local_decision_metrics(self):
        adj_nodes = nodes({0: 0.0, 1: 3.0, 2: 2.0, 3: 1.0})
        event = local_decision_event(
            0, [1, 2, 3], 2, {1: 1.0, 2: 3.0, 3: 2.0},
            adj_nodes, [0, 2], 1,
        )
        self.assertEqual(event["true_best_neighbor"], 1)
        self.assertEqual(event["predicted_best_neighbor"], 2)
        self.assertEqual(event["rank_of_true_best_neighbor"], 3)
        self.assertFalse(event["top1_action_accuracy"])
        self.assertTrue(event["top3_action_recall"])
        self.assertEqual(event["local_action_regret"], 1.0)
        self.assertAlmostEqual(event["local_spearman"], -0.5)

    def test_local_spearman_na_for_two_neighbors(self):
        event = local_decision_event(
            0, [1, 2], 1, {1: 2.0, 2: 1.0},
            nodes({0: 0.0, 1: 1.0, 2: 2.0}), [0, 1], 1,
        )
        self.assertIsNone(event["local_spearman"])

    def test_actionable_valley_metrics(self):
        adj = graph([(0, 1), (1, 2), (0, 3), (3, 2)])
        row_nodes = nodes({0: 1.0, 1: 0.5, 2: 3.0, 3: 0.0})
        task = {"start": 0, "start_fitness": 1.0, "target_fitness": 3.0}
        result = path_to_target_metrics(task, row_nodes, adj)
        self.assertTrue(result["reachable_under_strict_lineage"])
        self.assertTrue(result["requires_downhill"])
        self.assertEqual(result["minimum_number_downhill_steps"], 1)
        self.assertEqual(result["minimum_path_length_to_superior_target"], 2)
        self.assertEqual(result["required_horizon"], 2)
        self.assertAlmostEqual(result["minimum_true_downhill_depth"], 0.5)

    def test_actionable_valley_uses_specific_target(self):
        adj = graph([(0, 1), (1, 2), (0, 3)])
        row_nodes = nodes({0: 1.0, 1: 0.5, 2: 3.0, 3: 3.0})
        task = {"start": 0, "start_fitness": 1.0, "target_fitness": 3.0}
        any_equal = path_to_target_metrics(task, row_nodes, adj)
        original = path_to_target_metrics(task, row_nodes, adj, target_node=2)
        self.assertFalse(any_equal["requires_downhill"])
        self.assertTrue(original["requires_downhill"])

    def test_non_valley_target_has_zero_downhill_path(self):
        adj = graph([(0, 1), (1, 2)])
        result = path_to_target_metrics(
            {"start": 0, "start_fitness": 1.0, "target_fitness": 3.0},
            nodes({0: 1.0, 1: 2.0, 2: 3.0}), adj,
        )
        self.assertFalse(result["requires_downhill"])
        self.assertEqual(result["minimum_number_downhill_steps"], 0)
        self.assertIsNone(result["required_horizon"])

    def test_simple_maximal_path_count(self):
        adj = graph([(0, 1), (0, 2), (1, 3), (2, 4)])
        count, exact = count_simple_maximal_paths(0, adj, 100)
        self.assertTrue(exact)
        self.assertEqual(count, 2)

    def test_frozen_design_requires_expected_task_count(self):
        fields = [
            "assay", "task_id", "start_state", "initial_observation_nodes",
            "valley_required", "target_fitness", "start_fitness",
            "candidate_pool_size", "method", "seed", "requested_budget",
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tiny.csv"
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerow({
                    "assay": "A", "task_id": "t", "start_state": 0,
                    "initial_observation_nodes": json.dumps([0]),
                    "valley_required": False, "target_fitness": 1,
                    "start_fitness": 0, "candidate_pool_size": 2,
                    "method": "ADAPTIVE_GREEDY", "seed": 1,
                    "requested_budget": 10,
                })
            with self.assertRaises(ValueError):
                load_frozen_design(path)


if __name__ == "__main__":
    unittest.main()
