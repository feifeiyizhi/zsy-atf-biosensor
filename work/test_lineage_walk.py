#!/usr/bin/env python3
"""Invariant tests for the strict Phase 06 LINEAGE_WALK protocol."""
from __future__ import annotations

import inspect
import unittest
from collections import defaultdict

from run_lineage_walk import (
    budgeted_walk_oracle,
    fit_ridge,
    legal_neighbors,
    lineage_calibration_nodes,
    lineage_valley_metrics,
    plan_first_action,
    run_deployable,
    validate_lineage,
)


class LineageWalkTests(unittest.TestCase):
    def setUp(self):
        self.nodes = {
            0: {"fitness": 0.0},
            1: {"fitness": 0.2},
            2: {"fitness": -0.2},
            3: {"fitness": 1.0},
            4: {"fitness": 0.4},
            5: {"fitness": 1.5},
            6: {"fitness": 0.1},
        }
        self.adj = defaultdict(set)
        for left, right in ((0, 1), (0, 2), (1, 3), (2, 4), (4, 5), (1, 6)):
            self.adj[left].add(right)
            self.adj[right].add(left)
        self.mutations = {
            0: frozenset(),
            1: frozenset({"A1B"}),
            2: frozenset({"A2B"}),
            3: frozenset({"A1B", "A3B"}),
            4: frozenset({"A2B", "A4B"}),
            5: frozenset({"A2B", "A4B", "A5B"}),
            6: frozenset({"A1B", "A6B"}),
        }
        self.task = {
            "task_id": "toy:s1:t000",
            "start_node": 0,
            "best_reachable_fitness": 1.5,
            "valley_required": True,
        }

    def run_method(self, method, budgets=(3,), calibration_count=1):
        return run_deployable(
            method, "toy", self.task, self.nodes, self.adj,
            self.mutations, list(budgets), calibration_count, 1.0, 7, "test",
        )

    def test_calibration_preserves_current_legal_neighbors(self):
        calibration = lineage_calibration_nodes(0, set(self.nodes), self.adj, 3)
        self.assertEqual(calibration[0], 0)
        self.assertTrue(set(calibration[1:]).isdisjoint(self.adj[0]))
        self.assertEqual(legal_neighbors(0, set(calibration), self.adj), [1, 2])

    def test_every_executed_pair_is_a_graph_edge(self):
        for method in ("RANDOM_WALK", "ADAPTIVE_GREEDY", "ADAPTIVE_LOOKAHEAD"):
            row = self.run_method(method)[0]
            path = __import__("json").loads(row["executed_path"])
            self.assertEqual(validate_lineage(path, self.adj), (True, ""))

    def test_policy_cannot_jump_to_an_observed_branch(self):
        self.assertEqual(legal_neighbors(1, {0, 1, 2}, self.adj), [3, 6])
        self.assertNotIn(4, legal_neighbors(1, {0, 1, 2}, self.adj))

    def test_deployable_functions_do_not_accept_hidden_fitness(self):
        self.assertNotIn("nodes", inspect.signature(plan_first_action).parameters)
        self.assertNotIn("nodes", inspect.signature(fit_ridge).parameters)

    def test_ridge_is_retrained_after_every_reveal(self):
        row = self.run_method("ADAPTIVE_GREEDY")[0]
        events = __import__("json").loads(row["prequery_events"])
        self.assertEqual(
            [event["training_observation_count"] for event in events],
            list(range(1, len(events) + 1)),
        )

    def test_greedy_uses_horizon_one(self):
        row = self.run_method("ADAPTIVE_GREEDY")[0]
        self.assertEqual(row["planning_horizon"], 1)
        events = __import__("json").loads(row["prequery_events"])
        self.assertTrue(all(len(event["planned_path"]) == 2 for event in events))

    def test_lookahead_uses_horizon_three_and_executes_first_action(self):
        row = self.run_method("ADAPTIVE_LOOKAHEAD")[0]
        self.assertEqual(row["planning_horizon"], 3)
        path = __import__("json").loads(row["executed_path"])
        events = __import__("json").loads(row["prequery_events"])
        for index, event in enumerate(events):
            self.assertEqual(event["selected_node"], event["planned_path"][1])
            self.assertEqual(path[index + 1], event["selected_node"])
            self.assertLessEqual(len(event["planned_path"]) - 1, 3)

    def test_budgeted_oracle_respects_step_budget(self):
        fitness, path = budgeted_walk_oracle(0, 2, self.adj, self.nodes)
        self.assertEqual(fitness, 1.0)
        self.assertEqual(path, [0, 1, 3])
        self.assertLessEqual(len(path) - 1, 2)
        fitness, path = budgeted_walk_oracle(0, 3, self.adj, self.nodes)
        self.assertEqual(fitness, 1.5)
        self.assertEqual(path, [0, 2, 4, 5])

    def test_valley_success_requires_executed_downhill_and_recovery(self):
        success = lineage_valley_metrics([0, 2, 4, 5], self.nodes, 1.5, True, 1e-12)
        no_drop = lineage_valley_metrics([0, 1, 3], self.nodes, 1.0, True, 1e-12)
        no_recovery = lineage_valley_metrics([0, 2, 4], self.nodes, 1.5, True, 1e-12)
        self.assertTrue(success["strict_valley_crossing_success"])
        self.assertEqual(success["first_downhill_query"], 1)
        self.assertEqual(success["recovery_query"], 3)
        self.assertFalse(no_drop["strict_valley_crossing_success"])
        self.assertFalse(no_recovery["strict_valley_crossing_success"])

    def test_query_zero_success_is_initially_solved(self):
        task = {**self.task, "start_node": 5, "best_reachable_fitness": 1.5}
        rows = run_deployable(
            "RANDOM_WALK", "toy", task, self.nodes, self.adj,
            self.mutations, [1], 1, 1.0, 7, "test",
        )
        self.assertTrue(rows[0]["initially_solved"])
        self.assertFalse(rows[0]["top5_reached_after_calibration"])
        self.assertIsNone(rows[0]["queries_to_top5"])
        self.assertEqual(
            rows[0]["queries_to_top5_na_reason"],
            "threshold_reached_during_calibration",
        )

    def test_paired_policies_share_calibration(self):
        rows = [self.run_method(method)[0] for method in (
            "RANDOM_WALK", "ADAPTIVE_GREEDY", "ADAPTIVE_LOOKAHEAD",
        )]
        self.assertEqual(len({row["initial_observation_nodes"] for row in rows}), 1)
        self.assertEqual(len({row["initial_best_fitness"] for row in rows}), 1)


if __name__ == "__main__":
    unittest.main()
