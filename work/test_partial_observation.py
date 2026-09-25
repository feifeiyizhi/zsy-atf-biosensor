#!/usr/bin/env python3
"""Deterministic tests for the Phase 06 closed-pool benchmark."""
from __future__ import annotations

import csv
import inspect
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

from run_partial_observation import (
    append_registry,
    calibration_nodes,
    normalized_auc_best_curve,
    read_tasks,
    ridge_predict,
    run_policy,
    select_deployable_query,
    snapshot,
    strict_downhill_path_exists,
)


class PartialObservationTests(unittest.TestCase):
    def setUp(self):
        self.nodes = {
            0: {"fitness": 0.0}, 1: {"fitness": 0.2}, 2: {"fitness": -0.1},
            3: {"fitness": 1.0}, 4: {"fitness": 0.4}, 5: {"fitness": 1.5},
        }
        self.adj = defaultdict(set)
        for a, b in ((0, 1), (0, 2), (1, 3), (2, 4), (4, 5)):
            self.adj[a].add(b); self.adj[b].add(a)
        self.mutations = {
            0: frozenset(), 1: frozenset({"A1B"}), 2: frozenset({"A2B"}),
            3: frozenset({"A1B", "A3B"}), 4: frozenset({"A2B", "A4B"}),
            5: frozenset({"A2B", "A4B", "A5B"}),
        }
        self.task = {"task_id": "toy:s1:t000", "start_node": 0, "best_reachable_fitness": 1.0, "valley_required": True}

    def test_deployable_selector_has_no_hidden_nodes_argument(self):
        self.assertNotIn("nodes", inspect.signature(select_deployable_query).parameters)

    def test_ridge_uses_only_observed_fitness(self):
        predictions, r2, rho = ridge_predict({0: 0.0, 1: 0.2}, [2, 3], self.mutations, 1.0)
        self.assertEqual(set(predictions), {2, 3})
        self.assertTrue(all(isinstance(value, float) for value in predictions.values()))
        self.assertIsNotNone(r2)

    def test_calibration_is_deterministic(self):
        expected = calibration_nodes(0, set(self.nodes), self.adj, 4)
        self.assertEqual(expected, calibration_nodes(0, set(self.nodes), self.adj, 4))
        self.assertEqual(expected[0], 0)

    def test_budget_clipping_and_reproducibility(self):
        first = run_policy("RANDOM", "toy", self.task, self.nodes, self.adj, self.mutations, [2, 10], 2, 3, 4, 1.0, 7)
        second = run_policy("RANDOM", "toy", self.task, self.nodes, self.adj, self.mutations, [2, 10], 2, 3, 4, 1.0, 7)
        self.assertEqual([row["queried_nodes"] for row in first], [row["queried_nodes"] for row in second])
        self.assertEqual(first[-1]["requested_budget"], 10)
        self.assertEqual(first[-1]["actual_budget"], len(self.nodes) - 2)
        self.assertLessEqual(first[-1]["actual_budget"], first[-1]["requested_budget"])

    def test_policies_share_calibration(self):
        rows = []
        for method in ("RANDOM", "SURROGATE_GREEDY", "SURROGATE_BEAM", "STATIC_SURROGATE_LOOKAHEAD", "CLAIRVOYANT_REFERENCE"):
            rows.extend(run_policy(method, "toy", self.task, self.nodes, self.adj, self.mutations, [2], 3, 3, 4, 1.0, 1))
        self.assertEqual(len({row["initial_observation_nodes"] for row in rows}), 1)
        self.assertEqual(len({row["initial_observation_count"] for row in rows}), 1)

    def test_calibration_threshold_is_not_policy_success(self):
        task = {**self.task, "start_node": 5}
        row = snapshot(
            "RANDOM", "toy", task, 1, 10, [5], [], [1.5], self.nodes,
            set(self.nodes), self.adj, [], 0.0, 3, 4, 1.0, "test",
            "06_partial_observation_smoke",
        )
        self.assertTrue(row["calibration_top5_reached"])
        self.assertFalse(row["policy_top5_reached"])
        self.assertIsNone(row["queries_to_top5_after_calibration"])
        self.assertEqual(
            row["queries_to_top5_na_reason"],
            "threshold_reached_during_calibration",
        )

    def test_strict_downhill_path_requires_observed_drop_and_target(self):
        epsilon = 1e-12
        self.assertTrue(strict_downhill_path_exists(
            0, {0, 2, 4, 5}, self.adj, self.nodes, 1.5, epsilon, 3,
        ))
        self.assertFalse(strict_downhill_path_exists(
            0, {0, 1, 3, 5}, self.adj, self.nodes, 1.0, epsilon, 3,
        ))

    def test_normalized_auc_scales_calibration_improvement(self):
        self.assertAlmostEqual(
            normalized_auc_best_curve([0.0, 0.5, 1.0], 1.0, 1e-12),
            0.5,
        )
        self.assertIsNone(normalized_auc_best_curve([1.0, 1.0], 1.0, 1e-12))

    def test_stratified_task_selection_is_balanced_and_deterministic(self):
        rows = []
        for index in range(8):
            rows.append({
                "assay": "toy", "seed": 1, "task_id": f"toy:s1:t{index:03d}",
                "start_node": index, "best_reachable_fitness": 1.0,
                "valley_required": str(index < 4),
            })
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tasks.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader(); writer.writerows(rows)
            first = read_tasks(path, ["toy"], 4, "stratified_hash")
            second = read_tasks(path, ["toy"], 4, "stratified_hash")
        self.assertEqual(first, second)
        self.assertEqual(sum(row["valley_required"] for row in first), 2)
        self.assertEqual(len({row["start_node"] for row in first}), 4)

    def test_registry_append_is_idempotent(self):
        record = {
            "phase": "06_partial_observation_smoke", "assay": "toy", "task_id": "t0", "seed": 1,
            "method": "RANDOM", "requested_budget": 10, "code_version": "test", "normalized_regret": 0.5,
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "registry.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["dataset", "assay", "task_type"])
                writer.writeheader(); writer.writerow({"dataset": "ProteinGym", "assay": "old", "task_type": "oracle"})
            self.assertEqual(append_registry(path, [record]), 1)
            self.assertEqual(append_registry(path, [record]), 0)
            with path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[-1]["dataset"], "ProteinGym")
            self.assertEqual(rows[-1]["task_type"], "closed_pool_partial_observation")


if __name__ == "__main__":
    unittest.main()
