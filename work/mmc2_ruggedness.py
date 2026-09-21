#!/usr/bin/env python3
"""Reproduce mmc2-style normalized graph Dirichlet ruggedness."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd


def knn_adjacency(seqs: list[str], k: int) -> np.ndarray:
    n = len(seqs)
    # Hamming distance equals Euclidean distance up to a constant for one-hot.
    x = np.fromiter((ord(c) for s in seqs for c in s), dtype=np.int16).reshape(n, -1)
    dist = (x[:, None, :] != x[None, :, :]).sum(axis=2)
    a = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        order = np.argsort(dist[i], kind="stable")
        a[i, order[1:k + 1]] = 1.0
    # symmetric kNN graph as in the paper
    a = np.maximum(a, a.T)
    np.fill_diagonal(a, 0)
    return a


def energy(a: np.ndarray, y: np.ndarray) -> float:
    deg = a.sum(axis=1)
    lap = np.diag(deg) - a
    return float(y @ lap @ y / (len(y) * y @ y))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--k", type=int, default=34)
    p.add_argument("--replicates", type=int, default=1000)
    p.add_argument("--fraction", type=float, default=.95)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(args.input)
    results = []; summary = []
    rng = np.random.default_rng(args.seed)
    for dataset, g0 in d.groupby("dataset", sort=False):
        for length, g in g0.groupby("sequence_length", sort=False):
            if len(g) < max(args.k + 1, 10): continue
            g = g.reset_index(drop=True); seqs = g.aa_seq.tolist(); y = g.fitness.to_numpy(float)
            a = knn_adjacency(seqs, min(args.k, len(g)-1)); e = energy(a, y)
            vals = []
            m = max(args.k + 1, int(round(len(g) * args.fraction)))
            for _ in range(args.replicates):
                ix = rng.choice(len(g), size=m, replace=False)
                vals.append(energy(a[np.ix_(ix, ix)], y[ix]))
            summary.append({"dataset": dataset, "sequence_length": int(length), "n": len(g),
                            "k": args.k, "fraction": args.fraction, "full_energy": e,
                            "subsample_mean": float(np.mean(vals)), "subsample_sd": float(np.std(vals, ddof=1)),
                            "subsample_n": args.replicates})
            for i, row in g.iterrows():
                results.append({"dataset": dataset, "sequence_length": int(length),
                                "variant_id": row.variant_id, "fitness": row.fitness,
                                "degree": int(a[i].sum()), "energy_full": e})
    pd.DataFrame(results).to_csv(args.out / "mmc2_ruggedness_nodes.csv", index=False)
    pd.DataFrame(summary).to_csv(args.out / "mmc2_ruggedness_summary.csv", index=False)
    report = {"k": args.k, "fraction": args.fraction, "replicates": args.replicates,
              "seed": args.seed, "strata": summary}
    (args.out / "mmc2_ruggedness.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))

if __name__ == "__main__": main()
