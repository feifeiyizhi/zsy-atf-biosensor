#!/usr/bin/env python3
"""Prepare and audit the LacI/GalR fitness data shipped with mmc2.

The source files are kept outside the public repository.  This script writes
only derived tables and graph summaries, so it can be rerun after the source
location changes without copying raw data into git.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

AA = set("ACDEFGHIKLMNPQRSTVWY")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path, required=True,
                   help="directory containing anclaci/.../csvs and sequences")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--knn", type=int, default=34)
    return p.parse_args()


def find_source(root: Path, filename: str, subdir: str) -> Path:
    hits = list(root.glob(f"**/{subdir}/{filename}"))
    if not hits:
        raise FileNotFoundError(f"cannot find {subdir}/{filename} below {root}")
    return hits[0]


def clean_seq(x) -> str:
    return "".join(str(x).upper().split())


def mutations(seq: str, ref: str):
    if len(seq) != len(ref):
        return None
    return [f"{a}{i + 1}{b}" for i, (a, b) in enumerate(zip(ref, seq)) if a != b]


def load_dataset(source: Path, kind: str) -> pd.DataFrame:
    if kind == "dms":
        p = find_source(source, "DMS_log2_normalized_enrichment.csv", "csvs")
        d = pd.read_csv(p)
        ref = clean_seq(d.loc[d.lib_id.astype(str).str.lower().eq("wt"), "aa_seq"].iloc[0])
        d["dataset"] = "DMS"
        d["fitness"] = d["avg_log2_normalized_enrichment"]
        d["fitness_type"] = "log2_normalized_enrichment"
        d["fitness_sd"] = d["sd_log2_normalized_enrichment"]
        d["reference_sequence"] = ref
        d["mutation_list"] = d.aa_seq.map(lambda x: mutations(clean_seq(x), ref))
        d["reference_name"] = "EcLacI_DBD_WT"
    else:
        p = find_source(source, "LGF_log2_enrichment.csv", "csvs")
        d = pd.read_csv(p)
        d["dataset"] = "LGF"
        d["fitness"] = d["avg_log2_enrichment"]
        d["fitness_type"] = "log2_enrichment"
        d["fitness_sd"] = d["sd_log2_enrichment"]
        # LGF contains naturally varying DBD lengths (indels are present), so
        # a raw consensus coordinate is unsafe.  Preserve sequences exactly;
        # mutation_list remains empty for LGF until an explicit alignment is
        # supplied.  Graphs are built within equal-length strata below.
        seqs = d.aa_seq.map(clean_seq).tolist()
        d["reference_sequence"] = ""
        d["mutation_list"] = ""
        d["reference_name"] = "LGF_no_common_reference"
    d["aa_seq"] = d.aa_seq.map(clean_seq)
    d["sequence_length"] = d.aa_seq.str.len()
    d["mutation_count"] = d.mutation_list.map(lambda x: -1 if x is None else len(x))
    d["mutation_list"] = d.mutation_list.map(lambda x: "" if x is None else ";".join(x))
    d["variant_id"] = d["lib_id"].astype(str)
    d["replicate_count"] = 2
    cols = ["dataset", "variant_id", "aa_seq", "reference_sequence", "reference_name",
            "mutation_list", "mutation_count", "sequence_length", "fitness",
            "fitness_type", "fitness_sd", "replicate_count"]
    return d[cols]


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b)) if len(a) == len(b) else 10**9


def graph_summary(d: pd.DataFrame, k: int) -> pd.DataFrame:
    rows = []
    for dataset, g0 in d.groupby("dataset", sort=False):
        # Hamming distance is defined only within equal-length sequences.  LGF
        # therefore yields one graph component per observed length stratum.
        for seq_len, g in g0.groupby("sequence_length", sort=False):
            if len(g) < 2:
                continue
            seqs = g.aa_seq.tolist(); y = g.fitness.to_numpy(float)
            adj = defaultdict(set)
            for i, s in enumerate(seqs):
                ds = sorted((hamming(s, t), j) for j, t in enumerate(seqs) if i != j)
                for dist, j in ds[:min(k, len(ds))]:
                    adj[i].add(j); adj[j].add(i)
            deg = np.array([len(adj[i]) for i in range(len(seqs))])
            higher = np.array([sum(y[j] > y[i] for j in adj[i]) for i in range(len(seqs))])
            local_peak = (deg > 0) & (higher == 0)
            # Normalized Dirichlet energy on the undirected kNN graph.
            edge_sum = 0.0; edge_n = 0
            for i in range(len(seqs)):
                for j in adj[i]:
                    if j > i:
                        edge_sum += (y[i] - y[j]) ** 2; edge_n += 1
            energy = edge_sum / max(edge_n, 1)
            for i, row_id in enumerate(g.index):
                rows.append({"dataset": dataset, "sequence_length": int(seq_len),
                             "variant_id": g.loc[row_id, "variant_id"],
                             "knn_degree": int(deg[i]), "knn_higher_neighbors": int(higher[i]),
                             "local_peak": bool(local_peak[i]), "fitness": float(y[i]),
                             "graph_dirichlet_energy": float(energy), "knn": k})
    return pd.DataFrame(rows)


def main():
    args = parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    dms = load_dataset(args.source, "dms")
    lgf = load_dataset(args.source, "lgf")
    unified = pd.concat([dms, lgf], ignore_index=True)
    unified.to_csv(args.out / "mmc2_unified_variants.csv", index=False)
    dms.to_csv(args.out / "mmc2_dms_single_mutants.csv", index=False)
    lgf.to_csv(args.out / "mmc2_lgf_phylogenetic_variants.csv", index=False)
    graph = graph_summary(unified, args.knn)
    graph.to_csv(args.out / "mmc2_graph_metrics.csv", index=False)
    report = {
        "source": str(args.source), "rows": int(len(unified)),
        "datasets": {k: int(v) for k, v in unified.dataset.value_counts().items()},
        "bad_sequence_lengths": int((unified.sequence_length <= 0).sum()),
        "dms_non_single_mutants": int(((dms.mutation_count != 1) &
                                       (dms.variant_id.str.lower() != "wt")).sum()),
        "dms_wt_rows": int((dms.variant_id.str.lower() == "wt").sum()),
        "lgf_consensus_coordinate_only": True,
        "graph_knn": args.knn,
        "graph_local_peaks": {k: int(v) for k, v in graph[graph.local_peak].dataset.value_counts().items()},
        "fitness_ranges": {k: [float(v.fitness.min()), float(v.fitness.max())]
                           for k, v in unified.groupby("dataset")},
    }
    (args.out / "mmc2_validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
