#!/usr/bin/env python3
"""Audit real shared-volume fitness/structure assets without copying raw data."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path


def file_info(path: Path) -> dict:
    return {'path': str(path), 'exists': path.exists(), 'bytes': path.stat().st_size if path.exists() else 0}


def main():
    p = argparse.ArgumentParser(); p.add_argument('--root', type=Path, default=Path('/volume/schen04/ray/0724')); p.add_argument('--out', type=Path, required=True); a = p.parse_args()
    root = a.root
    candidates = {
        'laci_project': root / 'protein_graphwalks_structgw',
        'laci_deltaS': root / 'protein_graphwalks_structgw/out/variant_deltaS_activity.json',
        'laci_full_length_analysis': root / 'protein_graphwalks_structgw/out/fulllength_analysis.json',
        'laci_full_length_gain': root / 'protein_graphwalks_structgw/out/fulllength_gain.json',
        'laci_full_length_edges': root / 'protein_graphwalks_structgw/out/fulllength_edges.csv',
        'laci_full_length_manifest': root / 'protein_graphwalks_structgw/out/all_fulllength_manifest.json',
        'gb1_summary': root / 'paper_gw_bio/data/landscape_summary.json',
        'gb1_graph_stats': root / 'paper_gw_bio/datasets/protein_graphwalks/05_gb1_graph_stats.json',
        'gb1_data': root / 'paper_gw_bio/data/four_mutations_full_data.csv',
        'esmfold_weights': root / 'weights/esm/esmfold_3B_v1.pt',
        'laci_esmfold_structures': root / 'protein_graphwalks_structgw/out/fulllength_structures',
    }
    assets = {name: file_info(path) for name, path in candidates.items()}
    assets['laci_esmfold_structures']['file_count'] = len(list(candidates['laci_esmfold_structures'].glob('*.pdb'))) if candidates['laci_esmfold_structures'].exists() else 0
    report = {'data_root': str(root), 'public_repo_policy': 'Record paths and reviewed summaries only; do not copy raw data, PDB collections, or model weights into the public mirror.', 'assets': assets, 'scientific_boundaries': ['Full-length ESMFold local response is a structural observation, not evidence of improved held-out activity prediction.', 'GB1 measured and imputed landscapes must remain separately labeled.', 'LacI template contact features and per-variant ESMFold deltas must be evaluated as separate feature families.']}
    a.out.parent.mkdir(parents=True, exist_ok=True); a.out.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8'); print(json.dumps(report, indent=2))

if __name__ == '__main__': main()
