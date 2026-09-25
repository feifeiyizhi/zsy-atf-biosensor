#!/usr/bin/env python3
"""Collect public mutation-response datasets into a diffusion-training manifest.

This script intentionally keeps downloaded archives/cache outside the public repo.
The manifest preserves both reference `sequence` when available and the observed
`mutated_sequence`; ProteinGym substitution files commonly omit the wild-type
sequence and replicate uncertainty, so those fields remain empty rather than
being inferred or fabricated.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
import urllib.request
import zipfile
from pathlib import Path

DEFAULT_URL = "https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/DMS_ProteinGym_substitutions.zip"


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url}")
    with urllib.request.urlopen(url, timeout=60) as r, dest.open("wb") as f:
        shutil.copyfileobj(r, f)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def split_mutants(value: str) -> list[str]:
    value = (value or "").strip()
    if not value or value.upper() in {"WT", "WILDTYPE"}:
        return []
    return [x for x in value.replace(";", ":").replace(",", ":").split(":") if x]


def pick(row: dict[str, str], names: list[str]) -> str:
    lower = {k.lower(): v for k, v in row.items()}
    for name in names:
        if name.lower() in lower and lower[name.lower()] not in (None, ""):
            return str(lower[name.lower()])
    return ""


def split_name(protein_id: str) -> str:
    """Assign a protein deterministically; no protein appears in multiple splits."""
    bucket = int(hashlib.sha256(protein_id.encode("utf-8")).hexdigest()[:8], 16) % 100
    if bucket < 80:
        return "train"
    if bucket < 90:
        return "validation"
    return "test"


def normalize_csv(path: Path, source: str, writer: csv.DictWriter, limit: int | None, split_counts: dict[str, int], missing_counts: dict[str, int]) -> int:
    count = 0
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            protein = pick(row, ["protein_name", "protein_id", "DMS_id", "assay_id"])
            if not protein:
                protein = path.stem
            reference_sequence = pick(row, ["target_seq", "reference_sequence", "wild_type_sequence"])
            mutated_sequence = pick(row, ["mutated_sequence", "variant_sequence"])
            sequence = reference_sequence
            mutant = pick(row, ["mutant", "mutant_name", "mutation", "variant"])
            response = pick(row, ["DMS_score", "fitness", "score", "mean_score", "measurement"])
            assay = pick(row, ["assay", "DMS_description", "description", "selection"]) or path.stem
            if not mutant or not response:
                continue
            try:
                float(response)
            except ValueError:
                continue
            mutations = split_mutants(mutant)
            row_out = {
                "protein_id": protein,
                "sequence": sequence,
                "mutated_sequence": mutated_sequence,
                "mutant": mutant,
                "mutation_count": len(mutations),
                "response": response,
                "response_sd": pick(row, ["DMS_score_sd", "fitness_sd", "score_sd", "measurement_sd"]),
                "replicate_count": pick(row, ["replicate_count", "replicates", "n_replicates"]),
                "assay": assay,
                "split": split_name(protein),
                "source": source,
            }
            for field in ("sequence", "response_sd", "replicate_count"):
                if not row_out[field]:
                    missing_counts[field] += 1
            writer.writerow(row_out)
            split_counts[split_name(protein)] += 1
            count += 1
            if limit is not None and count >= limit:
                break
    return count


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default=DEFAULT_URL)
    ap.add_argument("--cache", type=Path, default=Path("/mnt/workspace/ray/mutation_response_cache"))
    ap.add_argument("--out", type=Path, default=Path("work/results/mutation_response_manifest.csv"))
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--keep-archive", action="store_true")
    args = ap.parse_args()

    archive = args.cache / Path(args.url).name
    if not archive.exists():
        try:
            download(args.url, archive)
        except Exception as e:
            raise SystemExit(
                f"Could not download {args.url}: {e}\n"
                "Pass --url to the current ProteinGym archive or a permitted DMS archive."
            )
    digest = sha256(archive)
    extract = args.cache / archive.stem
    if not extract.exists():
        extract.mkdir(parents=True)
        with zipfile.ZipFile(archive) as z:
            z.extractall(extract)

    files = sorted(extract.rglob("*.csv"))
    if not files:
        raise SystemExit(f"No CSV files found after extracting {archive}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["protein_id", "sequence", "mutated_sequence", "mutant", "mutation_count", "response", "response_sd", "replicate_count", "assay", "split", "source"]
    total = 0
    split_counts = {"train": 0, "validation": 0, "test": 0}
    missing_counts = {"sequence": 0, "response_sd": 0, "replicate_count": 0}
    with args.out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for path in files:
            remaining = None if args.limit is None else max(args.limit - total, 0)
            if remaining == 0:
                break
            total += normalize_csv(path, f"ProteinGym:{path.name}", writer, remaining, split_counts, missing_counts)
    report = {
        "source_url": args.url,
        "archive": str(archive),
        "archive_sha256": digest,
        "csv_files_seen": len(files),
        "rows_written": total,
        "split_counts": split_counts,
        "missing_counts": missing_counts,
        "manifest": str(args.out),
        "public_repo_policy": "Do not commit raw archive or downloaded source tables; commit only schema/code and reviewed derived summaries.",
    }
    args.out.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if not args.keep_archive:
        print(f"Archive retained in external cache: {archive}")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
