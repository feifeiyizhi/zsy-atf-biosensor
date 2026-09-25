#!/usr/bin/env python3
"""Compile a provenance-preserving, quality-gated recombined fitness table."""
from __future__ import annotations
import argparse, csv, hashlib, json, re
from collections import Counter
from pathlib import Path

AA = set('ACDEFGHIKLMNPQRSTVWY')
MUT_RE = re.compile(r'^([A-Z])([0-9]+)([A-Z])$')


def pick(row, names):
    low = {str(k).lower(): v for k, v in row.items()}
    for name in names:
        value = low.get(name.lower(), '')
        if value not in (None, ''):
            return str(value)
    return ''


def mutation_tokens(value):
    return [x.strip() for x in (value or '').replace(';', ':').replace(',', ':').split(':') if x.strip()]


def normalize(row, source_name, assay_scale):
    sequence = pick(row, ['sequence', 'reference_sequence', 'target_seq', 'aa_seq'])
    mutated = pick(row, ['mutated_sequence', 'variant_sequence', 'aa_seq'])
    mutant = pick(row, ['mutant', 'mutation', 'mutation_list'])
    variant_id = pick(row, ['variant_id', 'variant'])
    mutant = mutant or variant_id
    protein = pick(row, ['protein_id', 'uniprot_id', 'reference_name', 'dataset'])
    response = pick(row, ['response', 'DMS_score', 'fitness', 'score'])
    assay = pick(row, ['assay', 'source', 'dataset', 'reference_name']) or source_name
    parsed = mutation_tokens(mutant)
    valid_muts = []
    invalid = []
    for token in parsed:
        if token.lower() in ('wt', 'wildtype'):
            continue
        m = MUT_RE.match(token)
        if m and m.group(1) in AA and m.group(3) in AA:
            valid_muts.append(token)
        else:
            invalid.append(token)
    try:
        fitness = float(response)
        numeric = True
    except (ValueError, TypeError):
        fitness = ''
        numeric = False
    quality = []
    if not protein: quality.append('MISSING_PROTEIN_ID')
    if not mutated: quality.append('MISSING_MUTATED_SEQUENCE')
    if not numeric: quality.append('NONNUMERIC_FITNESS')
    if invalid: quality.append('INVALID_MUTATION_NOTATION')
    if not valid_muts and mutant.lower() not in ('', 'wt', 'wildtype'):
        quality.append('NO_PARSED_MUTATION')
    if sequence and mutated and len(sequence) != len(mutated): quality.append('SEQUENCE_LENGTH_MISMATCH')
    return {
        'row_id': hashlib.sha256(f'{source_name}|{protein}|{mutant}|{response}'.encode()).hexdigest()[:16],
        'source_dataset': source_name, 'protein_id': protein, 'assay': assay,
        'assay_scale': assay_scale, 'reference_sequence': sequence,
        'mutated_sequence': mutated, 'mutant': mutant,
        'mutation_count': len(valid_muts), 'response': fitness,
        'response_sd': pick(row, ['response_sd', 'fitness_sd', 'DMS_score_sd']),
        'replicate_count': pick(row, ['replicate_count', 'replicates']),
        'quality_status': 'PASS' if not quality else 'REJECT',
        'quality_reasons': ';'.join(quality),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', action='append', nargs=2, metavar=('CSV', 'SCALE'), required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    records = []
    for filename, scale in args.input:
        path = Path(filename)
        with path.open(newline='', encoding='utf-8-sig') as f:
            records.extend(normalize(row, path.name, scale) for row in csv.DictReader(f))
    fields = list(records[0]) if records else ['row_id']
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator='\n'); writer.writeheader(); writer.writerows(records)
    report = {'rows_in': len(records), 'rows_pass': sum(r['quality_status'] == 'PASS' for r in records), 'rows_reject': sum(r['quality_status'] == 'REJECT' for r in records), 'assays': len({r['assay'] for r in records}), 'proteins': len({r['protein_id'] for r in records}), 'scales': dict(Counter(r['assay_scale'] for r in records)), 'reject_reasons': dict(Counter(x for r in records for x in r['quality_reasons'].split(';') if x))}
    args.out.with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))

if __name__ == '__main__': main()
