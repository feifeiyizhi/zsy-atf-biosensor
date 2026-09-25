#!/usr/bin/env python3
"""Compile assay-level fitness landscape and GraphWalk eligibility statistics."""
from __future__ import annotations
import argparse, csv, json, math
from collections import Counter, defaultdict
from pathlib import Path


def parse_mutations(value: str) -> frozenset[str]:
    value = (value or '').strip()
    if not value or value.upper() in {'WT', 'WILDTYPE'}:
        return frozenset()
    return frozenset(x.strip() for x in value.replace(';', ':').replace(',', ':').split(':') if x.strip())


def stats(values: list[float]) -> tuple[float, float, float, float]:
    if not values:
        return (float('nan'),) * 4
    mean = sum(values) / len(values)
    var = sum((x - mean) ** 2 for x in values) / len(values)
    return min(values), max(values), mean, math.sqrt(var)


def csv_float(value: float) -> str:
    return '' if math.isnan(value) else f'{value:.8g}'


def compile_assay(assay: str, rows: list[dict[str, str]]) -> dict[str, object]:
    # Keep the best-defined row for a duplicate genotype; duplicate assay rows
    # are not additional graph nodes.
    nodes: dict[frozenset[str], dict[str, object]] = {}
    values: list[float] = []
    for row in rows:
        try:
            fitness = float(row['response'])
        except (KeyError, TypeError, ValueError):
            continue
        muts = parse_mutations(row.get('mutant', ''))
        old = nodes.get(muts)
        if old is None or len(row.get('mutated_sequence', '')) > len(str(old.get('mutated_sequence', ''))):
            nodes[muts] = {'fitness': fitness, 'sequence': row.get('mutated_sequence', ''), 'mutant': row.get('mutant', '')}
    values = [float(n['fitness']) for n in nodes.values()]
    depths = Counter(len(m) for m in nodes)
    edge_count = 0
    adjacency: dict[frozenset[str], set[frozenset[str]]] = defaultdict(set)
    for muts in nodes:
        for mutation in muts:
            parent = muts - {mutation}
            if parent in nodes:
                edge_count += 1
                adjacency[muts].add(parent)
                adjacency[parent].add(muts)
    seen: set[frozenset[str]] = set()
    largest_cc = 0
    for start in nodes:
        if start in seen:
            continue
        stack = [start]; seen.add(start); size = 0
        while stack:
            cur = stack.pop(); size += 1
            for nxt in adjacency[cur]:
                if nxt not in seen:
                    seen.add(nxt); stack.append(nxt)
        largest_cc = max(largest_cc, size)
    local_max = 0; local_min = 0; improving = 0
    for genotype, node in nodes.items():
        neigh = [nodes[x]['fitness'] for x in adjacency[genotype]]
        if not neigh:
            continue
        fit = float(node['fitness'])
        local_max += fit >= max(neigh)
        local_min += fit <= min(neigh)
        improving += any(x > fit for x in neigh)
    fmin, fmax, fmean, fstd = stats(values)
    n_multi = sum(len(m) > 1 for m in nodes)
    n_single = sum(len(m) == 1 for m in nodes)
    reasons = []
    if not values: reasons.append('INSUFFICIENT_FITNESS')
    if len(nodes) < 20: reasons.append('INSUFFICIENT_VARIANTS')
    if n_multi == 0: reasons.append('INSUFFICIENT_MULTIMUTANTS')
    if edge_count < 10: reasons.append('DISCONNECTED')
    if max(depths, default=0) < 2: reasons.append('INSUFFICIENT_DEPTH')
    if fstd == 0 or math.isnan(fstd): reasons.append('INSUFFICIENT_FITNESS')
    eligible = not reasons
    return {
        'assay_id': assay,
        'protein_id': rows[0].get('protein_id', '') if rows else '',
        'uniprot_id': rows[0].get('uniprot_id', '') if rows else '',
        'taxon': rows[0].get('taxon', '') if rows else '',
        'N_variants': len(nodes), 'N_unique_sequences': len({n['sequence'] for n in nodes.values() if n['sequence']}),
        'min_mutation_count': min(depths, default=0), 'max_mutation_count': max(depths, default=0),
        'mutation_count_distribution': json.dumps(dict(sorted(depths.items())), sort_keys=True),
        'fitness_nonnull_fraction': len(values) / len(rows) if rows else 0,
        'fitness_min': csv_float(fmin), 'fitness_max': csv_float(fmax), 'fitness_mean': csv_float(fmean), 'fitness_std': csv_float(fstd),
        'N_single_mutants': n_single, 'N_multi_mutants': n_multi, 'N_hamming1_edges': edge_count,
        'largest_connected_component_size': largest_cc,
        'largest_connected_component_fraction': largest_cc / len(nodes) if nodes else 0,
        'max_observed_mutation_depth': max(depths, default=0),
        'local_maxima': local_max, 'local_minima': local_min,
        'fraction_nodes_with_improving_neighbors': improving / len(nodes) if nodes else 0,
        'graphwalk_eligible': 'GRAPHWALK_ELIGIBLE' if eligible else 'PREDICTION_ONLY',
        'eligibility_reason': 'GRAPHWALK_ELIGIBLE' if eligible else ';'.join(reasons),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    with args.input.open(newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            groups[row.get('source', row.get('assay', 'unknown'))].append(row)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    records = [compile_assay(k, v) for k, v in sorted(groups.items())]
    fields = list(records[0]) if records else ['assay_id']
    with args.out.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(records)
    summary = Counter(r['graphwalk_eligible'] for r in records)
    report = {'input': str(args.input), 'assays': len(records), 'eligibility_counts': dict(summary), 'manifest': str(args.out)}
    args.out.with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))

if __name__ == '__main__': main()
