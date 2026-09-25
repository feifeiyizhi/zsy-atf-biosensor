#!/usr/bin/env python3
"""Build observed Hamming-1 mutation graphs for eligible fitness assays."""
from __future__ import annotations
import argparse, csv, json, re
from collections import defaultdict, deque
from pathlib import Path


def parse_mutations(value: str) -> frozenset[str]:
    value = (value or '').strip()
    if not value or value.upper() in {'WT', 'WILDTYPE'}:
        return frozenset()
    return frozenset(x.strip() for x in value.replace(';', ':').replace(',', ':').split(':') if x.strip())


def build(rows: list[dict[str, str]], out: Path) -> dict[str, object]:
    nodes_by_mut: dict[frozenset[str], dict[str, object]] = {}
    for row in rows:
        try:
            fitness = float(row['response'])
        except (KeyError, TypeError, ValueError):
            continue
        muts = parse_mutations(row.get('mutant', ''))
        if muts not in nodes_by_mut or len(row.get('mutated_sequence', '')) > len(str(nodes_by_mut[muts].get('sequence', ''))):
            nodes_by_mut[muts] = {
                'mutations': ':'.join(sorted(muts)),
                'sequence': row.get('mutated_sequence', ''),
                'fitness': fitness,
                'mutation_count': len(muts),
                'protein_id': row.get('protein_id', ''),
                'uniprot_id': row.get('uniprot_id', ''),
                'assay': row.get('assay', ''),
            }
    ordered = sorted(nodes_by_mut.items(), key=lambda x: (len(x[0]), tuple(sorted(x[0]))))
    index = {muts: i for i, (muts, _) in enumerate(ordered)}
    edges: list[tuple[int, int, str]] = []
    adjacency: dict[int, set[int]] = defaultdict(set)
    for child_muts, child_id in index.items():
        for mutation in child_muts:
            parent_muts = child_muts - {mutation}
            if parent_muts in index:
                parent_id = index[parent_muts]
                edges.append((parent_id, child_id, mutation))
                adjacency[parent_id].add(child_id)
                adjacency[child_id].add(parent_id)
    out.mkdir(parents=True, exist_ok=True)
    fields = ['node_id', 'mutations', 'sequence', 'fitness', 'mutation_count', 'protein_id', 'uniprot_id', 'assay']
    with (out / 'nodes.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator='\n'); w.writeheader()
        for node_id, (_, node) in enumerate(ordered):
            w.writerow({'node_id': node_id, **node})
    with (out / 'edges.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, lineterminator='\n'); w.writerow(['parent_id', 'child_id', 'mutation'])
        w.writerows(edges)
    connected = set()
    largest_cc = 0
    for start in range(len(ordered)):
        if start in connected or start not in adjacency:
            continue
        q = [start]; connected.add(start); size = 0
        while q:
            cur = q.pop(); size += 1
            for nxt in adjacency[cur]:
                if nxt not in connected: connected.add(nxt); q.append(nxt)
        largest_cc = max(largest_cc, size)
    local_max = local_min = improving = 0
    for node_id, (_, node) in enumerate(ordered):
        neigh = [float(ordered[n][1]['fitness']) for n in adjacency[node_id]]
        if not neigh: continue
        fit = float(node['fitness'])
        local_max += fit >= max(neigh)
        local_min += fit <= min(neigh)
        improving += any(v > fit for v in neigh)
    summary = {
        'assay_id': out.name, 'N_nodes': len(ordered), 'N_edges': len(edges),
        'N_connected_nodes': len(connected), 'largest_connected_component_size': largest_cc,
        'largest_connected_component_fraction': largest_cc / len(ordered) if ordered else 0,
        'max_mutation_depth': max((len(m) for m, _ in ordered), default=0),
        'local_maxima': local_max, 'local_minima': local_min,
        'fraction_nodes_with_improving_neighbors': improving / len(ordered) if ordered else 0,
        'fitness_min': min((float(n['fitness']) for _, n in ordered), default=None),
        'fitness_max': max((float(n['fitness']) for _, n in ordered), default=None),
    }
    (out / 'landscape_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--assay', action='append', required=True, help='Exact source suffix or substring; repeatable')
    ap.add_argument('--out-root', type=Path, required=True)
    args = ap.parse_args()
    wanted = set(args.assay)
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    with args.input.open(newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            source = row.get('source', '')
            for key in wanted:
                if source == key or key in source:
                    grouped[key].append(row); break
    if not grouped:
        raise SystemExit(f'No requested assays found: {sorted(wanted)}')
    reports = []
    for key, rows in grouped.items():
        safe = re.sub(r'[^A-Za-z0-9_.-]+', '_', key.replace('ProteinGym:', '').removesuffix('.csv'))
        reports.append(build(rows, args.out_root / safe))
    print(json.dumps(reports, indent=2))

if __name__ == '__main__': main()
