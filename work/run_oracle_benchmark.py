#!/usr/bin/env python3
"""Fair fixed-horizon oracle comparison on an observed fitness graph."""
from __future__ import annotations
import argparse, csv, json, random
from collections import defaultdict
from pathlib import Path


def load_graph(nodes_path: Path, edges_path: Path):
    nodes = {}
    by_depth = defaultdict(list)
    with nodes_path.open(newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            i = int(row['node_id']); row['fitness'] = float(row['fitness']); row['mutation_count'] = int(row['mutation_count'])
            nodes[i] = row; by_depth[row['mutation_count']].append(i)
    adj = defaultdict(set)
    with edges_path.open(newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            a, b = int(row['parent_id']), int(row['child_id']); adj[a].add(b); adj[b].add(a)
    return nodes, by_depth, adj


def best_path(start, horizon, nodes, adj):
    best = [start]
    def dfs(cur, path):
        nonlocal best
        if nodes[cur]['fitness'] > nodes[best[-1]]['fitness']:
            best = path[:]
        if len(path) - 1 >= horizon: return
        for nxt in sorted(adj[cur]):
            if nxt not in path: dfs(nxt, path + [nxt])
    dfs(start, [start]); return best


def greedy_path(start, horizon, nodes, adj):
    path = [start]
    for _ in range(horizon):
        options = [n for n in adj[path[-1]] if n not in path and nodes[n]['fitness'] > nodes[path[-1]]['fitness']]
        if not options: break
        path.append(max(options, key=lambda n: (nodes[n]['fitness'], -n)))
    return path


def random_path(start, horizon, nodes, adj, rng):
    path = [start]
    for _ in range(horizon):
        options = [n for n in adj[path[-1]] if n not in path]
        if not options: break
        path.append(rng.choice(sorted(options)))
    return path


def beam_path(start, horizon, width, nodes, adj):
    beam = [[start]]
    for _ in range(horizon):
        # Stopping is a valid action for every optimizing policy. Retain the
        # current beam as well as one-step extensions so paths are at most,
        # rather than exactly, ``horizon`` steps long.
        expanded = list(beam)
        for path in beam:
            for nxt in sorted(adj[path[-1]]):
                if nxt not in path: expanded.append(path + [nxt])
        expanded.sort(key=lambda p: (nodes[p[-1]]['fitness'], -len(p), -p[-1]), reverse=True)
        next_beam = expanded[:width]
        if next_beam == beam: break
        beam = next_beam
    return max(beam, key=lambda p: (nodes[p[-1]]['fitness'], -len(p), -p[-1]))


def lookahead_path(start, horizon, nodes, adj):
    # Exact bounded k-step lookahead: choose each move by the best reachable
    # terminal fitness in the remaining horizon, with no repeated genotype.
    # As for greedy and oracle, stopping before the horizon is allowed.
    def value(path, remaining):
        if remaining == 0: return nodes[path[-1]]['fitness'], path
        choices = [(nodes[path[-1]]['fitness'], path)]
        for nxt in sorted(adj[path[-1]]):
            if nxt not in path: choices.append(value(path + [nxt], remaining - 1))
        return max(choices, key=lambda x: (x[0], -len(x[1]), -x[1][-1]))
    path = [start]
    for _ in range(horizon):
        candidates = []
        for nxt in sorted(adj[path[-1]]):
            if nxt not in path:
                score, _ = value(path + [nxt], horizon - len(path))
                candidates.append((score, nxt))
        if not candidates: break
        score, nxt = max(candidates, key=lambda x: (x[0], -x[1]))
        if score <= nodes[path[-1]]['fitness']: break
        path.append(nxt)
    return path


def path_metrics(path, nodes):
    fits = [nodes[i]['fitness'] for i in path]
    deltas = [fits[i + 1] - fits[i] for i in range(len(fits) - 1)]
    return {'terminal': fits[-1], 'best_seen': max(fits), 'steps': len(path) - 1,
            'min_delta': min(deltas, default=0.0), 'downhill_steps': sum(d < 0 for d in deltas), 'path': path}


def run(nodes_path, edges_path, out, starts, horizon, seed, beam_width):
    nodes, by_depth, adj = load_graph(nodes_path, edges_path)
    rng = random.Random(seed); candidates = sorted(i for i in nodes if adj[i])
    starts_list = sorted(rng.sample(candidates, min(starts, len(candidates))))
    records = []
    for start in starts_list:
        paths = {
            'random': random_path(start, horizon, nodes, adj, rng),
            'greedy': greedy_path(start, horizon, nodes, adj),
            'beam': beam_path(start, horizon, beam_width, nodes, adj),
            'lookahead': lookahead_path(start, horizon, nodes, adj),
            'oracle': best_path(start, horizon, nodes, adj),
        }
        metrics = {name: path_metrics(path, nodes) for name, path in paths.items()}
        oracle_fit = metrics['oracle']['terminal']; start_fit = nodes[start]['fitness']
        # Target-specific valley definition: every non-repeating path to the
        # oracle target within horizon must contain a strict downhill edge.
        target = paths['oracle'][-1]; monotonic = False
        stack = [(start, [start])]
        while stack:
            cur, path = stack.pop()
            if cur == target: monotonic = True; break
            if len(path) - 1 >= horizon: continue
            for nxt in adj[cur]:
                if nxt not in path and nodes[nxt]['fitness'] >= nodes[cur]['fitness']:
                    stack.append((nxt, path + [nxt]))
        valley = oracle_fit > start_fit and not monotonic
        rec = {'start_state': start, 'start_fitness': start_fit, 'oracle_target': target, 'valley_required': valley}
        for name, m in metrics.items():
            for key in ('terminal', 'best_seen', 'steps', 'min_delta', 'downhill_steps'):
                rec[f'{name}_{key}'] = m[key]
            rec[f'{name}_path'] = json.dumps(m['path'], separators=(',', ':'))
            rec[f'{name}_regret'] = oracle_fit - m['terminal']
        records.append(rec)
    out.parent.mkdir(parents=True, exist_ok=True)
    csv_out = out.with_suffix('.csv')
    with csv_out.open('w', newline='', encoding='utf-8') as f:
        fields = list(records[0]) if records else ['start_state']; w = csv.DictWriter(f, fieldnames=fields, lineterminator='\n'); w.writeheader(); w.writerows(records)
    summary = {'nodes': len(nodes), 'edges': sum(map(len, adj.values())) // 2, 'tasks': len(records), 'horizon': horizon, 'beam_width': beam_width, 'seed': seed, 'valley_required_tasks': sum(r['valley_required'] for r in records)}
    for name in ('random', 'greedy', 'beam', 'lookahead', 'oracle'):
        vals = [r[f'{name}_terminal'] for r in records]
        summary[f'{name}_mean_terminal'] = sum(vals) / len(vals) if vals else None
        summary[f'{name}_mean_regret'] = sum(r[f'{name}_regret'] for r in records) / len(records) if records else None
        summary[f'{name}_valley_success'] = sum(r['valley_required'] and r[f'{name}_regret'] == 0 for r in records)
    out.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8'); print(json.dumps(summary, indent=2))


def main():
    p = argparse.ArgumentParser(); p.add_argument('--nodes', type=Path, required=True); p.add_argument('--edges', type=Path, required=True); p.add_argument('--out', type=Path, required=True); p.add_argument('--starts', type=int, default=500); p.add_argument('--horizon', type=int, default=4); p.add_argument('--beam-width', type=int, default=8); p.add_argument('--seed', type=int, default=42); a = p.parse_args(); run(a.nodes, a.edges, a.out, a.starts, a.horizon, a.seed, a.beam_width)
if __name__ == '__main__': main()
