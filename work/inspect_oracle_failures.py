#!/usr/bin/env python3
"""Summarize task-level planning failures for an assay."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
from collections import Counter


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    with args.input.open(newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    def fl(r, k): return float(r[k])
    failures = []
    for r in rows:
        delta = fl(r, 'lookahead_regret') - fl(r, 'greedy_regret')
        if delta > 0:
            failures.append({
                'start_state': r['start_state'], 'oracle_target': r['oracle_target'],
                'valley_required': r['valley_required'].lower() == 'true',
                'greedy_regret': fl(r, 'greedy_regret'), 'lookahead_regret': fl(r, 'lookahead_regret'),
                'delta_regret': delta, 'lookahead_downhill_steps': int(r['lookahead_downhill_steps']),
                'greedy_downhill_steps': int(r['greedy_downhill_steps']),
                'oracle_downhill_steps': int(r['oracle_downhill_steps']),
            })
    out = {
        'input': str(args.input), 'tasks': len(rows), 'failures_lookahead_vs_greedy': len(failures),
        'failure_fraction': len(failures) / len(rows) if rows else None,
        'failure_valley_tasks': sum(x['valley_required'] for x in failures),
        'failure_downhill_pattern': {f'{k[0]}->{k[1]}': v for k, v in Counter((x['lookahead_downhill_steps'], x['greedy_downhill_steps']) for x in failures).items()},
        'largest_regret_losses': sorted(failures, key=lambda x: x['delta_regret'], reverse=True)[:20],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(out, indent=2))

if __name__ == '__main__':
    main()
