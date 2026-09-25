#!/usr/bin/env python3
"""Oracle greedy versus bounded lookahead on an observed fitness graph."""
from __future__ import annotations
import argparse,csv,json,random,time
from collections import defaultdict, deque
from pathlib import Path

def run(nodes_path: Path, edges_path: Path, out: Path, starts: int, horizon: int, seed: int) -> dict:
    nodes={}; by_depth=defaultdict(list)
    with nodes_path.open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            i=int(r['node_id']); r['fitness']=float(r['fitness']); r['mutation_count']=int(r['mutation_count']); nodes[i]=r; by_depth[r['mutation_count']].append(i)
    adj=defaultdict(set)
    with edges_path.open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            a,b=int(r['parent_id']),int(r['child_id']); adj[a].add(b); adj[b].add(a)
    rng=random.Random(seed); candidates=[i for d in by_depth for i in by_depth[d] if adj[i]]
    starts_list=sorted(rng.sample(candidates,min(starts,len(candidates))))
    records=[]
    for s in starts_list:
        fit0=nodes[s]['fitness']
        # Enumerate simple paths up to horizon; this is exact for the sampled task.
        best=-float('inf'); best_path=[s]; monotonic_best=-float('inf'); mono_path=[s]
        valley_target=None
        def dfs(cur,path):
            nonlocal best,best_path,monotonic_best,mono_path
            val=nodes[cur]['fitness']
            if val>best: best,best_path=val,path[:]
            if val>=fit0 and val>monotonic_best: monotonic_best,mono_path=val,path[:]
            if len(path)-1>=horizon:return
            for nxt in adj[cur]:
                if nxt in path: continue
                dfs(nxt,path+[nxt])
        dfs(s,[s])
        # Greedy hill climb with deterministic tie-breaking.
        cur=s; gpath=[s]
        for _ in range(horizon):
            options=[n for n in adj[cur] if nodes[n]['fitness']>nodes[cur]['fitness']]
            if not options: break
            nxt=max(options,key=lambda n:(nodes[n]['fitness'],-n)); gpath.append(nxt); cur=nxt
        target=best_path[-1]
        target_fit=best
        # Valley-required is target-specific: no non-decreasing simple path to target.
        mono_to_target=False
        q=[(s,[s])]
        while q:
            cur,path=q.pop()
            if cur==target: mono_to_target=True; break
            if len(path)-1>=horizon: continue
            for nxt in adj[cur]:
                if nxt not in path and nodes[nxt]['fitness']>=nodes[cur]['fitness']:
                    q.append((nxt,path+[nxt]))
        valley_required=target_fit>fit0 and not mono_to_target
        min_drop=min((nodes[gpath[i+1]]['fitness']-nodes[gpath[i]]['fitness'] for i in range(len(gpath)-1)),default=0)
        min_drop_plan=min((nodes[best_path[i+1]]['fitness']-nodes[best_path[i]]['fitness'] for i in range(len(best_path)-1)),default=0)
        records.append({'start_state':s,'algorithm_greedy_terminal':nodes[gpath[-1]]['fitness'],'oracle_terminal':target_fit,'greedy_best_seen':max(nodes[i]['fitness'] for i in gpath),'oracle_best_seen':target_fit,'regret':target_fit-nodes[gpath[-1]]['fitness'],'greedy_steps':len(gpath)-1,'oracle_steps':len(best_path)-1,'valley_required':valley_required,'oracle_valley_crossed':min_drop_plan<0,'greedy_min_step_delta':min_drop,'oracle_min_step_delta':min_drop_plan,'target_state':target})
    out.parent.mkdir(parents=True,exist_ok=True); out.with_suffix('.csv').write_text('')
    with out.with_suffix('.csv').open('w',newline='',encoding='utf-8') as f:
        fields=list(records[0]) if records else ['start_state']; w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n'); w.writeheader(); w.writerows(records)
    summary={'nodes':len(nodes),'edges':sum(map(len,adj.values()))//2,'tasks':len(records),'horizon':horizon,'seed':seed,'valley_required_tasks':sum(r['valley_required'] for r in records),'greedy_mean_terminal':sum(r['algorithm_greedy_terminal'] for r in records)/len(records) if records else None,'oracle_mean_terminal':sum(r['oracle_terminal'] for r in records)/len(records) if records else None,'mean_regret':sum(r['regret'] for r in records)/len(records) if records else None,'valley_greedy_failures':sum(r['valley_required'] and r['regret']>0 for r in records)}
    out.write_text(json.dumps(summary,indent=2)+'\n'); print(json.dumps(summary,indent=2)); return summary

def main():
    p=argparse.ArgumentParser(); p.add_argument('--nodes',type=Path,required=True); p.add_argument('--edges',type=Path,required=True); p.add_argument('--out',type=Path,required=True); p.add_argument('--starts',type=int,default=500); p.add_argument('--horizon',type=int,default=4); p.add_argument('--seed',type=int,default=42); a=p.parse_args(); run(a.nodes,a.edges,a.out,a.starts,a.horizon,a.seed)
if __name__=='__main__':main()
