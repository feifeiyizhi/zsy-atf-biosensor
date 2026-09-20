#!/usr/bin/env python3
"""build_evolution_graph.py — 338lib 富集轨迹可达路径图(GraphWalks 真命题)。
边 ≠ Hamming-1 序列邻接(分析者构造),而是【突变阶相邻 nested subset+1】=真实进化步:
  child 的突变集 ⊃ parent 的突变集 且 |child|-|parent|=1 (加一个突变)。
  定向进化就是逐步累积突变,parent→child(多一个突变且包含全部旧突变)才是可观测的真进化边。
节点 = 可靠变体(reliable=True 或 sum_count>=阈值),适应度 = fitness_slope(逐轮富集斜率)。
核心分析:高适应度变体的可达路径是否破 greedy myopic ascent
  (greedy 每步选 df 最大的邻居;若存在"必须先降后升"才能到达的高峰 = valley/破greedy)。
产出: results/evolution_graph.json + evolution_graph_report.md
"""
import csv, numpy as np, re
from collections import defaultdict, Counter
D='./work'
rows=[r for r in csv.DictReader(open(f'{D}/results/genotype_17site_landscape.csv'))]
def parse_muts(r):
    m=r['mutations'].strip()
    if not m: return frozenset()
    out=[]
    for x in m.split(';'):
        if re.match(r'[A-Z]\d+[A-Z]',x.strip()): out.append(x.strip())
    return frozenset(out)

# 节点池:可靠变体(reliable=True)
pool=[]
for r in rows:
    if r['reliable']=='True' and r['fitness_slope']:
        pool.append(dict(muts=parse_muts(r), fit=float(r['fitness_slope']),
                         n_mut=int(r['n_mut']), sum_count=int(r['sum_count'])))
print(f"可靠节点: {len(pool)}", flush=True)
# 去重(同突变集)
seen={}; nodes=[]
for p in pool:
    k=p['muts']
    if k not in seen: seen[k]=len(nodes); nodes.append(p)
print(f"去重后节点: {len(nodes)}", flush=True)
mset2idx={n['muts']:i for i,n in enumerate(nodes)}

# 突变阶相邻边: 对每个节点(k突变), 去掉一个突变→(k-1)集合, 若在池中则连边(parent=k-1, child=k)
edges=[]  # (parent_idx, child_idx)
for i,n in enumerate(nodes):
    ms=n['muts']
    if len(ms)==0: continue
    for m in ms:
        parent=ms - {m}
        if parent in mset2idx:
            j=mset2idx[parent]
            edges.append((j,i))  # parent j → child i (child 多一个突变 m)
print(f"突变阶相邻边(nested subset+1): {len(edges)}", flush=True)

# 图统计
adj_out=defaultdict(list)  # parent→[children]
adj_in=defaultdict(list)   # child→[parents]
for j,i in edges:
    adj_out[j].append(i); adj_in[i].append(j)
deg_out=[len(adj_out[i]) for i in range(len(nodes))]
deg_in=[len(adj_in[i]) for i in range(len(nodes))]
n_connected=sum(1 for i in range(len(nodes)) if deg_out[i]+deg_in[i]>0)
print(f"有边节点: {n_connected}/{len(nodes)} ({100*n_connected/len(nodes):.1f}%)", flush=True)
print(f"出度: 中位{np.median(deg_out):.1f} 均值{np.mean(deg_out):.2f} 最大{max(deg_out)}", flush=True)
print(f"入度: 中位{np.median(deg_in):.1f} 均值{np.mean(deg_in):.2f} 最大{max(deg_in)}", flush=True)

# 沿边适应度变化 df = fit(child) - fit(parent)
dfs=[nodes[i]['fit']-nodes[j]['fit'] for j,i in edges]
dfs=np.array(dfs)
print(f"\n边适应度变化 df=fit(child)-fit(parent):", flush=True)
print(f"  上坡(df>0): {(dfs>0).sum()}/{len(dfs)} ({100*(dfs>0).mean():.1f}%)", flush=True)
print(f"  下坡(df<0): {(dfs<0).sum()}/{len(dfs)} ({100*(dfs<0).mean():.1f}%)", flush=True)
print(f"  df 中位{np.median(dfs):.3f} 范围[{dfs.min():.2f},{dfs.max():.2f}]", flush=True)

# 破 greedy 分析: 找有向路径(按突变阶递增,parent→child),
# 对每个 child, 到它的路径是否必须经过 df<0 的边(valley step)才能达到更高 fit
# 简化: 找"局部谷"节点 = 该节点 fit 低于其所有 parent, 但它有 child fit 更高(超过所有 parent)
valley_nodes=0; valley_examples=[]
for i,n in enumerate(nodes):
    parents=adj_in[i]; children=adj_out[i]
    if not parents or not children: continue
    max_parent_fit=max(nodes[j]['fit'] for j in parents)
    # 此节点是"谷": 比所有 parent 低
    if n['fit'] < max_parent_fit:
        # 但它能通向比 parent 更高的 child(过谷升天)
        best_child_fit=max(nodes[c]['fit'] for c in children)
        if best_child_fit > max_parent_fit:
            valley_nodes+=1
            if len(valley_examples)<8:
                # 找具体路径 parent→valley→child
                bp=max(parents,key=lambda j:nodes[j]['fit'])
                bc=max(children,key=lambda c:nodes[c]['fit'])
                added_p2v=list(nodes[i]['muts']-nodes[bp]['muts'])
                added_v2c=list(nodes[bc]['muts']-nodes[i]['muts'])
                valley_examples.append(dict(
                    parent_fit=round(nodes[bp]['fit'],3), valley_fit=round(n['fit'],3),
                    child_fit=round(nodes[bc]['fit'],3),
                    dip=round(n['fit']-nodes[bp]['fit'],3), rise=round(nodes[bc]['fit']-n['fit'],3),
                    n_mut=n['n_mut'], step_p2v=added_p2v[0] if added_p2v else '?',
                    step_v2c=added_v2c[0] if added_v2c else '?'))
print(f"\n=== 破 greedy 谷节点(fit<所有parent 但有child超越所有parent)===", flush=True)
print(f"  谷节点: {valley_nodes}/{n_connected} 有边节点", flush=True)
print(f"  含义: greedy(每步选fit升)会停在parent不走这步(因降fit),但过谷后能到更高峰", flush=True)
for e in valley_examples[:5]:
    print(f"    parent {e['parent_fit']} --{e['step_p2v']}--> 谷 {e['valley_fit']}(dip{e['dip']}) --{e['step_v2c']}--> {e['child_fit']}(rise{e['rise']})", flush=True)

# 多跳路径深度: 从最低阶(k_min)到最高阶的最长 subset 链
n_mut_range=(min(n['n_mut'] for n in nodes), max(n['n_mut'] for n in nodes))
print(f"\n节点突变阶范围: {n_mut_range}", flush=True)
# 最长 DAG 路径(按突变阶递增)
import functools
order=sorted(range(len(nodes)),key=lambda i:nodes[i]['n_mut'])
longest=[1]*len(nodes)
for i in order:
    for c in adj_out[i]:
        if longest[i]+1>longest[c]: longest[c]=longest[i]+1
print(f"最长进化链(subset递增多跳路径): {max(longest)} 跳", flush=True)

import json
json.dump({'n_nodes':len(nodes),'n_edges':len(edges),'n_connected':n_connected,
           'frac_connected':n_connected/len(nodes),
           'deg_out_mean':float(np.mean(deg_out)),'deg_in_mean':float(np.mean(deg_in)),
           'df_up_frac':float((dfs>0).mean()),'df_down_frac':float((dfs<0).mean()),
           'valley_nodes':valley_nodes,'longest_chain':int(max(longest)),
           'n_mut_range':list(n_mut_range),'valley_examples':valley_examples},
          open(f'{D}/results/evolution_graph.json','w'),indent=2)
print("\nwrote results/evolution_graph.json", flush=True)
