#!/usr/bin/env python3
"""推荐器 v5: 对接感知重排 (self-ESMFold pocket-aware)。
与 v4 的差别: 口袋距离改用【自身 ESMFold 结构】的 d_Q164 (self_site_distances.json),
替代 v4 的 3CDL 同源 d_pocket。只对 17 设计位点重排 (自身结构仅覆盖设计区)。
打分 = 0.45*esm_tol + 0.35*pocket_prox + 0.20*hth_allo
  esm_tol   : ESM2-3B 该 AA 容忍度 (位点内 min-max 归一)
  pocket_prox: 1/(1+d_Q164) 归一 —— 越靠近配体口袋越可能调特异性
  hth_allo   : 1/(1+|d_HTH - median|) ... 用 d_HTH 归一 —— 别构到 DNA 结合面
输出 CDCA / 7K-LCA 两个靶点的重排候选 (排除 Q164/R116 已固定)。
"""
import pandas as pd, numpy as np, os, json, ast

BASE = os.path.dirname(os.path.abspath(__file__))
R6   = os.path.join(BASE, "results/per_target_six")
RSST = os.path.join(BASE, "results/structure")
STR  = os.path.join(BASE, "structure/esmfold")
OUT  = os.path.join(BASE, "results/recommender")
AAS  = list("ACDEFGHIKLMNPQRSTVWY")
EXCLUDE_POS = {164, 116}

# --- 自身 ESMFold 距离 ---
dist = json.load(open(os.path.join(STR, "self_site_distances.json")))
dQ   = {d["pos"]: d["d_Q164"] for d in dist}
dHTH = {d["pos"]: d["d_HTH"]  for d in dist}
wt_at= {d["pos"]: d["aa"]     for d in dist}
SITES = sorted(dQ.keys())

# pocket proximity (越近越高)
pv = np.array([1.0/(1.0+dQ[p]) for p in SITES])
pv = (pv - pv.min())/(pv.max()-pv.min()+1e-9)
pocket_prox = dict(zip(SITES, pv))
# HTH allostery: 越靠近 HTH (d_HTH 小) 别构耦合越直接
hv = np.array([1.0/(1.0+dHTH[p]) for p in SITES])
hv = (hv - hv.min())/(hv.max()-hv.min()+1e-9)
hth_allo = dict(zip(SITES, hv))

# --- ESM2 全位点打分 ---
esm = pd.read_csv(os.path.join(RSST, "esm2_t36_3B_zeroshot_allpos.csv"))
esm_map = {(int(r.pos), r.alt): r.score for _, r in esm.iterrows()}
def esm_tol_at(p):
    arr = np.array([esm_map.get((p,a), -10.0) for a in AAS])
    n = (arr - arr.min())/(arr.max()-arr.min()+1e-9)
    return dict(zip(AAS, n))

# --- 各靶点实测突变 (标注 borrowed) ---
def observed(t):
    pf = pd.read_csv(os.path.join(R6, f"{t}_position_freq.csv"))
    obs = {}
    for _, r in pf.iterrows():
        ta = r.top_alt
        if isinstance(ta, str):
            try: ta = ast.literal_eval(ta)
            except: ta = None
        if isinstance(ta, tuple) and len(ta) >= 2 and r.nonwt_pct >= 5:
            obs[(int(r.pos), ta[0])] = float(r.nonwt_pct)
    return obs

# 全靶点观察(供 borrow 标注)
ALLT = ["CDCA","PCA","7K-LCA","23K-CDCA","LCA","LLDCA"]
obs_all = {t: observed(t) for t in ALLT}

# 靶点特异先验: 该靶点(及其结构分组同伴)实测富集过的位点 -> 邻域(序列±3)加权
GROUP = {  # 结构分组: 7位修饰组 / C23侧链组 (见 RESULTS_annotated §7)
    "7K-LCA": ["7K-LCA","LCA","LLDCA"],   # 7位修饰组
    "CDCA":   ["CDCA"],                    # 基准
}
def target_site_prior(target):
    """靶点(+同组)实测显著突变的位点, 及其序列邻域, 作为特异先验."""
    hot = set()
    for t in GROUP.get(target, [target]):
        for (p,a),pct in obs_all[t].items():
            if pct >= 5: hot.add(p)
    prior = {}
    for p in SITES:
        prior[p] = 1.0 if p in hot else (0.5 if any(abs(p-h)<=3 for h in hot) else 0.0)
    return prior

def rerank(target):
    rows = []
    own = obs_all[target]
    tprior = target_site_prior(target)
    for p in SITES:
        if p in EXCLUDE_POS: continue
        et = esm_tol_at(p)
        for a in AAS:
            if a == wt_at[p]: continue
            # 重平衡: 口袋近邻+靶点特异先验为主, ESM2 只作"可行性门"
            score = (0.30*pocket_prox[p] + 0.15*hth_allo[p]
                     + 0.35*tprior[p] + 0.20*et[a])
            borrowed = [t for t in ALLT if t!=target and (p,a) in obs_all[t]]
            rows.append({
                "target": target, "pos": p, "wt": wt_at[p], "alt": a,
                "mutation": f"{wt_at[p]}{p}{a}",
                "score": round(score,4),
                "esm_tol": round(et[a],3),
                "d_Q164": dQ[p], "d_HTH": dHTH[p],
                "pocket_prox": round(pocket_prox[p],3),
                "site_prior": tprior[p],
                "in_target_already": round(own.get((p,a),0.0),1),
                "borrowed_from": ",".join(borrowed) if borrowed else "-",
            })
    df = pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)
    return df

for t in ["CDCA","7K-LCA"]:
    df = rerank(t)
    df.to_csv(os.path.join(OUT, f"{t}_rec_v5_dockaware.csv"), index=False)
    print(f"\n===== {t} v5 对接感知 top15 (排除 Q164R/R116Q) =====")
    print(df.head(15)[["mutation","score","esm_tol","d_Q164","d_HTH","site_prior","in_target_already","borrowed_from"]].to_string(index=False))
print("\n[done] 写出 CDCA_rec_v5_dockaware.csv / 7K-LCA_rec_v5_dockaware.csv")
