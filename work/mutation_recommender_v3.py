#!/usr/bin/env python3
"""突变推荐v3: 在Q164R(通用开关)+R116Q(骨架)之外, 预测'再加哪个突变增强某胆酸响应'.
聚焦次级热点位置(LBD调谐位), 主信号=借鉴其他靶点在该位置已验证的AA."""
import pandas as pd, numpy as np, os, ast
R6="./work/results/per_target_six"
RSST="./work/results/structure"
OUT="./work/results/recommender"; os.makedirs(OUT,exist_ok=True)
TARGETS=["CDCA","PCA","7K-LCA","23K-CDCA","LCA","LLDCA"]
AAS=list("ACDEFGHIKLMNPQRSTVWY")
WT="MQKKLTRSQQKHLDIINAAKEEFIEFGFLAANMDRITSSAEVSKRTLYRHFESKEVLFESVLTIINDSVNESISYHFDPNKSTEEQLTEIAYKEIDVLYKTYGIALARTIVMEFLRQPEMAKTLIQNIYSIRAITQWFRSAIEAKRLKDADPKLMTDVYVSLFQGLFFWPQVMHLDLEPHGEELSQKIETLTTIFLQSYGVAE"
EXCLUDE_POS={116,164}  # Q164R/R116Q 已固定, 不再推荐该位置

# 1. 6靶点全突变(>=5%)
mut_by_target={t:{} for t in TARGETS}  # t -> {(pos,alt):freq}
pos_alts={}  # (pos,alt) -> {target:freq} 跨靶点
for t in TARGETS:
    pf=pd.read_csv(f"{R6}/{t}_position_freq.csv")
    for _,r in pf.iterrows():
        ta=r.top_alt
        if isinstance(ta,str):
            try: ta=ast.literal_eval(ta)
            except: ta=None
        if isinstance(ta,tuple) and len(ta)>=2 and r.nonwt_pct>=5:
            mut_by_target[t][(int(r.pos),ta[0])]=r.nonwt_pct
            pos_alts.setdefault((int(r.pos),ta[0]),{})[t]=r.nonwt_pct

# 次级热点位置(>=1靶点有>=5%突变, 排除Q164/R116)
hotspots=sorted(set(p for (p,a) in pos_alts if p not in EXCLUDE_POS))
print("次级响应调谐热点位置:",hotspots)
print("各位置已验证(target:freq):")
for p in hotspots:
    print(f"  aa{p}({WT[p-1]}): ", end="")
    for a in sorted(set(a2 for (p2,a2) in pos_alts if p2==p)):
        tg=pos_alts[(p,a)]
        print(f"{a}<-{'+'.join(f'{t}:{fr:.0f}%' for t,fr in tg.items())}",end="  ")
    print()

# 2. ESM2全位点
esm=pd.read_csv(f"{RSST}/esm2_zeroshot_allpos.csv")
esm_map={(int(r.pos),r.alt):r.score for _,r in esm.iterrows()}

# 3. 别构相关(3D距口袋, 17位点; 其他用热点计数)
dist=pd.read_csv(f"{RSST}/site_3d_distance.csv")
dpock={int(r.wt_pos):float(r.d_pocket) for _,r in dist.iterrows() if pd.notna(r.d_pocket)}
def allo(p):
    if p in dpock: return 1.0/(1.0+dpock[p])
    # 非设计位点: 在热点=有结构意义, 给0.5基准
    return 0.5 if p in hotspots else 0.0

def esm_norm_at(p):
    sub=[esm_map.get((p,a),-10) for a in AAS]
    return dict(zip(AAS,(np.array(sub)-min(sub))/(max(sub)-min(sub)+1e-9)))

print("\n"+"="*70); print("每靶点 下一步加什么突变增强响应(在Q164R/R116Q之外)"); print("="*70)
allrecs=[]
for t in TARGETS:
    known=set(mut_by_target[t].keys())  # 该靶点已有的(pos,alt)
    recs=[]
    for p in hotspots:
        if p in EXCLUDE_POS: continue
        wt_aa=WT[p-1]
        en=esm_norm_at(p)
        # 该位置其他靶点验证过的AA(借鉴)
        borrowed=set(a for (p2,a) in pos_alts if p2==p)
        for alt in AAS:
            if alt==wt_aa: continue
            if (p,alt) in known: continue  # 该靶点已有
            borrow_freq=max(pos_alts.get((p,alt),{}).get(ot,0) for ot in TARGETS if ot!=t)/100.0 if (p,alt) in pos_alts else 0
            # 综合分: 借鉴(0.45,主) + ESM2容忍(0.3) + 别构(0.25)
            score=0.45*borrow_freq+0.30*en[alt]+0.25*allo(p)
            recs.append(dict(target=t,pos=p,mutation=f"{wt_aa}{p}{alt}",score=round(score,3),
                borrowed_from=(f"{borrow_freq*100:.0f}%" if borrow_freq>0 else "-"),
                esm=round(esm_map.get((p,alt),-10),2),allosteric=round(allo(p),2)))
    df=pd.DataFrame(recs).sort_values("score",ascending=False)
    df.to_csv(f"{OUT}/{t}_rec_v3.csv",index=False)
    known_str=sorted(f"{WT[p-1]}{p}{a}({fr:.0f}%)" for (p,a),fr in mut_by_target[t].items() if p not in EXCLUDE_POS)
    print(f"\n### {t}  (已固定次级: {known_str}) ###")
    print(df.head(10)[["mutation","score","borrowed_from","esm","allosteric"]].to_string(index=False))
    allrecs.append(df.head(10))
pd.concat(allrecs).to_csv(f"{OUT}/all_rec_v3_top10.csv",index=False)
print("\n产物:",os.listdir(OUT))
