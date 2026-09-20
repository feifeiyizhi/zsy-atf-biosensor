#!/usr/bin/env python3
"""ESM-2 zero-shot 突变效应打分: 对17设计位点 × 20AA, 预测替换对序列似然的扰动
(高=模型认为可容忍/有利, 低=有害). 一次性前向(掩码语言模型, 每位输出全部AA logits).
与338lib真实负筛适应度交叉验证(突变敏感的AI先验 vs 真实DNA结合选择)."""
import torch, numpy as np, pandas as pd, os
from esm.pretrained import esm2_t30_150M_UR50D as esm2_xxx
import esm.pretrained as ep

WT="MQKKLTRSQQKHLDIINAAKEEFIEFGFLAANMDRITSSAEVSKRTLYRHFESKEVLFESVLTIINDSVNESISYHFDPNKSTEEQLTEIAYKEIDVLYKTYGIALARTIVMEFLRQPEMAKTLIQNIYSIRAITQWFRSAIEAKRLKDADPKLMTDVYVSLFQGLFFWPQVMHLDLEPHGEELSQKIETLTTIFLQSYGVAE"
SITES=[66,69,75,94,102,107,111,125,128,129,133,159,161,163,164,169,170]
AAS=list("ACDEFGHIKLMNPQRSTVWY")

OUT="./work/results/structure"; os.makedirs(OUT,exist_ok=True)

print("加载 ESM2 t33_650M ..."); import sys; sys.stdout.flush()
model,alphabet=esm2_xxx(); batch_converter=alphabet.get_batch_converter()
model.eval()
data=[("WT",WT)]
_,_,tokens=batch_converter(data)
toks=tokens  # ESM2: batch_converter 直接返回张量
print("tokens shape:",toks.shape)

with torch.no_grad():
    out=model(toks, repr_layers=[], return_contacts=False)
logits=out["logits"][0]  # [L+2, vocab]
logp=torch.log_softmax(logits,dim=-1)
# AA token id
aa_idx={aa:alphabet.get_idx(aa) for aa in AAS}

# 17位点 × 20AA: score = logp(aa) - logp(wt_aa)  (位置: token index = i, 因为BOS在0, 残基1在token1)
rows=[]
for p in SITES:
    ti=p  # 残基p(1-indexed) 对应 token index p (BOS在0)
    wt_aa=WT[p-1]
    wt_lp=logp[ti,aa_idx[wt_aa]].item()
    for aa in AAS:
        s=logp[ti,aa_idx[aa]].item()-wt_lp
        rows.append(dict(site=p,wt=wt_aa,alt=aa,score=s,
                         is_wt=(aa==wt_aa)))
df=pd.DataFrame(rows)
df.to_csv(f"{OUT}/esm2_zeroshot_17site.csv",index=False)
print(f"\n=== ESM2 zero-shot: 各位点 最佳非WT替换 (模型预测最可容忍) ===")
for p in SITES:
    sub=df[df.site==p].sort_values("score",ascending=False)
    top=sub[~sub.is_wt].head(3)
    wt_aa=WT[p-1]
    print(f"  aa{p}({wt_aa}): "+"  ".join(f"{r.alt}:{r.score:+.2f}" for _,r in top.iterrows()))

# === 与真实负筛交叉: 338lib 各位点 alt 的平均适应度 ===
g=pd.read_csv("./work/results/genotype_17site_landscape.csv")
spec=pd.read_csv("./work/results/mutated_position_spectrum.csv")
wt_at={int(r.aa_pos):r.wt_aa for _,r in spec.iterrows()}; wt_at={p:wt_at[p] for p in SITES}
sitecols=[f"p{p}_{wt_at[p]}" for p in SITES]
rel=g[g.reliable].copy()
print("\n=== ESM2预测 vs 真实负筛(每位点 alt 的平均适应度) 相关 ===")
from scipy.stats import spearmanr
for p in SITES:
    col=f"p{p}_{wt_at[p]}"
    # 该位点携带 alt aa 的可靠变体平均 fitness
    real_map=rel.groupby(col)["fitness_slope"].mean().to_dict()
    rows2=[]
    for aa in AAS:
        if aa==wt_at[p]: continue
        # real fitness for variants carrying aa at p
        if aa in real_map:
            esm=df[(df.site==p)&(df.alt==aa)].score.values[0]
            rows2.append(dict(aa=aa,esm=esm,real=real_map[aa],n=rel[rel[col]==aa].shape[0]))
    if len(rows2)>=3:
        d=pd.DataFrame(rows2)
        rho,pv=spearmanr(d.esm,d.real)
        print(f"  aa{p}({wt_at[p]}): rho={rho:+.2f} p={pv:.2f} n={len(d)} (alts: "+", ".join(f"{r.aa}({r.n})" for _,r in d.iterrows())+")")
d.to_csv(f"{OUT}/esm2_vs_realsel_per_site.csv",index=False) if len(d) else None
print("\n产物:",os.listdir(OUT))
