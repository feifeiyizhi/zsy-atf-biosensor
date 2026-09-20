#!/usr/bin/env python3
"""ProDy ANM 弹性网络: 在3CDL同源结构上做别构构象模式分析.
- 最慢模式波动谱 + 找"别构耦合"残基(口袋<->HTH).
- 扰动响应: 在配体口袋残基(3CDL 164=Q164)施加扰动, 看HTH(3CDL 34=α3)的响应.
- 把17位点投到波动谱上(谁最动态/谁刚性).
结构侧突变敏感性的轻量代理(ESMFold缺openfold时的兜底)."""
import numpy as np, pandas as pd, os, json
from prody import parsePDB, ANM, calcPerturbResponse, extendModel, LOGGER
LOGGER.verbosity='none'
ROOT="./work"
OUT=f"{ROOT}/results/structure"; os.makedirs(OUT,exist_ok=True)

# 17位点 -> 3CDL 残基号 映射
mp=json.load(open(f"{ROOT}/structure/wt_to_3cdl_map.json"))
mp={int(k):v for k,v in mp.items()}
SITES=[66,69,75,94,102,107,111,125,128,129,133,159,161,163,164,169,170]
# 关键: Q164=3CDL164(口袋), α3=3CDL34(HTH识别螺旋), R116=3CDL115
POCKET=mp[164]; HTHELIX=mp[35]; R116=mp[116]

struct=parsePDB(f"{ROOT}/structure/pdb/3CDL.pdb")
# 选链A的CA
chains=set(struct.getChids())
print("chains:",chains)
ca=struct.select('protein and name CA and chain A')
if ca is None or len(ca)<50: ca=struct.select('protein and name CA')  # 兜底
print("CA atoms:",len(ca)," chainids:",set(ca.getChids()))
resnums=ca.getResnums()
# 建3CDL残基号->CA index 映射
r2i={rn:i for i,rn in enumerate(resnums)}

# ANM
anm=ANM('3CDL')
anm.buildHessian(ca.getCoords(),cutoff=15.0)
anm.calcModes(n_modes=None)  # 全部
n=len(ca)
print("ANM built, n modes:",anm.numModes())

# 波动谱 (最慢20模式)
fluct=anm.getVariances()[:20].sum(axis=0) if False else None
# ProDy: 用 calcSqFlucts over slowest 20 modes
from prody import calcSqFlucts
slow=anm[:20]
sqf=calcSqFlucts(slow)  # per-residue
print("sqf len:",len(sqf)," min/max:",round(sqf.min(),3),round(sqf.max(),3))

# 17位点波动
rows=[]
for p in SITES:
    pr=mp.get(p)
    if pr is None or pr not in r2i: rows.append(dict(wt_pos=p,cdl_res=None,sqf=None)); continue
    rows.append(dict(wt_pos=p,cdl_res=pr,sqf=round(float(sqf[r2i[pr]]),4)))
pd.DataFrame(rows).to_csv(f"{OUT}/anm_17site_fluctuation.csv",index=False)
print("\n17位点 最慢20模式 波动(高=灵活/低=刚性):")
for r in rows: print(f"  WT aa{r['wt_pos']} (3CDL {r['cdl_res']}): sqf={r['sqf']}")

# 扰动响应: calcPerturbResponse 返回 PR[n,n] (PR[i,j]=扰动i对j的响应)
pi=r2i[POCKET]  # 口袋残基在CA中的index
PR=np.array(calcPerturbResponse(anm)[0])  # [n,n]
print("PR shape:",PR.shape)
r_pocket=PR[pi,:]   # 扰动口袋 -> 各残基响应
hi=r2i[HTHELIX]
r_hth=PR[hi,:]      # 扰动HTH -> 各残基响应
def getidx(pr): return r2i.get(pr)
print(f"\n扰动配体口袋(3CDL {POCKET}=Q164, idx{pi}) -> 各处响应:")
if HTHELIX in r2i: print(f"  HTH α3(3CDL {HTHELIX})响应: {r_pocket[r2i[HTHELIX]]:.3f}")
if R116 in r2i: print(f"  R116(3CDL {R116})响应: {r_pocket[r2i[R116]]:.3f}")
print("  17位点对口袋扰动的响应:")
for p in SITES:
    pr2=mp.get(p)
    if pr2 and pr2 in r2i: print(f"    aa{p}(3CDL {pr2}): {r_pocket[r2i[pr2]]:.3f}")
if POCKET in r2i:
    print(f"\n扰动HTH(3CDL {HTHELIX}) -> 口袋Q164(3CDL {POCKET})响应: {r_hth[r2i[POCKET]]:.3f}")
np.save(f"{OUT}/anm_sqf.npy",sqf); np.save(f"{OUT}/anm_PR.npy",PR)
np.save(f"{OUT}/anm_pocket_perturb_resp.npy",r_pocket)
print("\n产物:",os.listdir(OUT))
