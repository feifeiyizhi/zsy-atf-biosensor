#!/usr/bin/env python3
"""基线模型: 预测unseen突变体DNA结合活性(log2_enrichment).
特征: ESM2 t30_150M mean-pool嵌入(已缓存). 模型: Ridge/GBM/MLP.
数据: LGF景观(1158) + DMS单点(1122). 验证: 随机CV + 跨集迁移(LGF训->DMS测)."""
import pandas as pd, numpy as np, os, torch
from esm.pretrained import esm2_t30_150M_UR50D
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.linear_model import Ridge
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from scipy.stats import spearmanr
ROOT="./work/results/anclaci_model"
OUT=ROOT; os.makedirs(OUT,exist_ok=True)

lgf=pd.read_csv(f"{ROOT}/lgf_landscape.csv")
dms=pd.read_csv(f"{ROOT}/dms_landscape.csv")
y_lgf=lgf.avg_log2_enrichment.values
y_dms=dms.avg_log2_normalized_enrichment.values

print("加载 ESM2 t30_150M(缓存)...")
model,alphabet=esm2_t30_150M_UR50D()
model.eval()
bc=alphabet.get_batch_converter()

def embed(seqs,batch=8):
    embs=[]
    with torch.no_grad():
        for i in range(0,len(seqs),batch):
            chunk=[("s",s) for s in seqs[i:i+batch]]
            _,_,toks=bc(chunk)
            out=model(toks,repr_layers=[30],return_contacts=False)
            rep=out["representations"][30]  # [B,L,D]
            mask=(toks != alphabet.padding_idx)
            for j in range(rep.shape[0]):
                m=mask[j].float().unsqueeze(-1)
                embs.append((rep[j]*m).sum(0)/m.sum(0))  # mean-pool
    return torch.stack(embs).numpy()
print("嵌入 LGF..."); Elgf=embed(lgf.aa_seq.tolist())
print("嵌入 DMS..."); Edms=embed(dms.aa_seq.tolist())
print("LGF emb:",Elgf.shape," DMS emb:",Edms.shape)
np.save(f"{OUT}/esm2_emb_lgf.npy",Elgf); np.save(f"{OUT}/esm2_emb_dms.npy",Edms)

def eval_model(name,mdl,X,y,nf=5):
    kf=KFold(n_splits=nf,shuffle=True,random_state=42)
    pred=cross_val_predict(mdl,X,y,cv=kf)
    sp=spearmanr(pred,y).correlation
    mse=np.mean((pred-y)**2)
    # 功能性>0 的ranking (top-k precision)
    truth_func=y>0
    topk=np.argsort(-pred)[:max(5,int(truth_func.sum()*0.5))]
    prec=len([i for i in topk if truth_func[i]])/len(topk)
    print(f"  {name:18} CV Spearman={sp:+.3f} MSE={mse:.2f} func-topk-precision={prec:.2f}")
    return pred,sp

print("\n===== LGF 内部5折CV(预测unseen景观变体) =====")
Xs=StandardScaler().fit_transform(Elgf)
for m in [Ridge(alpha=1.0), GradientBoostingRegressor(n_estimators=200,max_depth=3,random_state=42),
          make_pipeline(StandardScaler(),MLPRegressor(hidden_layer_sizes=(128,64),max_iter=300,random_state=42))]:
    eval_model(type(m).__name__ if not hasattr(m,'named_steps') else "MLP",m,Xs,y_lgf)

print("\n===== 跨集迁移: LGF训 -> DMS测(单点突变) =====")
for nm,mdl in [("Ridge",Ridge(alpha=1.0)),("GBM",GradientBoostingRegressor(n_estimators=200,max_depth=3,random_state=42))]:
    mdl.fit(Elgf,y_lgf)
    p=mdl.predict(Edms)
    sp=spearmanr(p,y_dms).correlation
    # DMS里WT应~0, 功能性>0的少
    print(f"  {nm}: LGF->DMS Spearman={sp:+.3f} (DMS y range {y_dms.min():.1f}~{y_dms.max():.1f})")

print("\n===== DMS 内部5折CV =====")
for nm,mdl in [("Ridge",Ridge(alpha=1.0)),("GBM",GradientBoostingRegressor(n_estimators=200,max_depth=3,random_state=42))]:
    eval_model(nm,mdl,Edms,y_dms)
print("done")
