#!/usr/bin/env python3
"""别构(±小分子)维度模型: 用Hill数据(54克隆×8 IPTG剂量)预测
   unseen突变体的剂量响应曲线(GFP vs [IPTG]).
特征: ESM2 t30_150M mean-pool嵌入. 模型: Ridge多输出(每剂量一输出) + Hill参数拟合.
验证: LOO-CV. 这是v2"小分子有/无DNA结合变化"的核心模型. """
import pandas as pd, numpy as np, os, torch
from esm.pretrained import esm2_t30_150M_UR50D
from sklearn.linear_model import Ridge
from sklearn.model_selection import LeaveOneOut
from scipy.stats import spearmanr
ROOT="./work/results/anclaci_model"
OUT=ROOT
hill=pd.read_csv(f"{ROOT}/hill_allosteric.csv")  # 已聚合 seq×dose
doses=['X0','X0.1','X0.5','X1','X5','X25','X100','X1000']
iptg=[0,0.1,0.5,1,5,25,100,1000]
print("Hill聚合数据:",hill.shape,"  独特克隆:",hill.seq.nunique())
Y=hill[doses].values  # [54,8] GFP per dose
import re
seqs=[re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]','',s) for s in hill.seq.tolist()]
print("清洗后长度:",{len(s) for s in seqs})

# ESM2 嵌入(54条,快)
print("嵌入..."); model,alphabet=esm2_t30_150M_UR50D(); model.eval(); bc=alphabet.get_batch_converter()
embs=[]
with torch.no_grad():
    for i in range(0,len(seqs),8):
        chunk=[("s",s) for s in seqs[i:i+8]]; _,_,toks=bc(chunk)
        out=model(toks,repr_layers=[30],return_contacts=False); rep=out["representations"][30]
        mask=(toks!=alphabet.padding_idx).float().unsqueeze(-1)
        for j in range(rep.shape[0]): embs.append(((rep[j]*mask[j]).sum(0)/mask[j].sum(0)).numpy())
X=np.array(embs); print("emb:",X.shape)
np.save(f"{OUT}/esm2_emb_hill.npy",X)

# 多输出Ridge LOO-CV: 预测每剂量GFP
loo=LeaveOneOut(); pred=np.zeros_like(Y)
for tr,te in loo.split(X):
    m=Ridge(alpha=1.0).fit(X[tr],Y[tr]); pred[te]=m.predict(X[te])
# 每克隆预测曲线的Spearman(剂量-响应形状) + 数值MSE
sp_shape=[]; per_dose_err=[]
for i in range(len(Y)):
    s=spearmanr(pred[i],Y[i]).correlation; sp_shape.append(s if not np.isnan(s) else 0)
    per_dose_err.append(np.abs(pred[i]-Y[i]))
sp_shape=np.array(sp_shape); per_dose_err=np.array(per_dose_err)
print("\n===== LOO-CV: 剂量响应曲线预测 =====")
print(f"曲线形状Spearman(pred vs true, per clone): mean={sp_shape.mean():.2f} std={sp_shape.std():.2f}")
print(f"每剂量绝对误差: mean={per_dose_err.mean():.0f} GFP  (GFP范围 {Y.min():.0f}~{Y.max():.0f})")
# 别构shift(IPTG饱-无)预测
true_shift=Y[:,-1]-Y[:,0]; pred_shift=pred[:,-1]-pred[:,0]
print(f"别构shift预测: Spearman={spearmanr(pred_shift,true_shift).correlation:.2f}  (true {true_shift.min():.0f}~{true_shift.max():.0f})")
# Hill参数拟合(对true和pred各拟EC50/dynamic range)
from scipy.optimize import curve_fit
def hill_eq(x,bot,top,ec50,n):
    return bot+(top-bot)*x**n/(ec50**n+x**n)
def fit_params(y):
    try:
        p,_=curve_fit(hill_eq,iptg,y,p0=[y.min(),y.max(),10,1],maxfev=5000,bounds=([0,0,0.05,0.5],[1e6,1e6,2000,5]))
        return p
    except: return [np.nan]*4
true_p=np.array([fit_params(Y[i]) for i in range(len(Y))])
pred_p=np.array([fit_params(pred[i]) for i in range(len(Y))])
for k,nm in enumerate(["bot","top","EC50","n"]):
    m=~np.isnan(true_p[:,k])&~np.isnan(pred_p[:,k])
    if m.sum()>3: print(f"  Hill {nm}: pred-vs-true Spearman={spearmanr(pred_p[m,k],true_p[m,k]).correlation:+.2f} (n={m.sum()})")
np.save(f"{OUT}/hill_pred.npy",pred); np.save(f"{OUT}/hill_true.npy",Y)
print("done. 产物:",os.listdir(OUT))
