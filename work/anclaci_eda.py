#!/usr/bin/env python3
"""AncLacI 数据 EDA + 整合: 三套训练数据
1) LGF景观: 1158 DBD序列 -> DNA结合fitness(log2_enrichment), 无配体(X0)
2) DMS: 1122 EcLacI单点突变 -> fitness
3) Hill(别构): 205克隆 × 8 IPTG剂量 -> GFP(剂量响应=±小分子维度)
目标: 训模型预测unseen突变体的DNA结合活性(+别构剂量响应)."""
import pandas as pd, numpy as np, os
ROOT="./work/anclaci/anclaci-ref"
OUT="./work/results/anclaci_model"; os.makedirs(OUT,exist_ok=True)

# ---- 1. LGF 景观 ----
lgf=pd.read_csv(f"{ROOT}/csvs/LGF_log2_enrichment.csv")
print("===== LGF 景观 =====")
print("shape:",lgf.shape," cols:",list(lgf.columns))
print("seq长度:",lgf.aa_seq.str.len().value_counts().head().to_dict())
print("fitness(avg_log2_enrichment): mean=%.2f std=%.2f min=%.2f max=%.2f"%(
    lgf.avg_log2_enrichment.mean(),lgf.avg_log2_enrichment.std(),
    lgf.avg_log2_enrichment.min(),lgf.avg_log2_enrichment.max()))
print("功能性(enrichment>0):%d  非功能(<=-5):%d"%(
    (lgf.avg_log2_enrichment>0).sum(),(lgf.avg_log2_enrichment<=-5).sum()))
lgf.to_csv(f"{OUT}/lgf_landscape.csv",index=False)

# ---- 2. DMS ----
dms=pd.read_csv(f"{ROOT}/csvs/DMS_log2_normalized_enrichment.csv")
print("\n===== DMS =====")
print("shape:",dms.shape)
print("seq长度:",dms.aa_seq.str.len().value_counts().head().to_dict())
print("fitness: mean=%.2f std=%.2f min=%.2f max=%.2f"%(
    dms.avg_log2_normalized_enrichment.mean(),dms.avg_log2_normalized_enrichment.std(),
    dms.avg_log2_normalized_enrichment.min(),dms.avg_log2_normalized_enrichment.max()))
# 解析突变 (lib_id 形如 A15P 等)
dms["mut"]=dms.lib_id
print("示例lib_id:",dms.lib_id.head(8).tolist())
dms.to_csv(f"{OUT}/dms_landscape.csv",index=False)

# ---- 3. Hill 别构(±IPTG) ----
hill=pd.read_csv(f"{ROOT}/Hill_fitting/sequenced_clonal_data.csv")
print("\n===== Hill 别构(±IPTG剂量响应) =====")
print("shape:",hill.shape," cols:",list(hill.columns)[:15])
print("IPTG剂量列:",[c for c in hill.columns if c.startswith('X')])
print("独特克隆(seq):",hill.seq.nunique()," 独特anc_id:",hill.anc_id.nunique())
# 每克隆多行(多孔重复) -> 取每克隆每剂量均值
doses=['X0','X0.1','X0.5','X1','X5','X25','X100','X1000']
h=hill.groupby(['seq','anc_id'])[doses].mean().reset_index()
print("聚合后(克隆×剂量):",h.shape)
print("剂量下GFP均值:",[f"{c}={h[c].mean():.0f}" for c in doses])
# 别构效应: X1000(饱和IPTG) - X0(无IPTG) = IPTG引起的DNA结合(抑制)变化
h["allosteric_shift"]=h["X1000"]-h["X0"]
print("allosteric_shift(IPTG饱-无): mean=%.0f std=%.0f min=%.0f max=%.0f"%(
    h.allosteric_shift.mean(),h.allosteric_shift.std(),h.allosteric_shift.min(),h.allosteric_shift.max()))
h.to_csv(f"{OUT}/hill_allosteric.csv",index=False)

print("\n===== 跨数据集: 序列长度一致? =====")
print("LGF len:",lgf.aa_seq.iloc[0].__len__(), " DMS len:",dms.aa_seq.iloc[0].__len__(),
      " Hill seq len:",h.seq.iloc[0].__len__())
print("产物:",os.listdir(OUT))
