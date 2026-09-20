#!/usr/bin/env python3
"""
逐靶点(胆酸)群体突变分布 + 与17位点负筛景观交叉.
输入: work/tsm_results/{CDCA,PCA}_fastq (纳米孔646bp扩增子pool)
参考: 用户贴的 alt-codon WT(612bp/203aa, 去末尾*). 蛋白级数突变(绕开公司背景错+密码子歧义).
每条read: 定位基因(ATGCAGAAAA锚, 双向) -> 612bp -> 翻译 -> 对WT蛋白每位点AA.
纳米孔噪声压: 频率化(每站点非WT AA的比例), 只保留 >=THR 的真信号.
"""
import gzip, os, sys
from collections import Counter, defaultdict
import numpy as np, pandas as pd
from Bio.Seq import Seq

WT_DNA="ATGCAGAAAAAATTGACTCGTTCTCAGCAGAAACACTTGGACATCATCAACGCTGCTAAAGAAGAATTCATCGAATTCGGTTTCTTGGCTGCTAACATGGACCGTATCACTTCTTCTGCTGAAGTATCTAAACGTACTTTGTACCGTCACTTCGAATCTAAAGAAGTATTGTTCGAATCTGTATTGACTATCATCAACGACTCTGTAAACGAATCTATCTCTTACCACTTCGACCCGAACAAATCTACTGAAGAACAGTTGACTGAAATCGCTTACAAAGAAATCGACGTATTGTACAAAACTTACGGTATCGCTTTGGCTCGTACTATCGTAATGGAATTCTTGCGTCAGCCGGAAATGGCTAAAACTTTGATCCAGAACATCTACTCTATCCGTGCTATCACTCAGTGGTTCCGTTCTGCTATCGAAGCTAAACGTTTGAAAGACGCTGACCCGAAATTGATGACTGACGTATACGTATCTTTGTTCCAGGGTTTGTTCTTCTGGCCGCAGGTAATGCACTTGGACTTGGAACCGCACGGTGAAGAATTGTCTCAGAAAATCGAAACTTTGACTACTATCTTCTTGCAGTCTTACGGTGTAGCTGAATAA"
CODON={ 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L','ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P','ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E','TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
def tr(d): return "".join(CODON.get(d[i:i+3],'X') for i in range(0,len(d)-2,3))
def rc(s): return str(Seq(s).reverse_complement())
WT_PROT=tr(WT_DNA).rstrip('*')  # 203aa
WL=len(WT_PROT)
SITES=[66,69,75,94,102,107,111,125,128,129,133,159,161,163,164,169,170]
wt_at={p:WT_PROT[p-1] for p in SITES}
ANCHOR="ATGCAGAAAA"   # 两版WT共有前缀, 用于定位基因起点
ANCHOR_RC=rc(ANCHOR)
THR=10.0  # 非WT频率(%)阈值, 压纳米孔噪声

TSM="./work/tsm_results"
OUT="./work/results/per_target"; os.makedirs(OUT,exist_ok=True)

def iter_reads(path):
    op=gzip.open if path.endswith('.gz') else open
    fh=op(path,'rt')
    h=fh.readline()
    while h:
        s=fh.readline().strip(); fh.readline(); fh.readline()
        yield s
        h=fh.readline()

def call_read(read):
    """返回 (gene612, prot203, ok) 或 None."""
    # 定位锚点(双向)
    if ANCHOR in read:
        i=read.find(ANCHOR); gene=read[i:i+612]
    elif ANCHOR_RC in read:
        r=rc(read); i=r.find(ANCHOR); gene=r[i:i+612]
    else:
        return None
    if len(gene)!=612: return None
    prot=tr(gene).rstrip('*')
    if len(prot)!=WL or '*' in prot:  # 移码/截断/终止 -> 丢
        return None
    return gene,prot

targets=["CDCA","PCA"]
summary={}
for tgt in targets:
    path=os.path.join(TSM,f"{tgt}_fastq")
    n=0; ok=0; noanchor=0; badframe=0
    site_aa=defaultdict(Counter)       # pos -> Counter(aa)
    geno_cnt=Counter()                 # 17site genotype str -> count
    full_mut_cnt=Counter()            # full prot mutation-string -> count
    cons_gene=Counter()              # for consensus (per-nt, 略)
    for s in iter_reads(path):
        n+=1
        r=call_read(s)
        if r is None:
            if ANCHOR in s or ANCHOR_RC in s: badframe+=1
            else: noanchor+=1
            continue
        gene,prot=r; ok+=1
        # 每位点AA
        for i in range(WL): site_aa[i+1][prot[i]]+=1
        # 17位点基因型
        g="".join(prot[p-1] for p in SITES)
        geno_cnt[g]+=1
        # 全突变串
        muts=tuple(f"{WT_PROT[i]}{i+1}{prot[i]}" for i in range(WL) if prot[i]!=WT_PROT[i])
        full_mut_cnt[muts]+=1
    summary[tgt]=dict(n=n,ok=ok,noanchor=noanchor,badframe=badframe)
    print(f"\n===== {tgt}  reads={n:,}  有效翻译={ok:,}  无锚={noanchor:,}  移帧/含*={badframe:,} =====")
    # 位点级频率: 非WT AA
    rows=[]
    for pos in range(1,WL+1):
        c=site_aa[pos]; tot=sum(c.values())
        wt_aa=WT_PROT[pos-1]
        wt_c=c.get(wt_aa,0)
        non_wt_frac=(tot-wt_c)/tot*100 if tot else 0
        top_alt=c.most_common(3)
        # 去掉WT自己, 取最常alt
        alts=[(a,n) for a,n in top_alt if a!=wt_aa]
        rows.append(dict(pos=pos, wt=wt_aa, depth=tot, nonwt_pct=non_wt_frac,
                         top_alt=alts[0] if alts else None,
                         all_alt=";".join(f"{a}:{cnt}({cnt/tot*100:.1f}%)" for a,cnt in c.items() if a!=wt_aa and cnt/tot*100>=1)))
    posdf=pd.DataFrame(rows)
    posdf.to_csv(f"{OUT}/{tgt}_position_freq.csv",index=False)
    # 显著位点(非WT>=THR%)
    sig=posdf[posdf.nonwt_pct>=THR].sort_values("nonwt_pct",ascending=False)
    print(f"  显著突变位点(非WT>={THR}%): {len(sig)}")
    for _,r in sig.head(25).iterrows():
        insite= int(r.pos) in SITES
        print(f"    aa{int(r.pos):>3}({r.wt}{'*' if insite else ' '})  nonWT={r.nonwt_pct:5.1f}%  alt={r.top_alt}  [design-site:{insite}]")
    # 17位点基因型分布
    print(f"  17位点基因型数: {len(geno_cnt)}, 前8:")
    for g,c in geno_cnt.most_common(8):
        muts=";".join(f"{wt_at[SITES[i]]}{SITES[i]}{g[i]}" for i in range(17) if g[i]!=wt_at[SITES[i]]) or "WT"
        print(f"    {c:>6}  {muts}")
    geno_cnt_df=pd.DataFrame([{"genotype":g,"count":c,"pct":c/ok*100} for g,c in geno_cnt.most_common()])
    geno_cnt_df.to_csv(f"{OUT}/{tgt}_17site_genotype_dist.csv",index=False)
    # 顶级完整突变组合
    print(f"  顶级完整突变组合 前8:")
    for m,c in full_mut_cnt.most_common(8):
        print(f"    {c:>6}  {';'.join(m) if m else 'WT(0 mut)'}")
    full_mut_df=pd.DataFrame([{"mutations":";".join(m) if m else "WT","count":c,"pct":c/ok*100} for m,c in full_mut_cnt.most_common(50)])
    full_mut_df.to_csv(f"{OUT}/{tgt}_full_mutation_dist_top50.csv",index=False)

pd.DataFrame(summary).T.to_csv(f"{OUT}/read_summary.csv")
print("\n产物:",os.listdir(OUT))
