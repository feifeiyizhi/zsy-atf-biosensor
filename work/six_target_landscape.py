#!/usr/bin/env python3
"""6靶点(胆酸)群体突变分布 + 共享/特异景观.
CDCA+PCA(work/tsm_results) + 7K-LCA+23K-CDCA+LCA+LLDCA(work/tsm2).
每靶点纳米孔fastq -> 蛋白级突变(ALT密码子WT) -> 频率化 -> 6靶点共享/特异 + 负筛景观交叉."""
import gzip, os, json
from collections import Counter, defaultdict
import numpy as np, pandas as pd
from Bio.Seq import Seq

WT_DNA="ATGCAGAAAAAATTGACTCGTTCTCAGCAGAAACACTTGGACATCATCAACGCTGCTAAAGAAGAATTCATCGAATTCGGTTTCTTGGCTGCTAACATGGACCGTATCACTTCTTCTGCTGAAGTATCTAAACGTACTTTGTACCGTCACTTCGAATCTAAAGAAGTATTGTTCGAATCTGTATTGACTATCATCAACGACTCTGTAAACGAATCTATCTCTTACCACTTCGACCCGAACAAATCTACTGAAGAACAGTTGACTGAAATCGCTTACAAAGAAATCGACGTATTGTACAAAACTTACGGTATCGCTTTGGCTCGTACTATCGTAATGGAATTCTTGCGTCAGCCGGAAATGGCTAAAACTTTGATCCAGAACATCTACTCTATCCGTGCTATCACTCAGTGGTTCCGTTCTGCTATCGAAGCTAAACGTTTGAAAGACGCTGACCCGAAATTGATGACTGACGTATACGTATCTTTGTTCCAGGGTTTGTTCTTCTGGCCGCAGGTAATGCACTTGGACTTGGAACCGCACGGTGAAGAATTGTCTCAGAAAATCGAAACTTTGACTACTATCTTCTTGCAGTCTTACGGTGTAGCTGAATAA"
CODON={ 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L','ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P','ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E','TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
def tr(d): return "".join(CODON.get(d[i:i+3],'X') for i in range(0,len(d)-2,3))
WT_PROT=tr(WT_DNA).rstrip('*'); WL=len(WT_PROT)
SITES=[66,69,75,94,102,107,111,125,128,129,133,159,161,163,164,169,170]
wt_at={p:WT_PROT[p-1] for p in SITES}
ANCHOR="ATGCAGAAAA"
THR=10.0

DIRS=["./work/tsm_results","./work/tsm2"]
TARGETS=["CDCA","PCA","7K-LCA","23K-CDCA","LCA","LLDCA"]
def find_fq(t):
    for d in DIRS:
        for ext in [f"{t}_fastq",".fastq"]:
            p=os.path.join(d,f"{t}_fastq")
            if os.path.exists(p): return p
    return None

def iter_reads(p):
    op=gzip.open if p.endswith('.gz') else open
    fh=op(p,'rt'); h=fh.readline()
    while h:
        s=fh.readline().strip(); fh.readline(); fh.readline(); yield s; h=fh.readline()

def call_read(read):
    if ANCHOR in read:
        i=read.find(ANCHOR); g=read[i:i+612]
    elif Seq(read).reverse_complement().find(ANCHOR)>=0:
        r=str(Seq(read).reverse_complement()); i=r.find(ANCHOR); g=r[i:i+612]
    else: return None
    if len(g)!=612: return None
    p=tr(g).rstrip('*')
    if len(p)!=WL or '*' in p: return None
    return p

OUT="./work/results/per_target_six"; os.makedirs(OUT,exist_ok=True)
summary={}; tgt_site_aa={}; tgt_geno={}; tgt_mutsig={}; tgt_full={}
for t in TARGETS:
    p=find_fq(t)
    if not p: print(f"{t}: NO FASTQ"); continue
    n=0;ok=0; noa=0; bad=0
    site_aa=defaultdict(Counter); geno=Counter(); full=Counter()
    for s in iter_reads(p):
        n+=1; pr=call_read(s)
        if pr is None:
            if ANCHOR in s or Seq(s).reverse_complement().find(ANCHOR)>=0: bad+=1
            else: noa+=1
            continue
        ok+=1
        for i in range(WL): site_aa[i+1][pr[i]]+=1
        geno["".join(pr[p-1] for p in SITES)]+=1
        full[tuple(f"{WT_PROT[i]}{i+1}{pr[i]}" for i in range(WL) if pr[i]!=WT_PROT[i])]+=1
    summary[t]=dict(n=n,ok=ok,noanchor=noa,badframe=bad)
    tgt_site_aa[t]=site_aa; tgt_geno[t]=geno; tgt_full[t]=full
    # 位点频率
    rows=[]
    for pos in range(1,WL+1):
        c=site_aa[pos]; tot=sum(c.values()); wt_aa=WT_PROT[pos-1]
        nonwt=(tot-c.get(wt_aa,0))/tot*100 if tot else 0
        alts=[(a,cnt) for a,cnt in c.most_common(4) if a!=wt_aa and cnt/tot*100>=1]
        rows.append(dict(pos=pos,wt=wt_aa,depth=tot,nonwt_pct=nonwt,top_alt=alts[0] if alts else None))
    pd.DataFrame(rows).to_csv(f"{OUT}/{t}_position_freq.csv",index=False)
    sig=[r for r in rows if r["nonwt_pct"]>=THR]
    tgt_mutsig[t]=set(r["wt"]+str(r["pos"])+(r["top_alt"][0] if isinstance(r["top_alt"],tuple) else "") for r in sig)-{"R116Q"}
    print(f"\n=== {t} reads={n:,} 有效={ok:,} ===  显著位点(>=10%): {sorted(tgt_mutsig[t])}")
    for r in sorted(sig,key=lambda x:-x["nonwt_pct"])[:8]:
        insite=int(r["pos"]) in SITES
        print(f"  aa{int(r['pos']):>3}({r['wt']}{'*' if insite else ' '}) nonWT={r['nonwt_pct']:5.1f}% alt={r['top_alt']} [design:{insite}]")
    pd.DataFrame([dict(genotype=g,count=c,pct=c/ok*100) for g,c in geno.most_common(30)]).to_csv(f"{OUT}/{t}_17site_genotype_dist.csv",index=False)
pd.DataFrame(summary).T.to_csv(f"{OUT}/read_summary.csv")

# === 共享/特异 突变 (>=10%, 除骨架R116Q) ===
print("\n\n========== 6靶点 共享/特异 ==========")
allmuts=set().union(*tgt_mutsig.values())-{"R116Q"}
shared=set.intersection(*[tgt_mutsig[t] for t in TARGETS if tgt_mutsig[t]])-{"R116Q"}
print(f"各靶点显著突变(>=10%, 除R116Q):")
for t in TARGETS: print(f"  {t}: {sorted(tgt_mutsig.get(t,[]))}")
print(f"\n6靶点全共享: {sorted(shared) if shared else '∅'}")
for t in TARGETS:
    spec=tgt_mutsig.get(t,set())-set().union(*[tgt_mutsig.get(x,set()) for x in TARGETS if x!=t])
    print(f"{t}特异: {sorted(spec)}")
# 矩阵: 突变 x 靶点 频率
mrows=[]
for m in sorted(allmuts):
    row={"mutation":m}
    for t in TARGETS:
        # 找该突变在该靶点的频率
        pos=int(''.join(c for c in m[1:-1] if c.isdigit()))
        pf=pd.read_csv(f"{OUT}/{t}_position_freq.csv") if os.path.exists(f"{OUT}/{t}_position_freq.csv") else None
        if pf is not None and (pf.pos==pos).any():
            r=pf[pf.pos==pos].iloc[0]
            ta=r.top_alt
            alt_ch=ta[0] if isinstance(ta,tuple) else None
            if alt_ch==m[-1]: row[t]=round(r.nonwt_pct,1)
            else: row[t]=0.0
        else: row[t]=0.0
    mrows.append(row)
mm=pd.DataFrame(mrows).fillna(0)
mm.to_csv(f"{OUT}/mutation_by_target_matrix.csv",index=False)
print("\n突变×靶点 频率矩阵:"); print(mm.to_string(index=False))
print("\n产物:",os.listdir(OUT))
