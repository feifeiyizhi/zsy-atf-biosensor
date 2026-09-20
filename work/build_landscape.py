#!/usr/bin/env python3
"""
338lib1-5 负筛选(DNA结合正向选择) 数据处理 pipeline —— 纯Python, 忠实复刻 DATA_PREP 逻辑.
- 扩增子 = ref[175:530] = 355nt (锚点 AAAGCGTGCTGA @175). R1(250) + rc(R2)(250) 固定重叠145 -> window=R1[:250]+rc(R2)[145:250].
- 重建全长基因 = ref[:175]+window+ref[530:] (612nt) -> 翻译204aa -> 对WT找突变 -> 定位可变位点(17位点).
产物: results/ 下 逐lib DNA计数、蛋白变体×轮次矩阵、17位点突变表、位点谱、富集轨迹.
"""
import os, sys, gzip, subprocess, pickle
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd

WORK="./work"
FQDIR=os.path.join(WORK,"lib0_lib5negative_selection")
OUT=os.path.join(WORK,"results"); os.makedirs(OUT,exist_ok=True)
REF="ATGCAGAAAAAACTGACCCGCAGTCAGCAGAAACATCTGGATATTATTAACGCGGCGAAAGAAGAATTTATTGAATTTGGCTTTCTGGCGGCGAACATGGATCGCATTACGAGCAGCGCGGAAGTGAGCAAACGCACCCTGTATCGCCATTTTGAAAGCAAAGAAGTGCTGTTTGAAAGCGTGCTGACCATTATTAACGATAGCGTGAACGAAAGCATTAGCTATCATTTTGATCCGAACAAAAGCACCGAAGAACAGCTGACCGAAATTGCGTATAAAGAAATTGATGTGCTGTATAAAACCTATGGCATTGCGCTGGCGCGCACCATTGTGATGGAATTTCTGCGTCAGCCGGAAATGGCGAAAACCCTGATTCAGAACATTTATAGCATTCGCGCGATTACGCAGTGGTTTCGCAGCGCGATTGAAGCGAAACGCCTGAAAGATGCGGATCCGAAACTGATGACCGATGTGTATGTGAGCCTGTTTCAAGGCCTGTTCTTTTGGCCGCAAGTGATGCATCTGGATCTGGAACCGCATGGCGAAGAACTGAGTCAGAAAATTGAAACCCTGACCACCATTTTTCTGCAGAGCTATGGCGTGGCGGAATAA"
ANCHOR="AAAGCGTGCTGA"; ASTART=REF.find(ANCHOR); WLEN=355; OV=145; RLEN=250
WEND=ASTART+WLEN  # 530
PREFIX=REF[:ASTART]; SUFFIX=REF[WEND:]
LIBS=[1,2,3,4,5]
_comp=str.maketrans("ACGTN","TGCAN")
def rc(s): return s.translate(_comp)[::-1]
CODON={ 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L','ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P','ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E','TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
def translate(dna):
    return "".join(CODON.get(dna[i:i+3],'X') for i in range(0,len(dna)-2,3))
WT_PROT=translate(REF)  # 含末尾*

def process_lib(lib):
    r1=os.path.join(FQDIR,f"338lib{lib}",f"338lib{lib}_1.fq.gz")
    r2=os.path.join(FQDIR,f"338lib{lib}",f"338lib{lib}_2.fq.gz")
    p1=subprocess.Popen(["pigz","-dc",r1],stdout=subprocess.PIPE,text=True,bufsize=1<<20)
    p2=subprocess.Popen(["pigz","-dc",r2],stdout=subprocess.PIPE,text=True,bufsize=1<<20)
    cnt=Counter(); total=0; kept=0; badlen=0; noanchor=0
    f1,f2=p1.stdout,p2.stdout
    while True:
        h1=f1.readline()
        if not h1: break
        s1=f1.readline().strip(); f1.readline(); f1.readline()
        f2.readline(); s2=f2.readline().strip(); f2.readline(); f2.readline()
        total+=1
        if len(s1)!=RLEN or len(s2)!=RLEN: badlen+=1; continue
        window=s1[:RLEN]+rc(s2)[OV:RLEN]     # 355nt
        if len(window)!=WLEN: badlen+=1; continue
        if not window.startswith(ANCHOR): noanchor+=1; continue
        cnt[window]+=1; kept+=1
    p1.wait(); p2.wait()
    return lib,cnt,dict(total=total,kept=kept,badlen=badlen,noanchor=noanchor)

def main():
    print("WT protein(204aa):",WT_PROT.rstrip('*'))
    print(f"anchor@{ASTART}, window {ASTART}:{WEND} (len{WLEN}), prefix{len(PREFIX)} suffix{len(SUFFIX)}")
    stats={}; libcnt={}
    with ProcessPoolExecutor(max_workers=5) as ex:
        for lib,cnt,st in ex.map(process_lib,LIBS):
            libcnt[lib]=cnt; stats[lib]=st
            print(f"lib{lib}: total={st['total']:,} kept={st['kept']:,} "
                  f"badlen={st['badlen']:,} noanchor={st['noanchor']:,} unique={len(cnt):,}")
    pd.DataFrame(stats).T.to_csv(os.path.join(OUT,"read_stats.csv"))

    # ---- 蛋白级聚合: window(DNA) -> 全长基因 -> 蛋白 -> 突变集 ----
    # 逐 lib: 过滤 count>=3, 翻译, 聚到蛋白变体
    prot_lib=defaultdict(lambda: defaultdict(int))  # prot -> {lib:count}
    prot_muts={}
    for lib in LIBS:
        for win,c in libcnt[lib].items():
            if c<3: continue
            full=PREFIX+win+SUFFIX
            prot=translate(full).rstrip('*')
            prot_lib[prot][lib]+=c
    # 突变集 vs WT
    wt=WT_PROT.rstrip('*')
    for prot in prot_lib:
        if len(prot)!=len(wt):
            prot_muts[prot]=("LEN%d"%len(prot),)  # 移码/长度异常
            continue
        prot_muts[prot]=tuple(f"{wt[i]}{i+1}{prot[i]}" for i in range(len(wt)) if prot[i]!=wt[i])

    # 构建矩阵
    rows=[]
    for prot,libd in prot_lib.items():
        muts=prot_muts[prot]
        row={"protein":prot,"n_mut":(0 if muts==() else (999 if muts and str(muts[0]).startswith("LEN") else len(muts))),
             "mutations":";".join(muts)}
        for lib in LIBS: row[f"lib{lib}_count"]=libd.get(lib,0)
        rows.append(row)
    df=pd.DataFrame(rows)
    for lib in LIBS:
        tot=df[f"lib{lib}_count"].sum()
        df[f"lib{lib}_pct"]=df[f"lib{lib}_count"]/tot*100 if tot else 0
    df["sum_count"]=df[[f"lib{l}_count" for l in LIBS]].sum(axis=1)
    # 富集轨迹: lib5 vs lib1 (lib0随机池未测序, 暂用lib1做最早基线)
    eps=np.finfo(float).eps
    df["enrich_lib5_over_lib1"]=np.log2((df["lib5_pct"]+eps)/(df["lib1_pct"]+eps))
    df=df.sort_values(["sum_count"],ascending=False)
    df.to_csv(os.path.join(OUT,"protein_variant_by_round.csv"),index=False)

    # ---- 定位可变位点(=设计的17位点) ----
    wtlen=len(wt)
    pos_var=Counter()  # 位点 -> 携带该位点突变的(过滤后)变体数(按sum_count加权)
    aa_at=defaultdict(Counter)
    for prot,libd in prot_lib.items():
        if len(prot)!=wtlen: continue
        w=sum(libd.values())
        for i in range(wtlen):
            if prot[i]!=wt[i]:
                pos_var[i+1]+=1
                aa_at[i+1][prot[i]]+=w
    posdf=pd.DataFrame([{"aa_pos":p,"wt_aa":wt[p-1],"n_variants_mutated":pos_var[p],
                         "alt_aa_spectrum":";".join(f"{a}:{c}" for a,c in aa_at[p].most_common())}
                        for p in sorted(pos_var, key=lambda x:-pos_var[x])])
    posdf.to_csv(os.path.join(OUT,"mutated_position_spectrum.csv"),index=False)
    top_sites=sorted(pos_var.items(), key=lambda x:-x[1])[:25]
    print("\n== 变异最频繁的位点(候选17位点) ==")
    for p,n in top_sites:
        print(f"  aa{p} ({wt[p-1]}): {n} variants, alt={dict(aa_at[p].most_common(4))}")

    print(f"\n过滤后蛋白变体数(count>=3, 任一lib): {len(prot_lib):,}")
    print(f"其中长度正常(=204aa剪去*后{wtlen}): {sum(1 for p in prot_lib if len(p)==wtlen):,}")
    print("产物写入:",OUT)
    print(df.head(10).to_string())

if __name__=="__main__":
    main()
