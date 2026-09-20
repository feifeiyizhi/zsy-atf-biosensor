#!/usr/bin/env python3
# sequence_gate.py — 序列关(错误进管线即阻断,而非事后审计)
#
# 这道关专门拦上次乌龙的两类错误:
#   A. 解析错误 R116Q→位1(由 m[1:-3] 截错):解析器单测 R116Q→116/I73N→73/I142I→142
#   B. 参考系不一致:csv mutations 列 0810 是 A11 相对、7k-LCA/CDCA 是 WT 相对
#      → 全克隆 csv 突变列表与实际 aa 序列逐位比对,对不上即 exit 1(停止建模)
#
# 用法: python3 sequence_gate.py  # 通过=exit 0; 不通过=exit 1+详情
# 被建模脚本 import: from sequence_gate import parse_mutation, get_clone_sites, assert_gate
import csv, sys, os

BASE = "."
CSV = f"{BASE}/work/results/mutation_activity_full.csv"
WT_FA = f"{BASE}/work/ref/wt_reference_libbackbone.fasta"
WT_ALT_FA = f"{BASE}/work/ref/wt_reference_altcodon_userpasted.fasta"

CODON = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

def translate(dna):
    return ''.join(CODON.get(dna[i:i+3], '?') for i in range(0, len(dna)-2, 3))

def parse_mutation(m):
    """'R116Q' → ('R', 116, 'Q'); 'I142I'(silent) → ('I',142,'I'); '' → None.
    正确切片 m[1:-1](去首 WT aa + 尾 mut aa)= '116'。
    旧 bug m[1:-3] 截出 '1' → 这里单测必失败。"""
    m = m.strip()
    if not m or m == 'none':
        return None
    if len(m) < 3:
        raise ValueError(f"mutation too short: {m!r}")
    wt = m[0]; mut = m[-1]; pos_s = m[1:-1]
    if not pos_s.isdigit():
        raise ValueError(f"mutation pos not digit: {m!r} → pos={pos_s!r}")
    return wt, int(pos_s), mut

# ---- 单测:解析器正确性(拦 R116Q→1 类 bug) ----
def test_parser():
    cases = {'R116Q': ('R', 116, 'Q'), 'I73N': ('I', 73, 'N'), 'I142I': ('I', 142, 'I'),
             'Q136L': ('Q', 136, 'L'), 'D176N': ('D', 176, 'N'), 'L177S': ('L', 177, 'S')}
    for m, expect in cases.items():
        got = parse_mutation(m)
        assert got == expect, f"PARSE FAIL: {m} → {got} (expect {expect})"
        assert got[1] != 1 or m == 'X1X', f"PARSE BUG: {m} → pos=1 (the R116Q→1 regression!)"
    # 明确回归测试: m[1:-3] 旧 bug 会给 pos=1
    bad = 'R116Q'[1:-3]  # 旧 bug 写法
    assert bad == '1', f"regression fixture changed: {bad!r} (应为 '1', 证明旧 m[1:-3] 错)"
    print("[test_parser] OK — R116Q→116, I73N→73, I142I→142; m[1:-3] 回归测试就位")

def load_wt():
    fa = open(WT_FA).read().split('\n')
    dna = fa[1].strip()
    aa = translate(dna)
    # alt-codon 记录(旧密码子,蛋白等价但 reads 不匹配师兄库——单独记录,不混)
    alt = open(WT_ALT_FA).read().split('\n')[1].strip() if os.path.exists(WT_ALT_FA) else ''
    alt_aa = translate(alt) if alt else ''
    return dna, aa, alt, alt_aa

def get_clone_sites_aa(seq_aa, wt_aa):
    """从蛋白 seq 与 WT diff → 位点+突变(统一 WT 相对)。seq_aa 是蛋白。"""
    L = min(len(seq_aa), len(wt_aa))
    sites = []
    for i in range(L):
        a = seq_aa[i]; w = wt_aa[i]
        if a != w and a != '?' and w != '?':
            sites.append((w, i+1, a))  # (wt_aa, pos, mut_aa)
    return sites

def csv_muts_to_tuples(mut_str):
    """csv mutations 列 → set of (wt,pos,mut)。"""
    out = []
    for m in mut_str.split(';'):
        p = parse_mutation(m)
        if p: out.append(p)
    return set(out)

# ---- A11 origin(0810 的参考系) ----
A11_ORIGIN = {('R',116,'Q'), ('Q',136,'L'), ('Q',164,'R')}  # A11 = R116Q+Q136L+Q164R

def run_gate():
    test_parser()  # 拦解析 bug
    wt_dna, wt_aa, alt_dna, alt_aa = load_wt()
    assert len(wt_aa) == 204, f"WT len {len(wt_aa)}"
    rows = list(csv.DictReader(open(CSV)))

    print(f"\n[一致性检查] WT={len(wt_aa)}aa, 克隆={len(rows)}")
    print(f"  WT libbackbone: 旧密码子版(匹配师兄338lib reads)")
    print(f"  alt-codon: 蛋白{'等价' if alt_aa==wt_aa else '不等'}(reads 不匹配库,单独记录不用)")
    if alt_aa and alt_aa != wt_aa:
        print(f"  ⚠ alt-codon 蛋白与 libbackbone 不等!需查")

    blockers = []
    unified_ok = []
    for r in rows:
        seq_aa = (r.get('seq') or '').strip()
        csv_t = csv_muts_to_tuples(r['mutations'])
        seq_t = set(get_clone_sites_aa(seq_aa, wt_aa))
        # 直接比对(csv 应 == seq)
        if csv_t == seq_t:
            unified_ok.append((r['clone'], r['exp'], 'csv==seq', csv_t))
            continue
        # 不等:查是否 0810 的 A11 参考系(csv=A11相对,补 A11 后应==seq)
        if r['exp'] == '0810_23kPCA':
            csv_unified = csv_t | A11_ORIGIN
            if csv_unified == seq_t:
                unified_ok.append((r['clone'], r['exp'], '0810 A11-相对→补A11后一致', csv_t))
                continue
        # 仍不等:阻断
        blockers.append((r['clone'], r['exp'], csv_t, seq_t,
                          csv_t - seq_t, seq_t - csv_t))

    print(f"\n  一致: {len(unified_ok)}/{len(rows)}")
    for c,e,st,_ in unified_ok:
        print(f"    ✓ {c:<4}{e:<12} {st}")
    if blockers:
        print(f"\n  ✗ 阻断(对不上,停止建模): {len(blockers)}")
        for c,e,csvt,seqt,only_csv,only_seq in blockers:
            print(f"    {c:<4}{e:<12} csv有/seq无={only_csv}  seq有/csv无={only_seq}")
        print("\n  → 修复:统一参考系(0810 补 A11)或改 csv;不修不准建模。")
        return 1
    print("\n[GATE PASS] 全克隆 csv 突变 == seq aa 逐位一致(0810 补 A11 后),可建模。")
    return 0

if __name__ == '__main__':
    sys.exit(run_gate())
