#!/usr/bin/env python3
"""Extract mutation-sensitive structural features from paired PDB structures.

The script is deliberately conservative: it reports observed coordinate changes
only when both WT and mutant PDBs are supplied. It never treats a missing model
as a zero structural delta.
"""
from __future__ import annotations
import argparse,csv,json,math,re
from pathlib import Path

MUT_RE = re.compile(r'^([A-Z])([0-9]+)([A-Z])$')

def parse_ca(path):
    atoms = {}
    for line in path.read_text(errors='replace').splitlines():
        if not line.startswith(('ATOM  ','HETATM')) or line[12:16].strip() != 'CA': continue
        try:
            chain=line[21].strip() or '_'; res=int(line[22:26]); xyz=tuple(float(line[i:i+8]) for i in (30,38,46)); b=float(line[60:66])
        except (ValueError,IndexError): continue
        atoms[(chain,res)]={'xyz':xyz,'b_factor':b,'aa':line[17:20].strip()}
    return atoms

def dist(a,b): return math.sqrt(sum((x-y)**2 for x,y in zip(a,b)))

def main():
    p=argparse.ArgumentParser(); p.add_argument('--wt',type=Path,required=True); p.add_argument('--mutant',type=Path,required=True); p.add_argument('--mutation',required=True); p.add_argument('--out',type=Path,required=True); p.add_argument('--chain',default='A'); p.add_argument('--contact-cutoff',type=float,default=8.0); a=p.parse_args()
    m=MUT_RE.match(a.mutation);
    if not m: raise SystemExit('mutation must be like Q164R')
    pos=int(m.group(2)); wt=parse_ca(a.wt); mut=parse_ca(a.mutant); key=(a.chain,pos)
    if key not in wt or key not in mut: raise SystemExit(f'missing residue in both structures: {key}')
    changed=[]; common=set(wt)&set(mut)
    for k in sorted(common):
        delta=dist(wt[k]['xyz'],mut[k]['xyz'])
        if delta>0: changed.append((k,delta))
    wt_contacts=[k for k,v in wt.items() if k!=key and dist(v['xyz'],wt[key]['xyz'])<=a.contact_cutoff]
    mut_contacts=[k for k,v in mut.items() if k!=key and dist(v['xyz'],mut[key]['xyz'])<=a.contact_cutoff]
    row={'mutation':a.mutation,'wt_structure':str(a.wt),'mutant_structure':str(a.mutant),'status':'OBSERVED_STRUCTURE_DELTA','mutant_ca_displacement':dist(wt[key]['xyz'],mut[key]['xyz']),'n_common_ca':len(common),'n_residues_ca_delta_gt_0.5A':sum(d>0.5 for _,d in changed),'max_ca_delta':max((d for _,d in changed),default=0.0),'mean_ca_delta':sum(d for _,d in changed)/len(changed) if changed else 0.0,'wt_contact_count':len(wt_contacts),'mutant_contact_count':len(mut_contacts),'contact_gain':len(set(mut_contacts)-set(wt_contacts)),'contact_loss':len(set(wt_contacts)-set(mut_contacts)),'wt_b_factor_mutant_site':wt[key]['b_factor'],'mutant_b_factor_mutant_site':mut[key]['b_factor']}
    a.out.parent.mkdir(parents=True,exist_ok=True); a.out.write_text(json.dumps(row,indent=2)+'\n'); print(json.dumps(row,indent=2))
if __name__=='__main__':main()
