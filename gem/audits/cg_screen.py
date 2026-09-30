#!/usr/bin/env python3
import re, subprocess, os, collections
from Bio import SeqIO
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
WORK = TABLES / 'cg_work'; WORK.mkdir(exist_ok=True)   # DIAMOND db and hits, gitignored
SPP=[('auris_cladeI','C.auris_I'),('auris_cladeII','C.auris_II'),('auris_cladeIII','C.auris_III'),
     ('auris_cladeIV','C.auris_IV'),('haemulonii','C.haemulonii'),('duobushaemulonii','C.duobus'),
     ('pseudohaemulonii','C.pseudohae'),('parapsilosis','C.parapsilosis'),('lusitaniae','C.lusitaniae'),
     ('albicans','C.albicans')]
# pathway keyword -> gene family label (searched in albicans product descriptions)
GROUPS=[
 ('HSF1',['heat shock transcription factor','heat-shock transcription factor']),
 ('HSP90',['heat shock protein 90','hsp90']),
 ('HSP104_ClpB',['heat shock protein 104','hsp104','clpb','disaggregatase']),
 ('HSP70_SSA',['heat shock protein ssa','heat shock protein 70','hsp70','hsp71','ssa family']),
 ('HSP78',['heat shock protein 78']),
 ('HSP12',['heat shock protein 12']),
 ('Calcineurin_CNA_CNB',['calcineurin','phosphatase 2b','phosphatase 3']),
 ('CRZ1',['calcineurin-responsive','crz1']),
 ('AdenylateCyclase_CYR1',['adenylate cyclase']),
 ('PKA_TPK',['camp-dependent protein kinase','cyclic amp-dependent']),
 ('TPS1',['trehalose-6-phosphate synthase','alpha,alpha-trehalose-phosphate synthase']),
 ('TPS2',['trehalose-6-phosphate phosphatase','trehalose-phosphatase']),
 ('Trehalase_NTH',['trehalase']),
 ('ERG_sterol',['ergosterol','lanosterol','squalene','c-14 sterol','c-4 methylsterol','c-8 sterol','sterol 24-c','delta(24)','sterol reductase','sterol c-','methylsterol','erg']),
 ('OLE1_desaturase',['fatty acid desaturase','acyl-coa desaturase','delta-9']),
 ('FKS_glucansynthase',['1,3-beta-glucan synthase','1,3-beta-d-glucan synthase','glucan synthase']),
 ('CHS_chitinsynthase',['chitin synthase']),
 ('PKC1',['protein kinase c']),
 ('SOD',['superoxide dismutase']),
 ('Catalase',['catalase']),
 ('Glutaredoxin_GRX',['glutaredoxin']),
 ('FerricReductase',['ferric reductase','ferric-chelate reductase','ferric-chelate']),
 ('Siderophore_transp',['siderophore','ferrichrome']),
]
# 1) build albicans candidate queries
alb=list(SeqIO.parse(PROTEOMES / 'albicans.faa','fasta'))
def desc(r): return r.description.lower()
qrecs=[]; qgroup={}
used=set()
for gname,kws in GROUPS:
    for r in alb:
        d=desc(r)
        if any(k in d for k in kws) and r.id not in used:
            qrecs.append(r); qgroup[r.id]=gname; used.add(r.id)
with open(TABLES / 'cg_work' / 'cg_queries.faa','w') as fh:
    for r in qrecs: fh.write(f'>{r.id}\n{str(r.seq)}\n')
gc=collections.Counter(qgroup.values())
print('queries per group (from albicans):')
for g,_ in GROUPS: print(f'  {g:24s} {gc.get(g,0)}')
print('total queries:',len(qrecs))
# 2) build combined species-tagged db
with open(TABLES / 'cg_work' / 'cg_all.faa','w') as fh:
    for fpre,lab in SPP:
        for r in SeqIO.parse(PROTEOMES / f'{fpre}.faa','fasta'):
            fh.write(f'>{lab}|{r.id}\n{str(r.seq)}\n')
subprocess.run('diamond makedb --in cg_all.faa -d cg_all -p 4 --quiet',shell=True,check=True,cwd=WORK)
# 3) diamond blastp queries vs all
subprocess.run('diamond blastp -q cg_queries.faa -d cg_all -o cg_hits.tsv -p 4 --quiet '
               '--outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore '
               '--max-target-seqs 2000 --evalue 1e-20',shell=True,check=True,cwd=WORK)
print('blast done')
