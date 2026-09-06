#!/usr/bin/env python3
import subprocess, collections
from Bio import SeqIO
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
WORK = TABLES / 'cg_work'; WORK.mkdir(exist_ok=True)   # DIAMOND db and hits, gitignored
SPP=[('auris_cladeI','C.auris_I'),('auris_cladeII','C.auris_II'),('auris_cladeIII','C.auris_III'),
     ('auris_cladeIV','C.auris_IV'),('haemulonii','C.haemulonii'),('duobushaemulonii','C.duobus'),
     ('pseudohaemulonii','C.pseudohae'),('parapsilosis','C.parapsilosis'),('lusitaniae','C.lusitaniae'),
     ('albicans','C.albicans')]
GROUPS=[
 ('HSP90',['heat shock protein 90','hsp90 family chaperone','hsp90']),
 ('HSP70_SSA',['heat shock protein ssa','hsp70 family chaperone','hsp70 family atpase']),
 ('HSP12',['hsp12p']),
 ('Calcineurin_CNA_CNB',['calcineurin','phosphatase 2b','phosphatase 3']),
 ('AdenylateCyclase_CYR1',['adenylate cyclase']),
 ('PKA_TPK',['camp-dependent protein kinase','cyclic amp-dependent']),
 ('TPS1',['trehalose-6-phosphate synthase']),
 ('TPS2',['trehalose-6-phosphate phosphatase','trehalose-phosphatase']),
 ('Trehalase_NTH',['trehalase']),
 ('ERG_sterol',['ergosterol','lanosterol','squalene','c-14 sterol','c-4 methylsterol','c-8 sterol','c-5 sterol','c-22 sterol','sterol 24-c','sterol reductase','methylsterol']),
 ('OLE1_desaturase',['stearoyl-coa','fatty acid desaturase','acyl-coa desaturase']),
 ('FKS_glucansynthase',['1,3-beta-glucan synthase','glucan synthase']),
 ('CHS_chitinsynthase',['chitin synthase']),
 ('PKC1',['protein kinase c']),
 ('SOD',['superoxide dismutase']),
 ('Catalase',['catalase']),
 ('Glutaredoxin_GRX',['glutaredoxin']),
 ('FerricReductase',['ferric reductase','ferric-chelate']),
 ('Siderophore_transp',['siderophore','ferrichrome']),
]
alb=list(SeqIO.parse(PROTEOMES / 'albicans.faa','fasta'))
qrecs=[]; qgroup={}; used=set()
for gname,kws in GROUPS:
    for r in alb:
        d=r.description.lower()
        if any(k in d for k in kws) and r.id not in used:
            qrecs.append(r); qgroup[r.id]=gname; used.add(r.id)
open(TABLES / 'cg_work' / 'cg_queries.faa','w').write(''.join(f'>{r.id}\n{r.seq}\n' for r in qrecs))
subprocess.run('diamond blastp -q cg_queries.faa -d cg_all -o cg_hits.tsv -p 4 --quiet '
 '--outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore --max-target-seqs 3000 --evalue 1e-20',
 shell=True,check=True,cwd=WORK)
# best group per subject
best={}
for ln in open(TABLES / 'cg_work' / 'cg_hits.tsv'):
    q,s,pid,length,qlen,slen,ev,bit=ln.split('\t')
    pid=float(pid); length=int(length); qlen=int(qlen); slen=int(slen); bit=float(bit)
    if pid<30: continue
    if length/qlen<0.5 or length/slen<0.5: continue
    if s not in best or bit>best[s][1]: best[s]=(qgroup[q],bit)
cnt=collections.defaultdict(lambda: collections.Counter())
for s,(g,b) in best.items():
    sp=s.split('|')[0]; cnt[g][sp]+=1
labs=[l for _,l in SPP]
glabs=[g for g,_ in GROUPS]
# print matrix
w=22
hdr=['auris_I','auris_II','auris_III','auris_IV','hae','duo','pseudo','para','lus','alb']
print('gene_family'.ljust(w)+''.join(h.rjust(9) for h in hdr))
for g in glabs:
    row=g.ljust(w)
    for l in labs:
        row+=str(cnt[g].get(l,0)).rjust(9)
    print(row)
