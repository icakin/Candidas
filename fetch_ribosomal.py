import os, sys, time, json, urllib.parse, urllib.request, csv
import pandas as pd
OUT='phylo/ribo'; os.makedirs(OUT,exist_ok=True)
M=pd.read_csv('phylo/inputs/thermal_manifest.csv')
M=M[M.proteome_source.astype(str).str.lower().eq('uniprot')]
def binom(s):
    w=[x for x in str(s).replace('(','').replace(')','').split() if x[0].isalpha()]
    return ' '.join(w[:2])
def get(url,tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url,timeout=90) as r: return r.read().decode()
        except Exception as e:
            if i==tries-1: return None
            time.sleep(3*(i+1))
log=open('phylo/ribo/_log.csv','a',newline='')
W=csv.writer(log)
done={f[:-4] for f in os.listdir(OUT) if f.endswith('.faa')}
for n,r in enumerate(M.itertuples()):
    sp=binom(r.species); slug=sp.replace(' ','_')
    if slug in done: continue
    q=urllib.parse.quote(f'organism_name:"{sp}"')
    j=get(f'https://rest.uniprot.org/proteomes/search?query={q}&format=json&size=5')
    up=None
    if j:
        try:
            for e in json.loads(j).get('results',[]):
                if e.get('proteomeType','').startswith(('Reference','Other')) or True:
                    up=e['id']; break
        except Exception: pass
    if not up:
        W.writerow([sp,'NO_PROTEOME','','']); log.flush(); continue
    fq=urllib.parse.quote(f'(proteome:{up}) AND (protein_name:"ribosomal protein")')
    fa=get(f'https://rest.uniprot.org/uniprotkb/stream?query={fq}&format=fasta')
    nseq=fa.count('>') if fa else 0
    if nseq>=20:
        open(f'{OUT}/{slug}.faa','w').write(fa)
        W.writerow([sp,up,nseq,r.tmax])
    else:
        W.writerow([sp,up,f'ONLY_{nseq}',r.tmax])
    log.flush()
    if n%25==0: print(f'{n}/{len(M)}  {sp}  {nseq} seqs',flush=True)
    time.sleep(0.3)
print('done')
