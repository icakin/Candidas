#!/usr/bin/env bash
# Run Seq2Tm on ribosomal proteins per species.
# Run from the repo root: bash run_ribo_seq2tm.sh
# Output: ribo_tm_predictions.csv (id=species|gene, pred_tm)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GEM="$HERE/gem"
export TORCH_HOME="${TORCH_HOME:-$GEM/external/torch_home}"
CKPT="$GEM/external/large_model_pth/model_tm_window=3_r2=0.76.pth"
CODE="$GEM/external/Seq2Topt/code"
OUT="$HERE/ribo_tm_predictions.csv"

echo "Extracting ribosomal sequences..."
python3 - <<'EOF'
import os, re, pandas as pd
os.chdir(os.environ.get('HERE', '.'))

def fasta(p):
    h,s=None,[]
    for l in open(p):
        if l.startswith('>'):
            if h: yield h,''.join(s)
            h,s=l[1:].strip(),[]
        else: s.append(l.strip())
    if h: yield h,''.join(s)

Oh = dict(pd.read_csv('gem/tables/rbh/pid_haemulonii.tsv',sep='\t',header=None,names=['a','b','p'])[['a','b']].values)
Od = dict(pd.read_csv('gem/tables/rbh/pid_duobushaemulonii.tsv',sep='\t',header=None,names=['a','b','p'])[['a','b']].values)
HA = {h.split()[0]:s for h,s in fasta('phylo/proteomes/haemulonii.faa')}
DU = {h.split()[0]:s for h,s in fasta('phylo/proteomes/duobushaemulonii.faa')}
AU_ribo = {h.split('|')[1]:seq for h,seq in fasta('gem/tables/machinery_ribosome.faa') if h.split('|')[0]=='auris'}
shared = [g for g in AU_ribo if g in Oh and Oh[g] in HA and g in Od and Od[g] in DU]

with open('/tmp/ribo_seq2tm.faa','w') as f:
    for g in shared:
        f.write(f'>auris|{g}\n{AU_ribo[g]}\n')
        f.write(f'>haemulonii|{Oh[g]}\n{HA[Oh[g]]}\n')
        f.write(f'>duobushaemulonii|{Od[g]}\n{DU[Od[g]]}\n')
    for h,seq in fasta('phylo/proteomes/parapsilosis.faa'):
        if re.search(r'\b(40S|60S) ribosomal protein\b', h, re.I):
            f.write(f'>parapsilosis|{h.split()[0]}\n{seq}\n')

n = sum(1 for l in open('/tmp/ribo_seq2tm.faa') if l.startswith('>'))
print(f"  {n} sequences written (124 orthologs × 3 species + parapsilosis 40S/60S)")
EOF

echo "Running Seq2Tm..."
python3 "$GEM/15_run_seq2tm.py" /tmp/ribo_seq2tm.faa "$OUT" --ckpt "$CKPT" --code "$CODE"
echo "Done: $OUT"
