import re
from pathlib import Path
SPS = {'auris':'auris_cladeI','haemulonii':'haemulonii',
       'duobushaemulonii':'duobushaemulonii','parapsilosis':'parapsilosis'}
SETS = {
 'ribosome':  re.compile(r'ribosomal protein|ribosomal.*subunit', re.I),
 'chaperone': re.compile(r'chaperon|heat shock|HSP\d|prefoldin|co-chaperone', re.I),
 'translate': re.compile(r'elongation factor|initiation factor|release factor|'
                         r'aminoacyl|tRNA synthetase|tRNA ligase', re.I),
}
EXCL = re.compile(r'ribosomal RNA|rRNA (methyl|pseudo)|biogenesis|assembly|processing', re.I)
P = Path('phylo/proteomes'); OUT = Path('gem/tables'); OUT.mkdir(exist_ok=True)
def fasta(p):
    h, s = None, []
    for line in open(p):
        if line.startswith('>'):
            if h: yield h, ''.join(s)
            h, s = line[1:].strip(), []
        else: s.append(line.strip())
    if h: yield h, ''.join(s)
for name, pat in SETS.items():
    rows, counts = [], {}
    for sp, f in SPS.items():
        n = 0
        for h, seq in fasta(P / f'{f}.faa'):
            if pat.search(h) and not EXCL.search(h) and 30 <= len(seq) <= 3000:
                rows.append((f"{sp}|{h.split()[0]}", seq)); n += 1
        counts[sp] = n
    out = OUT / f'machinery_{name}.faa'
    with open(out, 'w') as fh:
        for i, s in rows: fh.write(f">{i}\n{s}\n")
    print(f"{name:10s} {counts}  -> {out}")
