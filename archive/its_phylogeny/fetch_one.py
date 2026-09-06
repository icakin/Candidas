#!/usr/bin/env python3
"""Download one proteome by exact accession.
   python3 fetch_one.py GCF_000182765.1 parapsilosis
"""
import io, os, sys, urllib.request, zipfile

acc, tag = sys.argv[1], sys.argv[2]
OUT = os.path.expanduser("~/Desktop/Projects/Candidas/phylo/proteomes")
os.makedirs(OUT, exist_ok=True)
url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{acc}"
       f"/download?include_annotation_type=PROT_FASTA&filename={tag}.zip")
req = urllib.request.Request(url, headers={"User-Agent": "candidas-phylo/1.0"})
blob = urllib.request.urlopen(req, timeout=180).read()
with zipfile.ZipFile(io.BytesIO(blob)) as z:
    names = [n for n in z.namelist() if n.endswith("protein.faa")]
    if not names:
        sys.exit(f"FAIL: no protein.faa in {acc}")
    data = z.open(names[0]).read()
p = os.path.join(OUT, f"{tag}.faa")
open(p, "wb").write(data)
n = sum(1 for l in open(p) if l.startswith(">"))
print(f"OK  {acc} -> {tag}.faa  ({n} proteins, {os.path.getsize(p)/1e6:.1f} MB)")
