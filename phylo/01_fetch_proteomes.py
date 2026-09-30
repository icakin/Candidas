#!/usr/bin/env python3
"""
fetch_proteomes.py — download annotated proteomes via the NCBI Datasets v2 API.

The old Assembly database (db=assembly) was retired in 2024, which is why the
esearch-based script returned nothing. This uses the supported REST API.

Run on YOUR Mac (needs internet):
    python3 ~/Desktop/Projects/Candidas/phylo/fetch_proteomes.py

Writes protein FASTAs to phylo/proteomes/ and prints exactly what it picked.
Uses only the standard library.
"""
import io, json, os, re, sys, time, urllib.parse, urllib.request, zipfile

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "proteomes")
os.makedirs(OUT, exist_ok=True)

# (tag, NCBI taxon name, strain hint or None)
TARGETS = [
    ("auris_cladeI",     "Candidozyma auris",             "B8441"),
    ("auris_cladeII",    "Candidozyma auris",             "B11220"),
    ("auris_cladeIII",   "Candidozyma auris",             "B11221"),
    ("auris_cladeIV",    "Candidozyma auris",             "B11245"),
    ("haemulonii",       "Candidozyma haemuli",           None),
    ("duobushaemulonii", "Candidozyma duobushaemuli",     None),
    ("pseudohaemulonii", "Candidozyma pseudohaemuli",     None),
    ("parapsilosis",     "Candida parapsilosis",          "CDC317"),
    ("albicans",         "Candida albicans",              "SC5314"),
    ("lusitaniae",       "Clavispora lusitaniae",         None),
]

def get(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "candidas-phylo/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(3)

def reports(taxon):
    """All annotated assemblies for a taxon."""
    url = (f"{API}/genome/taxon/{urllib.parse.quote(taxon)}/dataset_report"
           f"?filters.has_annotation=true&page_size=200")
    try:
        data = json.loads(get(url))
    except Exception as e:
        print(f"     ! API error for {taxon}: {e}")
        return []
    return data.get("reports", [])

def strain_of(rep):
    inf = (rep.get("assembly_info", {}) or {}).get("biosample", {}) or {}
    bits = []
    for k in ("strain", "isolate"):
        v = ((rep.get("assembly_info", {}) or {}).get("infraspecific_names", {}) or {}).get(k)
        if v: bits.append(v)
    for a in inf.get("attributes", []) or []:
        if a.get("name") in ("strain", "isolate") and a.get("value"):
            bits.append(a["value"])
    name = (rep.get("assembly_info", {}) or {}).get("assembly_name", "")
    bits.append(name)
    return " ".join(str(b) for b in bits)

def n_proteins(rep):
    st = (rep.get("annotation_info", {}) or {}).get("stats", {}) or {}
    gc = st.get("gene_counts", {}) or {}
    return gc.get("protein_coding") or 0

def pick(reps, hint):
    if not reps: return None, []
    scored = []
    for r in reps:
        acc = r.get("accession", "")
        s = 0
        if hint and hint.lower() in strain_of(r).lower(): s += 100
        if acc.startswith("GCF_"): s += 20            # RefSeq annotation preferred
        cat = (r.get("assembly_info", {}) or {}).get("refseq_category", "")
        if cat: s += 10
        if n_proteins(r): s += 5
        scored.append((s, acc, r))
    scored.sort(key=lambda x: (-x[0], x[1]))
    if hint and scored[0][0] < 100:
        return None, scored[:8]                        # hint not matched -> report options
    return scored[0][2], scored[:8]

def download_prot(acc, tag):
    url = (f"{API}/genome/accession/{acc}/download"
           f"?include_annotation_type=PROT_FASTA&filename={tag}.zip")
    blob = get(url)
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        names = [n for n in z.namelist() if n.endswith("protein.faa")]
        if not names:
            return None
        with z.open(names[0]) as fh:
            data = fh.read()
    path = os.path.join(OUT, f"{tag}.faa")
    with open(path, "wb") as f:
        f.write(data)
    return path

print(f"Downloading proteomes -> {OUT}\n")
ok = fail = 0
for tag, taxon, hint in TARGETS:
    print(f"[{tag}] {taxon}" + (f"  strain~{hint}" if hint else ""))
    reps = reports(taxon)
    if not reps:
        print("     FAIL: no annotated assemblies returned\n"); fail += 1; continue
    rep, options = pick(reps, hint)
    if rep is None:
        print(f"     FAIL: strain '{hint}' not matched. Candidates:")
        for s, acc, r in options:
            print(f"        {acc:20} {n_proteins(r):>6} prot  {strain_of(r)[:60]}")
        print(); fail += 1; continue
    acc = rep["accession"]
    try:
        path = download_prot(acc, tag)
    except Exception as e:
        print(f"     FAIL download {acc}: {e}\n"); fail += 1; continue
    if not path:
        print(f"     FAIL: {acc} has no protein.faa\n"); fail += 1; continue
    n = sum(1 for l in open(path) if l.startswith(">"))
    print(f"     OK  {acc}  {n} proteins  [{strain_of(rep)[:50]}]\n")
    ok += 1
    time.sleep(0.5)

print("=" * 60)
print(f"downloaded {ok}, failed {fail}")
for f in sorted(os.listdir(OUT)):
    p = os.path.join(OUT, f)
    print(f"  {f:26} {os.path.getsize(p)/1e6:7.1f} MB")
