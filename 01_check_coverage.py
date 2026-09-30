#!/usr/bin/env python3
"""01_check_coverage.py -- how many yeasts actually have a growth temperature?

This is the gating question for the yeast-wide approach. If only a few dozen yeast
species have both a recorded growth temperature and a sequenced genome, the approach is
no better than the four species we already have and we stop. Run it before anything else.

Source: Engqvist (2018) BMC Microbiology 18:177, growth temperatures for 21,498
microorganisms, deposited at Zenodo record 1175608. Roughly 2% are eukaryotes.

    python3 01_check_coverage.py

Writes growth_temperatures.tsv and prints the counts. Needs internet; run it in your own
Terminal (the sandboxes I use are blocked from Zenodo).
"""
import json, os, sys, urllib.request, csv, collections

REC = "https://zenodo.org/api/records/1175608"

def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "candidas-yeastwide/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read()

def main():
    print("fetching the Zenodo record ...")
    rec = json.loads(get(REC))
    files = rec.get("files", [])
    print(f"record: {rec.get('metadata',{}).get('title','?')}")
    for f in files:
        print(f"   file: {f.get('key')}  {f.get('size')} bytes")
    # prefer a tab-delimited file over the xml
    cand = [f for f in files if f.get("key","").lower().endswith((".tsv", ".txt", ".csv"))]
    if not cand:
        sys.exit("no tabular file found; open the record by hand and tell me the file name")
    f = max(cand, key=lambda x: x.get("size", 0))
    url = f.get("links", {}).get("self") or f.get("links", {}).get("download")
    print(f"\ndownloading {f['key']} ...")
    raw = get(url)
    open(f["key"], "wb").write(raw)
    print(f"saved {f['key']} ({len(raw)} bytes)")

    # parse: sniff the delimiter, find the columns
    txt = raw.decode("utf-8", "replace").splitlines()
    delim = "\t" if txt[0].count("\t") >= txt[0].count(",") else ","
    rows = list(csv.DictReader(txt, delimiter=delim))
    cols = rows[0].keys() if rows else []
    print(f"\n{len(rows)} rows, columns: {list(cols)}")

    def col(*names):
        for n in names:
            for c in cols:
                if n.lower() in c.lower(): return c
        return None
    c_org  = col("organism", "species", "name")
    c_temp = col("temperature", "ogt", "topt")
    c_dom  = col("superkingdom", "domain", "kingdom")
    c_lin  = col("lineage", "taxonomy", "phylum", "class")
    print(f"using: organism={c_org!r} temperature={c_temp!r} domain={c_dom!r} lineage={c_lin!r}")

    # how many are fungi / yeasts?
    YEAST = ("saccharomyc", "candida", "candidozyma", "pichia", "kluyveromyc", "debaryomyc",
             "clavispora", "yarrowia", "ogataea", "komagataella", "lodderomyces", "meyerozyma",
             "torulaspora", "zygosaccharomyces", "lachancea", "nakaseomyces", "cryptococc",
             "rhodotorula", "malassezia", "schizosaccharomyc", "hanseniaspora", "wickerhamomyc",
             "starmerella", "brettanomyces", "dekkera", "metschnikowia", "cyberlindnera")
    n_euk = n_fungi = n_yeast = 0
    genera = collections.Counter(); cand_rows = []
    for r in rows:
        blob = " ".join(str(r.get(c, "")) for c in (c_org, c_dom, c_lin) if c).lower()
        if "eukary" in blob: n_euk += 1
        if "fungi" in blob or "ascomycot" in blob or "basidiomycot" in blob: n_fungi += 1
        if any(k in blob for k in YEAST):
            n_yeast += 1; cand_rows.append(r)
            genera[str(r.get(c_org, "")).split()[0] if r.get(c_org) else "?"] += 1
    print(f"\neukaryotes: {n_euk}   fungi: {n_fungi}   yeast-like by genus: {n_yeast}")
    print("top yeast genera present:")
    for g, n in genera.most_common(20): print(f"   {g:<24}{n}")

    # do OUR four species appear?
    print("\nour species in the dataset:")
    for want in ("auris", "haemulonii", "duobushaemulonii", "parapsilosis", "albicans"):
        hits = [r for r in rows if want in str(r.get(c_org, "")).lower()]
        print(f"   {want:<20}{len(hits):>3}" + (f"   e.g. {hits[0].get(c_org)} = {hits[0].get(c_temp)}" if hits else ""))

    if cand_rows:
        with open("yeast_growth_temperatures.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(cols), delimiter="\t"); w.writeheader(); w.writerows(cand_rows)
        print(f"\nwrote yeast_growth_temperatures.tsv ({len(cand_rows)} rows)")
    print("\nVERDICT: if 'yeast-like by genus' is under ~80 distinct species, the yeast-wide")
    print("approach is not worth building and we should say so.")

if __name__ == "__main__":
    main()
