#!/usr/bin/env python3
"""13_transfer_orthologs.py - build a draft GEM for a species that has NO curated
model (C. haemulonii, C. duobushaemulonii) by orthology transfer from a curated
template, instead of de-novo carving (CarveMe/gapseq are bacterially oriented).

    python3 gem/13_transfer_orthologs.py \\
        --template gem/models/auris_iRV973.xml \\
        --template_faa phylo/proteomes/auris_cladeI.faa \\
        --target_faa   phylo/proteomes/haemulonii.faa \\
        --out gem/models/haemulonii_draft.xml

Method: reciprocal best hits (DIAMOND) between template and target proteomes; a
template reaction is KEPT if every gene in its GPR has a target ortholog (or if it
has no GPR - spontaneous/transport - flagged). Reactions dropped for missing
orthologs are logged. Localization/compartments are inherited from the template.

Needs: cobra, diamond in PATH. Requires manual gap-filling afterwards - a dropped
reaction may be a real absence OR a missed ortholog; DO NOT claim species-unique
capability from this step alone (plan Step 1 / review).
"""
import argparse, subprocess, tempfile, os, collections
from pathlib import Path
import cobra
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

def rbh(faa_a, faa_b):
    """reciprocal best hits a<->b -> dict{gene_a: gene_b}"""
    td = Path(tempfile.mkdtemp())
    def dmnd(q, s, out):
        subprocess.run(["diamond","makedb","--in",s,"-d",str(td/"db"),"--quiet"],check=True)
        subprocess.run(["diamond","blastp","-q",q,"-d",str(td/"db"),"-o",out,
                        "--max-target-seqs","1","--evalue","1e-10","--quiet",
                        "--outfmt","6","qseqid","sseqid","bitscore"],check=True)
    a2b_f, b2a_f = str(td/"a2b.tsv"), str(td/"b2a.tsv")
    dmnd(faa_a, faa_b, a2b_f); dmnd(faa_b, faa_a, b2a_f)
    def best(f):
        d={}
        for line in open(f):
            q,s,bs=line.split("\t"); bs=float(bs)
            if q not in d or bs>d[q][1]: d[q]=(s,bs)
        return {k:v[0] for k,v in d.items()}
    ab, ba = best(a2b_f), best(b2a_f)
    return {q:s for q,s in ab.items() if ba.get(s)==q}

def gene_ids_in_faa(faa):
    return {l[1:].split()[0] for l in open(faa) if l.startswith(">")}

def main():
    p=argparse.ArgumentParser()
    for a in ("--template","--template_faa","--target_faa","--out"):
        p.add_argument(a, required=True)
    args=p.parse_args()

    tmpl = cobra.io.read_sbml_model(args.template)
    ortho = rbh(args.template_faa, args.target_faa)     # template_gene -> target_gene
    tgt_have = set(ortho)                                # template genes WITH a target ortholog

    kept, dropped, nogpr = [], [], []
    m = tmpl.copy()
    for r in list(m.reactions):
        genes = {g.id for g in r.genes}
        if not genes:
            nogpr.append(r.id); continue                # transport/spontaneous - keep, flag
        if genes <= tgt_have: kept.append(r.id)
        else: dropped.append(r.id); m.remove_reactions([r.id])

    m.id = Path(args.out).stem
    cobra.io.write_sbml_model(m, args.out)
    log = Path(args.out).with_suffix(".transfer.log")
    log.write_text(
        f"template={args.template}\ntarget={args.target_faa}\n"
        f"orthologs(RBH)={len(ortho)}\nreactions_kept={len(kept)}\n"
        f"reactions_dropped(missing ortholog)={len(dropped)}\n"
        f"reactions_no_GPR(kept, flagged)={len(nogpr)}\n\n"
        "DROPPED (verify: real absence vs missed ortholog before any claim):\n" +
        "\n".join(dropped) + "\n\nNO_GPR (spontaneous/transport, kept):\n" +
        "\n".join(nogpr) + "\n")
    print(f"wrote {args.out}")
    print(f"orthologs {len(ortho)} | kept {len(kept)} | dropped {len(dropped)} | "
          f"no-GPR {len(nogpr)}  (details in {log.name})")
    print("NEXT: gap-fill against YMS, then diff vs template. Do NOT claim unique "
          "capability from dropped reactions without gene/orthology/localization evidence.")

if __name__ == "__main__":
    main()
