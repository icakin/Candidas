#!/usr/bin/env python3
"""fetch_kegg_smiles.py - get SMILES for the KEGG compounds the DLKcat run needs.
Run on the Mac (internet):  python3 gem/fetch_kegg_smiles.py
Reads gem/kegg_compounds_needed.txt, writes gem/kegg_smiles.tsv (kegg<TAB>SMILES).
Needs rdkit (conda install -c conda-forge rdkit) + requests. ~1000 compounds, ~10 min.
"""
import time, requests
from pathlib import Path
from rdkit import Chem
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
GEM=Path(__file__).resolve().parent
cids=[l.strip() for l in open(INPUTS / "kegg_compounds_needed.txt") if l.strip()]
out=open(INPUTS / "kegg_smiles.tsv","w"); ok=miss=0
for i,c in enumerate(cids,1):
    try:
        mol=requests.get(f"https://rest.kegg.jp/get/cpd:{c}/mol",timeout=30).text
        m=Chem.MolFromMolBlock(mol)
        if m: out.write(f"{c}\t{Chem.MolToSmiles(m)}\n"); ok+=1
        else: miss+=1
    except Exception: miss+=1
    if i%100==0: print(f"  {i}/{len(cids)}  ok={ok} miss={miss}",flush=True); out.flush()
    time.sleep(0.15)
out.close(); print(f"done: {ok} SMILES, {miss} missing -> gem/kegg_smiles.tsv")
