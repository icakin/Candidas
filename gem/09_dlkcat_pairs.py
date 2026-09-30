#!/usr/bin/env python3
"""09_dlkcat_pairs.py -- every (reaction, enzyme, substrate) triple DLKcat must predict a kcat for.

DLKcat (Li et al. 2022, Nat Catal 5:662) predicts kcat from one protein sequence and one
substrate SMILES. For each species this walks the model, takes every reaction that has a
gene-protein-reaction rule, every gene in that rule, and every reactant (left-hand side,
stoichiometric coefficient < 0) that carries a KEGG compound annotation and is not a
currency metabolite. Reversible reactions contribute their left-hand side only: the models
are written in the KEGG forward direction and the kcat is charged per unit flux whichever
way it runs.

Currency metabolites are excluded because DLKcat was trained on the reaction's principal
substrate and its estimate for ATP or NADH as "the substrate" is not meaningful:
    water, ATP/ADP/AMP, NAD(P)(H), O2, Pi, PPi, CoA, CO2, NH3, UTP/UDP/UMP, FAD, SAM,
    GTP/GDP, H+, IDP, CDP, dATP, dADP
This list was fixed before any kcat was looked at and is the one the committed tables use
(regenerating them from this script gives the identical 4042/3974/3871/4667 triples).

Outputs
    gem/tables/dlkcat_pairs_<species>.tsv     reaction, gene, substrate_kegg
    gem/inputs/kegg_compounds_needed.txt      the union of substrates, for 12_fetch_kegg_smiles.py
    gem/tables/dlkcat_input_<species>.tsv     DLKcat's input format (Substrate Name, Substrate
                                              SMILES, Protein Sequence), one row per triple that
                                              has a SMILES; Substrate Name is gene|reaction|kegg
                                              so 11_aggregate_kcat.py can read the triple back.
                                              Written only when gem/inputs/kegg_smiles.tsv exists.
Run:  python3 gem/09_dlkcat_pairs.py        then 12_fetch_kegg_smiles.py (once), then re-run
      this to write the DLKcat inputs, then 10_run_dlkcat.sh.
"""
import pandas as pd
import cobra
from gempaths import *  # GEM, INPUTS, MODELS, TABLES, SP, load_proteome

CURRENCY = {'C00001', 'C00002', 'C00003', 'C00004', 'C00005', 'C00006', 'C00007', 'C00008',
            'C00009', 'C00010', 'C00011', 'C00013', 'C00014', 'C00015', 'C00016', 'C00019',
            'C00020', 'C00035', 'C00044', 'C00075', 'C00080', 'C00104', 'C00105', 'C00112',
            'C00131', 'C00206'}

def kegg_id(met):
    k = met.annotation.get('kegg.compound')
    return k[0] if isinstance(k, list) else k

def triples(model):
    rows = []
    for r in model.reactions:
        if not r.gene_reaction_rule:
            continue
        subs = sorted({kegg_id(m) for m, c in r.metabolites.items()
                       if c < 0 and kegg_id(m) and kegg_id(m) not in CURRENCY})
        for g in sorted(x.id for x in r.genes):
            for k in subs:
                rows.append((r.id, g, k))
    return pd.DataFrame(rows, columns=['reaction', 'gene', 'substrate_kegg'])

def main():
    smiles_f = INPUTS / 'kegg_smiles.tsv'
    smiles = (dict(pd.read_csv(smiles_f, sep='\t', header=None).values) if smiles_f.exists() else None)
    needed = set()
    for sp, (xml, _) in SP.items():
        model = cobra.io.read_sbml_model(str(MODELS / xml))
        df = triples(model)
        df.to_csv(TABLES / f'dlkcat_pairs_{sp}.tsv', sep='\t', index=False)
        needed |= set(df.substrate_kegg)
        msg = f'{sp:18s} {df.reaction.nunique():5d} reactions  {df.gene.nunique():4d} genes  {len(df):5d} triples'
        if smiles is not None:
            seqs = load_proteome(sp)
            d = df[df.substrate_kegg.isin(smiles)].copy()
            d['Substrate Name']   = d.gene + '|' + d.reaction + '|' + d.substrate_kegg
            d['Substrate SMILES'] = d.substrate_kegg.map(smiles)
            d['Protein Sequence'] = d.gene.map(seqs)
            d = d.dropna(subset=['Protein Sequence'])
            d[['Substrate Name', 'Substrate SMILES', 'Protein Sequence']].to_csv(
                TABLES / f'dlkcat_input_{sp}.tsv', sep='\t', index=False)
            msg += f'  -> dlkcat_input_{sp}.tsv ({len(d)} rows with SMILES)'
        print(msg)
    with open(INPUTS / 'kegg_compounds_needed.txt', 'w') as fh:
        fh.write('\n'.join(sorted(needed)) + '\n')
    print(f'{len(needed)} distinct substrates -> inputs/kegg_compounds_needed.txt')

if __name__ == '__main__':
    main()
