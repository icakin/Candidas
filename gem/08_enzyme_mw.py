#!/usr/bin/env python3
"""08_enzyme_mw.py -- molecular weight of every enzyme in each species' model.

The sMOMENT protein-pool constraint charges each reaction MW_r / kcat_r per unit flux, so
every gene in every model needs a molecular weight. This takes the gene list from the
SBML model, looks each gene up in the species' proteome and computes the average-isotopic
molecular weight of the sequence (Biopython ProtParam, same as ExPASy Compute pI/Mw).

Inputs
    gem/models/<model>.xml                     gene ids in GPRs
    phylo/proteomes/auris_cladeI.faa           C. auris B8441 RefSeq (XP_ ids)
    phylo/proteomes/haemulonii.faa             C. haemulonii RefSeq (XP_ ids)
    phylo/proteomes/duobushaemulonii.faa       C. duobushaemulonii RefSeq (XP_ ids)
    gem/inputs/parap_uniprot.tsv               C. parapsilosis CDC317 UniProt proteome
                                               (UP000005221) as a table: Entry, ordered-locus
                                               gene name (CPAR2_ tags, the model's gene ids),
                                               Sequence. The FASTA of the same proteome lacks
                                               a GN= field for 149 model genes, so the TSV is
                                               the source.
Outputs
    gem/tables/enzyme_mw_<species>.csv         gene, length_aa, MW_kDa

Run:  python3 gem/08_enzyme_mw.py            (~10 s; needs cobra, biopython)
Genes without a sequence are listed on stderr and omitted; the modelling scripts treat a
missing MW as the species median, so this list should be empty and is checked.
"""
import sys
import pandas as pd
import cobra
from Bio.SeqUtils import molecular_weight
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP, PROTEOMES

def main():
    for sp, (xml, _) in SP.items():
        model = cobra.io.read_sbml_model(str(MODELS / xml))
        seqs = load_proteome(sp)
        rows, missing = [], []
        for g in sorted(model.genes, key=lambda x: x.id):
            s = seqs.get(g.id)
            if s is None:
                missing.append(g.id); continue
            rows.append(dict(gene=g.id, length_aa=len(s),
                             MW_kDa=round(molecular_weight(s, 'protein') / 1000, 2)))
        # keep the model's own gene order rather than sorted, to match the original tables
        order = {g.id: i for i, g in enumerate(model.genes)}
        df = pd.DataFrame(rows).sort_values('gene', key=lambda c: c.map(order)).reset_index(drop=True)
        out = TABLES / f'enzyme_mw_{sp}.csv'
        df.to_csv(out, index=False)
        print(f'{sp:18s} {len(model.genes):4d} genes, {len(df):4d} with sequence -> {out.name}')
        if missing:
            print(f'   no sequence for {len(missing)} gene(s): {missing[:8]}{" ..." if len(missing) > 8 else ""}',
                  file=sys.stderr)

if __name__ == '__main__':
    main()
