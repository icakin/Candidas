"""gempaths.py -- the one place the etcGEM layer's directories are defined.

Every script under gem/ and gem/audits/ imports from here, so the layout can be
described once:

    gem/
      inputs/     hand-curated inputs: media maps, KO annotations, KEGG SMILES, KO->reaction list
      models/     the four SBML models (two curated, two drafts built by 06_build_drafts.py)
      tables/     everything a script writes: predictor outputs, kcat/MW tables, calibration,
                  counterfactual sweeps, audit results
      external/   third-party code and weights (DLKcat, Seq2Topt, ESM2, KOfam profiles, RefSeq
                  downloads). Not tracked by git; fetch with gem/fetch_external.sh
      notes/      the written record of each modelling decision and audit

Usage in a script at gem/ root:      from gempaths import *
Usage in gem/audits/:                 import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1])); from gempaths import *
CANDIDAS_ROOT, if set, overrides the repository root (for running against a copy).
"""
import os
from pathlib import Path

GEM = Path(__file__).resolve().parent
ROOT = Path(os.environ['CANDIDAS_ROOT']).resolve() if os.environ.get('CANDIDAS_ROOT') else GEM.parent
if os.environ.get('CANDIDAS_ROOT'):
    GEM = ROOT / 'gem'

INPUTS   = GEM / 'inputs'
MODELS   = GEM / 'models'
TABLES   = GEM / 'tables'
EXTERNAL = GEM / 'external'
NOTES    = GEM / 'notes'
RESULTS_TABLES  = ROOT / 'results' / 'tables'
RESULTS_FIGURES = ROOT / 'results' / 'figures'
PROTEOMES       = ROOT / 'phylo' / 'proteomes'

# the four species and their model + medium files, shared by every simulation script
SP = {'auris':            ('auris_iRV973_rekeyed.xml',    'medium_iRV973_auris.csv'),
      'haemulonii':       ('haemulonii_draft.xml',        'medium_iRV973_auris.csv'),
      'duobushaemulonii': ('duobushaemulonii_draft.xml',  'medium_iRV973_auris.csv'),
      'parapsilosis':     ('parapsilosis_iDC1003.xml',    'medium_iDC1003_parapsilosis.csv')}

# where each species' protein sequences come from, keyed by the gene ids the models use
PROTEOME = {'auris':            PROTEOMES / 'auris_cladeI.faa',        # B8441 RefSeq, XP_ ids
            'haemulonii':       PROTEOMES / 'haemulonii.faa',          # RefSeq, XP_ ids
            'duobushaemulonii': PROTEOMES / 'duobushaemulonii.faa',    # RefSeq, XP_ ids
            'parapsilosis':     INPUTS / 'parap_uniprot.tsv'}          # UniProt UP000005221; CPAR2_ locus tags

def load_proteome(sp):
    """gene id -> amino-acid sequence for one species (terminal '*' and 'X' removed).
    RefSeq FASTA is keyed by record id; the UniProt TSV by every name in its
    'Gene Names (ordered locus)' column, first entry winning."""
    import pandas as pd
    path = PROTEOME[sp]; out = {}
    if path.suffix == '.tsv':
        t = pd.read_csv(path, sep='\t', dtype=str).fillna('')
        for _, r in t.iterrows():
            for g in r['Gene Names (ordered locus)'].split():
                out.setdefault(g, r['Sequence'])
    else:
        from Bio import SeqIO
        for rec in SeqIO.parse(str(path), 'fasta'):
            out.setdefault(rec.id, str(rec.seq))
    return {k: v.rstrip('*').replace('X', '') for k, v in out.items()}

for _d in (TABLES,):
    _d.mkdir(parents=True, exist_ok=True)

__all__ = ['GEM', 'ROOT', 'INPUTS', 'MODELS', 'TABLES', 'EXTERNAL', 'NOTES',
           'RESULTS_TABLES', 'RESULTS_FIGURES', 'PROTEOMES', 'PROTEOME', 'SP', 'load_proteome']
