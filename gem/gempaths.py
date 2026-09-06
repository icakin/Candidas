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

for _d in (TABLES,):
    _d.mkdir(parents=True, exist_ok=True)

__all__ = ['GEM', 'ROOT', 'INPUTS', 'MODELS', 'TABLES', 'EXTERNAL', 'NOTES',
           'RESULTS_TABLES', 'RESULTS_FIGURES', 'PROTEOMES', 'SP']
