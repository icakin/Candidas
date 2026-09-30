#!/usr/bin/env bash
# =============================================================================
# 10_run_dlkcat.sh -- predict kcat for every (enzyme, substrate) pair with DLKcat
# =============================================================================
#   bash gem/10_run_dlkcat.sh                 all four species (~20-40 min on a laptop CPU)
#   bash gem/10_run_dlkcat.sh haemulonii      one species
#
# DLKcat: Li, F. et al. (2022) Deep learning-based kcat prediction enables improved
# enzyme-constrained model reconstruction. Nat Catal 5:662-672. Code and the trained
# model come from https://github.com/SysBioChalmers/DLKcat (gem/fetch_external.sh clones
# it to gem/external/DLKcat). The trained weights ship in the repository under
# DeeplearningApproach/Results/output/, and the token dictionaries in Data/input.zip.
#
# Input   gem/tables/dlkcat_input_<species>.tsv     from 09_dlkcat_pairs.py
# Output  gem/tables/kcat_<species>.tsv             DLKcat's own output format
#           Substrate Name (gene|reaction|kegg), Substrate SMILES, Protein Sequence,
#           Kcat value (1/s)
#
# The upstream script hard-codes '../../Data/input/*.pickle', '../../Results/output/...'
# and './output.tsv', so it has to run from DeeplearningApproach/Code/example. Pairs whose
# SMILES contain a '.' (salts) are skipped by DLKcat itself; the model does not see them.
#
# Environment: python with torch, rdkit, numpy, requests (gem/requirements.txt). CPU is
# fine; the model is small. The prediction is deterministic for a given checkpoint.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DLK="$HERE/external/DLKcat/DeeplearningApproach"
TAB="$HERE/tables"
PY="${PYTHON:-python3}"

[ -d "$DLK" ] || { echo "!! $DLK not found. Run gem/fetch_external.sh first."; exit 1; }
if [ ! -d "$DLK/Data/input" ]; then
  echo "unzipping DLKcat token dictionaries (Data/input.zip)"
  (cd "$DLK/Data" && unzip -q input.zip)
fi

SPECIES=("$@"); [ ${#SPECIES[@]} -eq 0 ] && SPECIES=(auris haemulonii duobushaemulonii parapsilosis)
for sp in "${SPECIES[@]}"; do
  in="$TAB/dlkcat_input_${sp}.tsv"; out="$TAB/kcat_${sp}.tsv"
  [ -f "$in" ] || { echo "!! $in missing: run 09_dlkcat_pairs.py after 12_fetch_kegg_smiles.py"; exit 1; }
  echo "== $sp: $(($(wc -l < "$in") - 1)) pairs"
  ( cd "$DLK/Code/example" && "$PY" prediction_for_input.py "$in" > "$TAB/dlkcat_${sp}.log" 2>&1 \
      && mv output.tsv "$out" )
  echo "   -> $out ($(($(wc -l < "$out") - 1)) predictions)"
done
echo "next: python3 gem/11_aggregate_kcat.py"
