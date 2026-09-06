#!/bin/bash
# =============================================================================
# 10_fetch_curated_models.sh - get the two curated Candida GEMs. Run on the Mac.
#   bash gem/10_fetch_curated_models.sh
#
# Verified sources (2026-09):
#  * C. parapsilosis  iDC1003 / Viana2022  -> BioModels MODEL2404250002
#      paper: doi:10.3390/genes13020303  (PMID 35205348)   [directly downloadable]
#  * C. auris         iRV973              -> FEMS Yeast Res doi:10.1093/femsyr/foad045
#      NOT in BioModels; get the SBML from the paper's supplementary data (or the
#      authors). No stable direct URL - hence the manual step below.
#
# Do NOT substitute a CarveMe/gapseq draft for these (bacterially oriented).
# =============================================================================
set -u
DIR="$(cd "$(dirname "$0")" && pwd)/models"
mkdir -p "$DIR"

echo "== C. parapsilosis (BioModels MODEL2404250002) =="
# BioModels REST download endpoint (returns the model file(s)):
if curl -fsSL "https://www.ebi.ac.uk/biomodels/model/download/MODEL2404250002" \
        -o "$DIR/parapsilosis_iDC1003.download" 2>/dev/null; then
  # response may be a zip or an xml; detect and normalise
  if file "$DIR/parapsilosis_iDC1003.download" | grep -qi zip; then
    unzip -o "$DIR/parapsilosis_iDC1003.download" -d "$DIR/parap_tmp" >/dev/null
    x=$(find "$DIR/parap_tmp" -iname '*.xml' | head -1)
    [ -n "$x" ] && cp "$x" "$DIR/parapsilosis_iDC1003.xml"
    rm -rf "$DIR/parap_tmp"
  else
    mv "$DIR/parapsilosis_iDC1003.download" "$DIR/parapsilosis_iDC1003.xml"
  fi
  rm -f "$DIR/parapsilosis_iDC1003.download"
  echo "  -> $DIR/parapsilosis_iDC1003.xml"
else
  echo "  automated pull failed - open https://www.ebi.ac.uk/biomodels/MODEL2404250002"
  echo "  and download the SBML by hand to $DIR/parapsilosis_iDC1003.xml"
fi

echo
echo "== C. auris iRV973 (manual: journal supplement) =="
echo "  https://doi.org/10.1093/femsyr/foad045  -> supplementary SBML"
echo "  save as  $DIR/auris_iRV973.xml"
echo
echo "(optional reference) C. albicans iRV781 -> $DIR/albicans_iRV781.xml"
echo
echo "Then:  python3 gem/11_inspect_models.py"
