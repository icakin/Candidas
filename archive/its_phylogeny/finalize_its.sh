#!/bin/bash
# =============================================================================
# finalize_its.sh - fetch the two corrected records by exact accession and
# reorganise raw/ into the final 9-taxon set for the phylogeny.
# Run: bash ~/Desktop/Projects/Candidas/phylo/finalize_its.sh
# =============================================================================
set -u
OUT="$HOME/Desktop/Projects/Candidas/phylo/raw"
EU="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
mkdir -p "$OUT"

get_acc () {  # get_acc <accession> <outfile>
  curl -s "${EU}/efetch.fcgi?db=nucleotide&id=${1}&rettype=fasta&retmode=text" > "${OUT}/${2}.fasta"
  sleep 0.5
  if grep -q "^>" "${OUT}/${2}.fasta"; then
    echo "OK  ${2}"; head -1 "${OUT}/${2}.fasta"
  else
    echo "FAIL ${2} (${1})"
  fi
}

# corrected records
get_acc NR_125332.1 albicans            # C. albicans CBS 562, type (was a chromosome)
get_acc NR_202147.1 haemulonii_vulnera  # C. haemuli var. vulneris, type (context tip)
get_acc NR_130669.1 haemulonii          # C. haemuli CBS 5149, type  <-- sensu stricto

echo ""
echo "=== FINAL SET (${OUT}) ==="
for f in "${OUT}"/*.fasta; do
  printf "%-26s %6s bp  " "$(basename "$f")" "$(grep -v '^>' "$f" | tr -d '\n' | wc -c | tr -d ' ')"
  head -1 "$f"
done
