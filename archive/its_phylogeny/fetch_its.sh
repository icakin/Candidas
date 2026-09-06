#!/bin/bash
# =============================================================================
# fetch_its.sh - download curated ITS reference sequences (NCBI) for the
# Candidas phylogeny figure. Run from anywhere:
#     bash ~/Desktop/Projects/Candidas/phylo/fetch_its.sh
# Writes FASTA files into  ~/Desktop/Projects/Candidas/phylo/raw/
#
# Per taxon it tries, in order:
#   1) RefSeq Targeted-Loci ITS record (NR_..., curated, usually type strain)
#   2) any ITS sequence flagged "sequence from type"
#   3) any ITS sequence for the organism (last resort)
# and saves the first hit. Sources are NCBI E-utilities (esearch/efetch).
# =============================================================================
set -u
OUT="$HOME/Desktop/Projects/Candidas/phylo/raw"
mkdir -p "$OUT"
EU="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

fetch_one () {
  local tag="$1"; shift
  local ok=""
  for org in "$@"; do
    for extra in "%20AND%20refseq%5Bfilter%5D" "%20AND%20sequence_from_type%5Bfilter%5D" ""; do
      local term="%22${org// /%20}%22%5BOrganism%5D%20AND%20%22internal%20transcribed%20spacer%22${extra}"
      local ids
      ids=$(curl -s "${EU}/esearch.fcgi?db=nucleotide&term=${term}&retmax=1" \
            | grep -o "<Id>[0-9]*</Id>" | head -1 | grep -o "[0-9]*")
      sleep 0.5
      if [ -n "${ids}" ]; then
        curl -s "${EU}/efetch.fcgi?db=nucleotide&id=${ids}&rettype=fasta&retmode=text" \
          > "${OUT}/${tag}.fasta"
        sleep 0.5
        if grep -q "^>" "${OUT}/${tag}.fasta"; then
          echo "OK  ${tag}  <- ${org}  (uid ${ids})"
          head -1 "${OUT}/${tag}.fasta"
          ok=1; break 2
        fi
      fi
    done
  done
  [ -z "$ok" ] && echo "FAIL ${tag} - no ITS record found"
}

echo "Fetching ITS references into ${OUT}"
fetch_one auris            "Candidozyma auris" "Candida auris"
fetch_one haemulonii       "Candidozyma haemuli" "Candidozyma haemulonii" "Candida haemulonii"
fetch_one duobushaemulonii "Candidozyma duobushaemulonis" "Candidozyma duobushaemulonii" "Candida duobushaemulonii"
fetch_one pseudohaemulonii "Candidozyma pseudohaemuli" "Candida pseudohaemulonii"
fetch_one parapsilosis     "Candida parapsilosis"
fetch_one albicans         "Candida albicans"
fetch_one lusitaniae       "Clavispora lusitaniae"
fetch_one cerevisiae       "Saccharomyces cerevisiae"
echo "Done. Files in ${OUT}:"
ls -la "$OUT"
