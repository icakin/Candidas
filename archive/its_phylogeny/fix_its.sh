#!/bin/bash
# =============================================================================
# fix_its.sh - repair the two bad pulls from fetch_its.sh
#   * albicans  : got a whole chromosome (NC_032096) instead of ITS
#   * haemulonii: got var. vulneris instead of C. haemulonii sensu stricto
# Adds a sequence-length filter (200-1200 bp) so genome records can't match,
# and lists the top candidates so we can pick the right one.
# Run:  bash ~/Desktop/Projects/Candidas/phylo/fix_its.sh
# =============================================================================
set -u
OUT="$HOME/Desktop/Projects/Candidas/phylo/raw"
CAND="$HOME/Desktop/Projects/Candidas/phylo/candidates"
mkdir -p "$OUT" "$CAND"
EU="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

list_candidates () {
  local tag="$1"; local term="$2"
  echo ""
  echo "=== candidates for ${tag} ==="
  local ids
  ids=$(curl -s "${EU}/esearch.fcgi?db=nucleotide&term=${term}&retmax=6" \
        | grep -o "<Id>[0-9]*</Id>" | grep -o "[0-9]*")
  sleep 0.5
  if [ -z "${ids}" ]; then echo "  (no hits)"; return; fi
  : > "${CAND}/${tag}_candidates.fasta"
  for id in ${ids}; do
    curl -s "${EU}/efetch.fcgi?db=nucleotide&id=${id}&rettype=fasta&retmode=text" \
      >> "${CAND}/${tag}_candidates.fasta"
    sleep 0.4
  done
  grep "^>" "${CAND}/${tag}_candidates.fasta" | nl
}

# length-filtered ITS searches (200-1200 bp excludes genomes/chromosomes)
ALB='%22Candida%20albicans%22%5BOrganism%5D%20AND%20%22internal%20transcribed%20spacer%22%20AND%20200%3A1200%5BSLEN%5D%20AND%20refseq%5Bfilter%5D'
HAE='%22Candidozyma%20haemuli%22%5BOrganism%5D%20AND%20%22internal%20transcribed%20spacer%22%20AND%20200%3A1200%5BSLEN%5D'
HAE2='%22Candida%20haemulonii%22%5BOrganism%5D%20AND%20%22internal%20transcribed%20spacer%22%20AND%20200%3A1200%5BSLEN%5D%20AND%20sequence_from_type%5Bfilter%5D'

list_candidates albicans     "${ALB}"
list_candidates haemulonii   "${HAE}"
list_candidates haemulonii2  "${HAE2}"

echo ""
echo "Candidate FASTAs written to: ${CAND}"
echo "Review the headers above, then tell Claude which to use."
