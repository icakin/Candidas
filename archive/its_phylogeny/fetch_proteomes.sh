#!/bin/bash
# =============================================================================
# fetch_proteomes.sh - download annotated proteomes (protein FASTA) for the
# phylogenomics + thermal-gene census.
#
#   Run on YOUR Mac (needs internet):
#     bash ~/Desktop/Projects/Candidas/phylo/fetch_proteomes.sh
#
# Resolves each strain via NCBI Assembly esearch/esummary (so we don't rely on
# hard-coded accessions), then pulls <asm>_protein.faa.gz from the RefSeq FTP
# path (falling back to GenBank). Writes to phylo/proteomes/.
# =============================================================================
set -u
OUT="$HOME/Desktop/Projects/Candidas/phylo/proteomes"
mkdir -p "$OUT"
EU="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

get_proteome () {   # get_proteome <tag> <esearch term>
  local tag="$1"; local term="$2"
  local uid ftp base url
  uid=$(curl -s "${EU}/esearch.fcgi?db=assembly&term=${term}&retmax=1" \
        | grep -o "<Id>[0-9]*</Id>" | head -1 | grep -o "[0-9]*")
  sleep 0.4
  if [ -z "${uid}" ]; then echo "FAIL ${tag}: no assembly hit"; return; fi

  local sum; sum=$(curl -s "${EU}/esummary.fcgi?db=assembly&id=${uid}")
  sleep 0.4
  ftp=$(echo "$sum" | grep -o "<FtpPath_RefSeq>[^<]*" | sed 's/<FtpPath_RefSeq>//' | head -1)
  [ -z "$ftp" ] && ftp=$(echo "$sum" | grep -o "<FtpPath_GenBank>[^<]*" | sed 's/<FtpPath_GenBank>//' | head -1)
  local name; name=$(echo "$sum" | grep -o "<AssemblyName>[^<]*" | sed 's/<AssemblyName>//' | head -1)
  local acc;  acc=$(echo "$sum"  | grep -o "<AssemblyAccession>[^<]*" | sed 's/<AssemblyAccession>//' | head -1)
  local org;  org=$(echo "$sum"  | grep -o "<Organism>[^<]*" | sed 's/<Organism>//' | head -1)

  if [ -z "$ftp" ]; then echo "FAIL ${tag}: no FTP path (${acc})"; return; fi
  ftp=${ftp/ftp:\/\//https://}
  base=$(basename "$ftp")
  url="${ftp}/${base}_protein.faa.gz"

  curl -s -f -o "${OUT}/${tag}.faa.gz" "$url"
  if [ $? -ne 0 ] || [ ! -s "${OUT}/${tag}.faa.gz" ]; then
    echo "FAIL ${tag}: no protein.faa.gz (${acc} ${name}) - assembly may be unannotated"
    rm -f "${OUT}/${tag}.faa.gz"; return
  fi
  gunzip -f "${OUT}/${tag}.faa.gz"
  local n; n=$(grep -c "^>" "${OUT}/${tag}.faa")
  printf "OK   %-22s %-18s %-34s %6d proteins\n" "$tag" "$acc" "$org" "$n"
}

echo "Downloading proteomes -> ${OUT}"
echo ""

# --- the four C. auris clade reference genomes (Muñoz et al. 2018) ---
get_proteome auris_cladeI   'B8441%5BAll%20Fields%5D%20AND%20auris%5BOrganism%5D'
get_proteome auris_cladeII  'B11220%5BAll%20Fields%5D%20AND%20auris%5BOrganism%5D'
get_proteome auris_cladeIII 'B11221%5BAll%20Fields%5D%20AND%20auris%5BOrganism%5D'
get_proteome auris_cladeIV  'B11245%5BAll%20Fields%5D%20AND%20auris%5BOrganism%5D'

# --- the haemulonii complex ---
get_proteome haemulonii       'haemulonii%5BOrganism%5D%20NOT%20duobus%5BAll%20Fields%5D%20NOT%20pseudo%5BAll%20Fields%5D%20AND%20latest%5Bfilter%5D'
get_proteome duobushaemulonii 'duobushaemulonii%5BOrganism%5D%20AND%20latest%5Bfilter%5D'
get_proteome pseudohaemulonii 'pseudohaemulonii%5BOrganism%5D%20AND%20latest%5Bfilter%5D'

# --- outgroup / context species ---
get_proteome parapsilosis 'CDC317%5BAll%20Fields%5D%20AND%20parapsilosis%5BOrganism%5D'
get_proteome albicans     'SC5314%5BAll%20Fields%5D%20AND%20%22Candida%20albicans%22%5BOrganism%5D%20AND%20latest%5Bfilter%5D'
get_proteome lusitaniae   '%22Clavispora%20lusitaniae%22%5BOrganism%5D%20AND%20latest%5Bfilter%5D'

echo ""
echo "=== RESULT ==="
ls -la "$OUT"
echo ""
echo "Total size:"; du -sh "$OUT"
