#!/usr/bin/env bash
# =============================================================================
# fetch_external.sh -- third-party code, weights and databases the etcGEM layer needs
# =============================================================================
#   bash gem/fetch_external.sh            everything (~9 GB, most of it KOfam profiles)
#   bash gem/fetch_external.sh --no-kofam skip KOfam (only needed to redo 04_kofam_annotate.sh;
#                                         the KO tables it produces are tracked in gem/inputs/kofam/)
#
# Everything lands in gem/external/, which git ignores. What is fetched and why:
#
#   DLKcat/            SysBioChalmers/DLKcat (kcat predictor, trained weights included in the
#                      repo under DeeplearningApproach/Results/output/). Used by 10_run_dlkcat.sh.
#   Seq2Topt/          SizheQiu/Seq2Topt v1.0.0 (Topt and Tm predictors). Used by 14_ and 15_.
#   large_model_pth/   the two Seq2Topt release checkpoints:
#                        model_topt_window=3_r2=0.57.pth   model_tm_window=3_r2=0.76.pth
#   torch_home/        ESM2-t6-8M weights (esm2_t6_8M_UR50D.pt + contact-regression head), in
#                      the layout torch.hub expects; 14_run_seq2topt_seq2tm.sh sets TORCH_HOME here.
#   kofam/             KOfam HMM profiles + ko_list (ftp.genome.jp/pub/db/kofam). Used by 04_.
#   b8441_refseq/      NCBI datasets download of GCF_002759435.1 (C. auris B8441), GFF only;
#                      the protein FASTA the pipeline uses is tracked at phylo/proteomes/.
#
# Curated models (iRV973, iDC1003) are fetched by 01_fetch_curated_models.sh, not here.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXT="$HERE/external"; mkdir -p "$EXT"; cd "$EXT"
NOKOFAM=0; [ "${1:-}" = "--no-kofam" ] && NOKOFAM=1

[ -d DLKcat ]   || git clone --depth 1 https://github.com/SysBioChalmers/DLKcat.git
[ -d Seq2Topt ] || git clone --depth 1 --branch v1.0.0 https://github.com/SizheQiu/Seq2Topt.git

# release assets are named model_tm_window.3_r2.0.76.pth; the upstream code loads them as
# model_tm_window=3_r2=0.76.pth, so save under the name the code expects
mkdir -p large_model_pth
for f in "model_topt_window=3_r2=0.57.pth" "model_tm_window=3_r2=0.76.pth"; do
  asset="$(echo "$f" | sed 's/window=3_r2=/window.3_r2./')"
  [ -f "large_model_pth/$f" ] || curl -fL -o "large_model_pth/$f" \
      "https://github.com/SizheQiu/Seq2Topt/releases/download/v1.0.0/$asset"
done

mkdir -p torch_home/hub/checkpoints
for f in esm2_t6_8M_UR50D.pt esm2_t6_8M_UR50D-contact-regression.pt; do
  [ -f "torch_home/hub/checkpoints/$f" ] || curl -fL -o "torch_home/hub/checkpoints/$f" \
      "https://dl.fbaipublicfiles.com/fair-esm/models/$f"
done

if [ "$NOKOFAM" = 0 ]; then
  mkdir -p kofam && cd kofam
  [ -f ko_list ]      || { curl -fL -o ko_list.gz https://www.genome.jp/ftp/db/kofam/ko_list.gz && gunzip ko_list.gz; }
  [ -d profiles ]     || { curl -fL -o profiles.tar.gz https://www.genome.jp/ftp/db/kofam/profiles.tar.gz && tar xzf profiles.tar.gz; }
  cd ..
fi

echo "external/ ready:"; du -sh "$EXT"/* 2>/dev/null
