#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Render the manuscript against the REGENERATED (runs/C1_reproduction) figures, without
# touching manuscript/draft/*.qmd and without writing into results/.
#
# _body_main.qmd and _body_supp.qmd reference figures as
#     ![](../../results/figures/manuscript/FIG1_decoupling.png)
# i.e. relative to manuscript/draft. So we build a STAGING COPY of the Quarto
# project whose sibling ../../results/figures/manuscript is a symlink to
# runs/C1_reproduction/figures/manuscript. Same .qmd bytes, different figures.
#
#   bash reports/tools/render_with_C1_figures.sh
#   -> runs/C1_reproduction/manuscript_build/draft/_output/{manuscript,supplementary,
#                                                 manuscript_combined}.{docx,pdf}
# ---------------------------------------------------------------------------
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RES="${CANDIDAS_RESULTS:-$ROOT/runs/C1_reproduction}"
STAGE="$RES/manuscript_build"

rm -rf "$STAGE"
mkdir -p "$STAGE/manuscript/draft" "$STAGE/results/figures"
cp "$ROOT"/manuscript/draft/*.qmd \
   "$ROOT"/manuscript/draft/_quarto.yml \
   "$ROOT"/manuscript/draft/header.tex \
   "$ROOT"/manuscript/draft/references.bib \
   "$ROOT"/manuscript/draft/nature-communications.csl "$STAGE/manuscript/draft/"
ln -s "$RES/figures/manuscript" "$STAGE/results/figures/manuscript"

echo "staging project : $STAGE/manuscript/draft"
echo "figures resolved: $(readlink "$STAGE/results/figures/manuscript")"
echo "figures present : $(ls "$RES/figures/manuscript" | wc -l | tr -d ' ')"
echo

# prove every referenced figure resolves BEFORE quarto silently drops it
miss=0
grep -ho '!\[\](\([^)]*\))' "$STAGE/manuscript/draft"/_body_*.qmd \
  | sed 's/.*(\(.*\)).*/\1/' | sed 's/{.*//' | sort -u | while read -r f; do
    p="$STAGE/manuscript/draft/$f"
    if [ -e "$p" ]; then echo "  OK      $f"; else echo "  MISSING $f"; miss=1; fi
  done
echo

cd "$STAGE/manuscript/draft" && quarto render
echo
echo "outputs:"; ls -la "$STAGE/manuscript/draft/_output"
