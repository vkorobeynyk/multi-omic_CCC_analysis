#!/usr/bin/env bash
# Render the whole pipeline in order, stopping at the first failure.
#
# Usage (from anywhere):
#   ./render_all.sh          # render everything, 01 -> 11
#   ./render_all.sh 08       # start at 08 (skip 01-07), e.g. after fixing a later script
#
# .qmd files are rendered with Quarto (R). .ipynb notebooks are executed with
# `quarto render --execute`, using the kernel saved in each notebook
# (lianaPlus_localenv), so that kernel must be registered with Jupyter.
# Each render's full log goes to render_logs/<script>.log.

set -euo pipefail

cd "$(dirname "$0")"   # run from the project folder so relative output/ paths work

SCRIPTS=(
  01_processing_subset_to_HPC.qmd
  02_Integration_clustering.qmd
  03a_RCTD_spotdeconvolution.qmd
  03b_runLiana_scRNAseq.ipynb
  03c_crossModality_plots.qmd
  04_KEGG_PathwayConvergence.qmd
  05_STRING_NetworkAnalysis.qmd
  06_Synthesis_FunctionBySpatialPair.qmd
  07_STARmapPLUS_ControlSpatialObject.qmd
  08_STARmapPLUS_CCC_GenePanelCheck.qmd
  09_STARmapPLUS_BivariateMoransI.ipynb
  10_STARmapPLUS_BivariateMoransI_AgeComparison.ipynb
  11_Proteomics_LR_AgeTrends.qmd
)

START_FROM="${1:-}"
LOG_DIR="render_logs"
mkdir -p "$LOG_DIR"

started=false
[[ -z "$START_FROM" ]] && started=true
total_start=$(date +%s)

for f in "${SCRIPTS[@]}"; do
  if ! $started; then
    [[ "$f" == "$START_FROM"* ]] && started=true || continue
  fi

  if [[ ! -f "$f" ]]; then
    echo "!! $f not found -- stopping." >&2
    exit 1
  fi

  log="$LOG_DIR/${f%.*}.log"
  echo "==> [$(date +%H:%M:%S)] rendering $f  (log: $log)"
  t0=$(date +%s)

  if [[ "$f" == *.ipynb ]]; then
    # Notebooks carry no Quarto YAML, so embed figures into the .html here (the .qmd
    # files already set embed-resources: true) -- otherwise Quarto writes a <name>_files/ folder
    cmd=(quarto render "$f" --execute -M embed-resources:true)
  else
    cmd=(quarto render "$f")
  fi

  if ! "${cmd[@]}" > "$log" 2>&1; then
    echo "!! $f FAILED after $(( $(date +%s) - t0 ))s. Last lines of the log:" >&2
    tail -n 30 "$log" >&2
    echo "Fix it, then resume with: ./render_all.sh ${f%%_*}" >&2
    exit 1
  fi

  echo "    done in $(( $(date +%s) - t0 ))s"
done

if ! $started; then
  echo "!! No script matches '$START_FROM'." >&2
  exit 1
fi

echo "All done in $(( ($(date +%s) - total_start) / 60 )) min."
