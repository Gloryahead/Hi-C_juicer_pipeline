#!/usr/bin/env bash
# run_pipeline.sh
# Convenience wrapper to launch the full Snakemake pipeline on SLURM.
# Run this after setup is complete (environments built, refs + Juicer installed).
#
# Usage (from repo root):
#   bash scripts/run_pipeline.sh              # full pipeline
#   bash scripts/run_pipeline.sh --dry-run    # check DAG without submitting
#   bash scripts/run_pipeline.sh --unlock     # unlock stale Snakemake directory
#   bash scripts/run_pipeline.sh --report     # generate HTML report
#   bash scripts/run_pipeline.sh --target results/hic1/aligned/inter_30.hic

set -euo pipefail

export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
export MAMBA_EXE=/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba
MC="${MAMBA_EXE}"

DRY_RUN=false
UNLOCK=false
REPORT=false
JOBS=100
TARGET=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run|-n)   DRY_RUN=true; shift ;;
    --unlock)       UNLOCK=true; shift ;;
    --report)       REPORT=true; shift ;;
    --target)       TARGET="$2"; shift 2 ;;
    *) shift ;;
  esac
done

[[ -f "workflow/Snakefile" ]]  || { echo "ERROR: run from the repo root (workflow/Snakefile not found)."; exit 1; }
[[ -f "config/config.yaml" ]]  || { echo "ERROR: config/config.yaml not found."; exit 1; }
[[ -f "config/samples.tsv" ]]  || { echo "ERROR: config/samples.tsv not found."; exit 1; }

if $UNLOCK; then
  "${MC}" run -n snakemake_env snakemake --unlock --snakefile workflow/Snakefile
  echo "Directory unlocked."
  exit 0
fi

CMD=(
  snakemake
    --snakefile workflow/Snakefile
    --profile   workflow/profiles/slurm
    --jobs      "${JOBS}"
    --rerun-incomplete
    --keep-going
    --latency-wait 60
    --printshellcmds
)

$DRY_RUN && CMD+=(--dry-run)
$REPORT  && CMD+=(--report results/snakemake_report.html)
[[ -n "${TARGET}" ]] && CMD+=("${TARGET}")

echo "=== Launching Hi-C (Juicer) pipeline ==="
echo "  Dry run : ${DRY_RUN}"
echo "  Target  : ${TARGET:-<rule all>}"
echo "  Command : ${CMD[*]}"
echo ""

"${MC}" run -n snakemake_env "${CMD[@]}"

if ! $DRY_RUN && ! $REPORT; then
  echo ""
  echo "=== Pipeline finished === (see results/qc/hic_qc_summary.txt and results/*/aligned/inter_30_qc.txt)"
fi
