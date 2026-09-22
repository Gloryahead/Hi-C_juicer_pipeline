#!/usr/bin/env bash
# build_all_sifs.sh
# Builds the Apptainer SIF container from its .def file.
# Run this ONCE before launching the Snakemake pipeline with --use-singularity.
# Requires: apptainer >= 1.2 (module load apptainer on most HPCs)
#
# Usage:
#   cd <repo_root>
#   bash containers/build_all_sifs.sh
#
# To build directly:
#   apptainer build containers/hic_juicer.sif containers/hic_juicer.def
#
# After the SIF is built, run the pipeline:
#   snakemake --use-singularity \
#             --singularity-args "--bind /xdisk,/groups" \
#             --profile workflow/profiles/slurm
#
# SIF size (approximate): hic_juicer.sif ~1.2 GB (toolchain) + juicer_tools.jar (~400 MB)

set -euo pipefail

CONTAINER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${CONTAINER_DIR}")"

module load apptainer 2>/dev/null || true  # HPC module; skip if not needed

build_sif() {
  local def="$1"
  local sif="$2"
  local name
  name=$(basename "${sif}" .sif)

  echo ""
  echo "═══════════════════════════════════════════════════"
  echo "  Building: ${name}.sif"
  echo "═══════════════════════════════════════════════════"

  if [[ -f "${sif}" ]]; then
    echo "  → Already exists; skipping. Remove to rebuild."
    return
  fi

  # Build from repo root so %files paths resolve correctly
  (cd "${REPO_ROOT}" && apptainer build "${sif}" "${def}")
  echo "  ✓ ${sif} built ($(du -sh "${sif}" | cut -f1))"
}

build_sif containers/hic_juicer.def containers/hic_juicer.sif

echo ""
echo "═══════════════════════════════════════════════════"
echo "  SIF built. Now run:"
echo "  snakemake --use-singularity \\"
echo "    --singularity-args \"--bind /xdisk,/groups\" \\"
echo "    --profile workflow/profiles/slurm"
echo "═══════════════════════════════════════════════════"
