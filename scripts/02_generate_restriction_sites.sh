#!/usr/bin/env bash
# 02_generate_restriction_sites.sh
# Generates the restriction-site position file Juicer needs for fragment-level
# deduplication and normalization, using Juicer's own misc/generate_site_positions.py
# (the tutorial's hand-rolled version of this script is not used here — this
# repo prefers the upstream, field-tested implementation; see docs/KEY_CONCEPTS.md).
#
# Requires scripts/00_download_refs.sh and scripts/01_setup_juicer.sh to have
# already run (needs both hg38.fa and the cloned juicer repo).
#
# Usage:
#   REF_DIR=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline/ref \
#   TOOLS_DIR=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline/tools \
#   ENZYME=MboI \
#     bash scripts/02_generate_restriction_sites.sh

set -euo pipefail

REF_DIR="${REF_DIR:-ref}"
TOOLS_DIR="${TOOLS_DIR:-tools}"
ENZYME="${ENZYME:-MboI}"
GENOME_ID="${GENOME_ID:-hg38}"

GEN_SCRIPT="${TOOLS_DIR}/juicer/misc/generate_site_positions.py"
[[ -f "${GEN_SCRIPT}" ]] || { echo "ERROR: ${GEN_SCRIPT} not found — run scripts/01_setup_juicer.sh first."; exit 1; }
[[ -f "${REF_DIR}/hg38.fa" ]] || { echo "ERROR: ${REF_DIR}/hg38.fa not found — run scripts/00_download_refs.sh first."; exit 1; }

OUT_FILE="${REF_DIR}/${GENOME_ID}_${ENZYME}.txt"
if [[ -f "${OUT_FILE}" ]]; then
  echo "  → ${OUT_FILE} already present, skipping"
  exit 0
fi

cd "${REF_DIR}"
python3 "${GEN_SCRIPT}" "${ENZYME}" "${GENOME_ID}" "$(pwd)/hg38.fa"

# generate_site_positions.py writes <genomeID>_<enzyme>.txt into the cwd —
# fail loudly rather than silently pointing config.yaml at a missing file.
[[ -f "${GENOME_ID}_${ENZYME}.txt" ]] || { echo "ERROR: expected output ${GENOME_ID}_${ENZYME}.txt was not created."; exit 1; }

echo "Done: ${OUT_FILE}"
echo "Point config/config.yaml's genome.restriction_sites at this path if you used custom REF_DIR/ENZYME/GENOME_ID."
