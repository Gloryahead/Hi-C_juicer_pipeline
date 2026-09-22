#!/usr/bin/env bash
# 00_download_refs.sh
# Downloads the hg38 reference FASTA (UCSC, same source as the tutorial) and
# builds the BWA index + chrom.sizes file Juicer needs. ~1 GB download,
# ~1 hour for the BWA index; one-time.
#
# Usage:
#   REF_DIR=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline/ref bash scripts/00_download_refs.sh

set -euo pipefail

REF_DIR="${REF_DIR:-ref}"
mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

if [[ -f hg38.fa ]]; then
  echo "  → hg38.fa already present, skipping download"
else
  echo "  ↓ hg38.fa.gz (UCSC)"
  wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  gunzip hg38.fa.gz
fi

echo "Building BWA index (one-time, ~1 hour)..."
if [[ ! -f hg38.fa.bwt ]]; then
  bwa index hg38.fa
else
  echo "  → BWA index already present, skipping"
fi

echo "Indexing FASTA + writing chrom.sizes..."
samtools faidx hg38.fa
cut -f1,2 hg38.fa.fai > hg38.chrom.sizes

echo ""
echo "Done. Point config/config.yaml's genome.* paths at ${REF_DIR}/ if you used a custom REF_DIR."
echo "Next: scripts/01_setup_juicer.sh, then scripts/02_generate_restriction_sites.sh"
