#!/usr/bin/env bash
# 03_download_test_data.sh
# Downloads the tutorial's own worked-example dataset: SRR1658570 (see
# data/accessions/01_hic_test_data.txt for provenance) as paired FASTQs
# named per config/samples.tsv's "hic1" row.
#
# Usage:
#   DATA_DIR=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline/data/fastq \
#     bash scripts/03_download_test_data.sh

set -euo pipefail

DATA_DIR="${DATA_DIR:-data/fastq}"
ACCESSION="${ACCESSION:-SRR1658570}"

mkdir -p "${DATA_DIR}"
cd "${DATA_DIR}"

if [[ -f "${ACCESSION}_1.fastq.gz" && -f "${ACCESSION}_2.fastq.gz" ]]; then
  echo "  → ${ACCESSION} FASTQs already present, skipping"
  exit 0
fi

echo "  ↓ prefetch ${ACCESSION}"
prefetch "${ACCESSION}"

echo "  → fasterq-dump (paired)"
fasterq-dump "${ACCESSION}" --split-files --threads 8

gzip "${ACCESSION}_1.fastq" "${ACCESSION}_2.fastq"
rm -rf "${ACCESSION}"   # prefetch's .sra cache dir, no longer needed once dumped

echo ""
echo "Done: ${DATA_DIR}/${ACCESSION}_1.fastq.gz, ${DATA_DIR}/${ACCESSION}_2.fastq.gz"
echo "Now fill in config/samples.tsv's hic1 row: fastq_r1/fastq_r2 -> these two paths."
