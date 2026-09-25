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

# A real SRR1658570 mate file is many GB gzipped (~202M read pairs); a
# leftover from an interrupted/failed prefetch (a partial or 0-byte stub)
# is nowhere close. Checking existence alone let exactly that kind of
# leftover silently "skip" every retry forever — hit this for real after a
# prefetch TLS failure still left small files behind. Require a real size
# AND a valid gzip stream (catches a truncated-but-large-enough file too)
# before trusting a file as complete; anything that fails either check gets
# removed so the download below runs for real instead of skipping.
MIN_BYTES=$((1 * 1024 * 1024 * 1024))   # 1 GB floor — real files are ~10-20 GB each
_looks_complete() {
  local f="$1"
  [[ -f "${f}" ]] || return 1
  local size
  size=$(stat -c%s "${f}" 2>/dev/null || stat -f%z "${f}" 2>/dev/null || echo 0)
  [[ "${size}" -ge "${MIN_BYTES}" ]] || return 1
  gzip -t "${f}" 2>/dev/null || return 1
  return 0
}

if _looks_complete "${ACCESSION}_1.fastq.gz" && _looks_complete "${ACCESSION}_2.fastq.gz"; then
  echo "  → ${ACCESSION} FASTQs already present and look complete, skipping"
  exit 0
elif [[ -f "${ACCESSION}_1.fastq.gz" || -f "${ACCESSION}_2.fastq.gz" || -d "${ACCESSION}" ]]; then
  # ${ACCESSION}/ is prefetch's OWN working directory (${ACCESSION}/${ACCESSION}.sra) —
  # a run that fails before reaching the final `rm -rf "${ACCESSION}"` below leaves
  # it behind, and prefetch will silently reuse whatever partial/stale .sra is in
  # there on the next attempt ("is found locally") instead of fetching fresh. Hit
  # this for real: a stale dir here is exactly what let a 25GB accession quietly
  # resolve to a stale local copy once prefetch's own size cap skipped the real
  # fresh fetch. Must go, not just the final fastq(.gz) outputs.
  echo "  → found an incomplete/corrupt leftover from a prior failed run — removing and re-downloading"
  rm -f "${ACCESSION}_1.fastq.gz" "${ACCESSION}_2.fastq.gz" "${ACCESSION}_1.fastq" "${ACCESSION}_2.fastq"
  rm -rf "${ACCESSION}"
fi

echo "  ↓ prefetch ${ACCESSION}"
# --max-size: prefetch's default cap is 20GB; SRR1658570's real file is ~25GB and
# gets silently "skipped" (falling back to whatever's cached, if anything) without
# this — hit that for real. --force yes: ignore/overwrite anything prefetch itself
# finds already cached, so behavior doesn't depend on what a prior failed attempt
# left behind (the cleanup above already handles OUR side of that).
prefetch --max-size 30g --force yes "${ACCESSION}"

echo "  → fasterq-dump (paired)"
fasterq-dump "${ACCESSION}" --split-files --threads 8

gzip "${ACCESSION}_1.fastq" "${ACCESSION}_2.fastq"
rm -rf "${ACCESSION}"   # prefetch's .sra cache dir, no longer needed once dumped

echo ""
echo "Done: ${DATA_DIR}/${ACCESSION}_1.fastq.gz, ${DATA_DIR}/${ACCESSION}_2.fastq.gz"
echo "Now fill in config/samples.tsv's hic1 row: fastq_r1/fastq_r2 -> these two paths."
