#!/bin/bash
# --------------------------------------------------------------
### One-off recovery: re-run only Juicer's "final" stage (hic maps,
### normalization, statistics, Arrowhead) on an existing aligned/ directory.
### Use after a run that completed alignment and de-duplication but lost or
### mis-built the .hic files (job 24212676: addNorm was OOM-killed because the
### job had less memory than juicer.sh's hard-coded 150000 MB Java heap, and
### Snakemake then deleted inter_30.hic over a wrongly declared output file).
### Takes about 1-2 h instead of repeating the ~14 h alignment.
### Run it by hand; do NOT use it while Snakemake is also running this sample.
# --------------------------------------------------------------
#SBATCH --job-name=hic_juicer_final
#SBATCH --output=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline/logs/slurm/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline/logs/slurm/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=40            # Puma gives 5 GB per CPU: 40 CPUs = 200 GB, above the 150000 MB heap
#SBATCH --time=12:00:00

set -euo pipefail

SAMPLE="${SAMPLE:-hic1}"
DATA_BASE=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline
TOPDIR="${DATA_BASE}/results/${SAMPLE}"

export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
eval "$(/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba shell hook --shell bash)"
micromamba activate hic_juicer_env

echo "========== Juicer final stage started at: $(date) =========="

for f in merged_dedup.bam merged1.txt merged30.txt; do
  [[ -s "${TOPDIR}/aligned/${f}" ]] || { echo "ERROR: ${TOPDIR}/aligned/${f} missing; the final stage needs the aligned/ output of an earlier run."; exit 1; }
done

# Run from xdisk: juicer's samtools/java temporaries use the current directory.
cd "${TOPDIR}"

bash "${DATA_BASE}/tools/juicer/scripts/common/juicer.sh" \
    -S final \
    -D "${DATA_BASE}/tools/juicer" \
    -d "${TOPDIR}" \
    -g hg38 \
    -s MboI \
    -p "${DATA_BASE}/ref/hg38.chrom.sizes" \
    -y "${DATA_BASE}/ref/hg38_MboI.txt" \
    -z "${DATA_BASE}/ref/hg38.fa" \
    -t 16

echo "========== Juicer final stage finished at: $(date) =========="
ls -la "${TOPDIR}/aligned/inter_30.hic" "${TOPDIR}/aligned/inter_30.txt"
echo "Next: bash scripts/run_pipeline.sh --dry-run (juicer_pipeline should no longer be listed), then sbatch slurm/03_pipeline.slurm."
