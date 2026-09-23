#!/bin/bash
# --------------------------------------------------------------
### Download the tutorial's test dataset (SRR1658570, ~202M read
### pairs / ~40GB — see data/accessions/01_hic_test_data.txt). Not
### login-node-appropriate; run as its own job between slurm/01_setup.slurm
### and slurm/03_pipeline.slurm.
# --------------------------------------------------------------
#SBATCH --job-name=hic_download_test_data
#SBATCH --output=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline/logs/slurm/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline/logs/slurm/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16gb
#SBATCH --time=6:00:00

echo "========== Test data download started at: $(date) =========="

SCRIPT_BASE=/home/u11/maarowosegbe/Hi-C_juicer_pipeline
DATA_BASE=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline

export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
export MAMBA_EXE=/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba
MC="${MAMBA_EXE}"

cd "${SCRIPT_BASE}"

DATA_DIR="${DATA_BASE}/data/fastq" "${MC}" run -n hic_juicer_env bash scripts/03_download_test_data.sh

echo "========== Test data download finished at: $(date) =========="
echo "Next: fill in config/samples.tsv's hic1 row with:"
echo "  ${DATA_BASE}/data/fastq/SRR1658570_1.fastq.gz"
echo "  ${DATA_BASE}/data/fastq/SRR1658570_2.fastq.gz"
echo "Then submit slurm/03_pipeline.slurm."
