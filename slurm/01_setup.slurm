#!/bin/bash
# --------------------------------------------------------------
### One-time setup: build environments, download refs, install Juicer
# --------------------------------------------------------------
#SBATCH --job-name=hic_pipeline_setup
#SBATCH --output=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline/logs/slurm/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline/logs/slurm/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32gb
#SBATCH --time=12:00:00

echo "========== Setup started at: $(date) =========="

SCRIPT_BASE=/home/u11/maarowosegbe/Hi-C_juicer_pipeline
DATA_BASE=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline

export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
export MAMBA_EXE=/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba
MC="${MAMBA_EXE}"

cd "${SCRIPT_BASE}"

echo "--- Building micromamba environments ---"
bash environments/setup_all_envs.sh

echo "--- Downloading hg38 reference + BWA index ---"
REF_DIR="${DATA_BASE}/ref" "${MC}" run -n hic_juicer_env bash scripts/00_download_refs.sh

echo "--- Installing Juicer (pinned commit) + juicer_tools.jar ---"
TOOLS_DIR="${DATA_BASE}/tools" bash scripts/01_setup_juicer.sh

echo "--- Generating MboI restriction-site file ---"
REF_DIR="${DATA_BASE}/ref" TOOLS_DIR="${DATA_BASE}/tools" ENZYME=MboI \
  "${MC}" run -n hic_juicer_env bash scripts/02_generate_restriction_sites.sh

echo "========== Setup finished at: $(date) =========="
echo "Next: download test data (scripts/03_download_test_data.sh), fill in"
echo "config/samples.tsv's fastq_r1/fastq_r2, then submit slurm/02_pipeline.slurm."
