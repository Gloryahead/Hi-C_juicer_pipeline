#!/bin/bash
# --------------------------------------------------------------
### PART 1: Request resources
# --------------------------------------------------------------
#SBATCH --job-name=hic_juicer_pipeline
#SBATCH --output=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline/logs/slurm/%x-%j.out
#SBATCH --error=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline/logs/slurm/%x-%j.err
#SBATCH --account=haining
#SBATCH --partition=standard

### The SLURM job itself only runs the Snakemake controller (low CPU/RAM).
### Snakemake submits individual rules as their own child SLURM jobs, each
### with the resources declared in workflow/rules/*.smk.
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=240:00:00       # controller must stay alive for the whole run

# --------------------------------------------------------------
### PART 2: Run Snakemake pipeline
# --------------------------------------------------------------

echo "========== Pipeline started at: $(date) =========="
echo "  Node    : $(hostname)"
echo "  Job ID  : ${SLURM_JOB_ID}"

SCRIPT_BASE=/home/u11/maarowosegbe/Hi-C_juicer_pipeline
DATA_BASE=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline

export MAMBA_ROOT_PREFIX=/groups/haining/maarowosegbe/micromamba
export MAMBA_EXE=/opt/ohpc/pub/apps/micromamba/2.0.2-2/bin/micromamba
MC="${MAMBA_EXE}"

cd "${SCRIPT_BASE}"

echo ""
echo "--- Snakemake dry run (DAG check) ---"
"${MC}" run -n snakemake_env \
  snakemake \
    --snakefile "${SCRIPT_BASE}/workflow/Snakefile" \
    --configfile "${SCRIPT_BASE}/config/config.yaml" \
    --profile "${SCRIPT_BASE}/workflow/profiles/slurm" \
    --dry-run --quiet

echo ""
echo "--- Full pipeline run ---"
"${MC}" run -n snakemake_env \
  snakemake \
    --snakefile "${SCRIPT_BASE}/workflow/Snakefile" \
    --configfile "${SCRIPT_BASE}/config/config.yaml" \
    --profile "${SCRIPT_BASE}/workflow/profiles/slurm" \
    --rerun-incomplete --keep-going --latency-wait 60 --jobs 100 --printshellcmds \
  2>&1 | tee "${DATA_BASE}/logs/slurm/snakemake_${SLURM_JOB_ID}.log"

echo "========== Pipeline finished at: $(date) =========="
