#!/usr/bin/env bash
# setup_all_envs.sh
# Creates both micromamba environments for the pipeline. Run once.
#
# Usage:
#   bash environments/setup_all_envs.sh
#
# Requires micromamba on PATH (or MAMBA_EXE set — see scripts/run_pipeline.sh
# for the exact HPC paths this repo assumes).

set -euo pipefail

ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MC="${MAMBA_EXE:-micromamba}"

for yml in "${ENV_DIR}"/*.yml; do
  name=$(grep '^name:' "${yml}" | awk '{print $2}')
  echo ""
  echo "═══════════════════════════════════════════════════"
  echo "  Creating: ${name}  (from $(basename "${yml}"))"
  echo "═══════════════════════════════════════════════════"
  if "${MC}" env list | grep -q "  ${name}$\| ${name} \*"; then
    echo "  → Already exists; skipping. Remove with '${MC} env remove -n ${name}' to rebuild."
    continue
  fi
  "${MC}" env create -f "${yml}" --yes
done

echo ""
echo "All environments created. Next: scripts/00_download_refs.sh, then"
echo "scripts/01_setup_juicer.sh to install Juicer itself into hic_juicer_env's reach."
