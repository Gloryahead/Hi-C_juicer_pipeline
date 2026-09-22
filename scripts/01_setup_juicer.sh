#!/usr/bin/env bash
# 01_setup_juicer.sh
# Installs Juicer (aidenlab/juicer) and juicer_tools.jar into TOOLS_DIR.
#
# This is a host-side convenience for --use-conda runs (`use_conda: true`) or
# for running scripts/02_generate_restriction_sites.sh. It is NOT needed for
# --use-singularity runs: containers/hic_juicer.def bakes in the exact same
# pinned commit + jar version at build time, so the container never depends
# on GitHub/S3 being reachable from a compute node. The two install paths are
# kept version-identical on purpose (see docs/KEY_CONCEPTS.md) — if you bump
# one, bump the other.
#
# Usage:
#   TOOLS_DIR=/xdisk/haining/maarowosegbe/Hi-C_juicer_pipeline/tools bash scripts/01_setup_juicer.sh

set -euo pipefail

TOOLS_DIR="${TOOLS_DIR:-tools}"
JUICER_COMMIT="177eb610397e4207fc56db1df169b2d08d06d43a"   # aidenlab/juicer main, 2025-08-31 — see containers/hic_juicer.def
JUICER_TOOLS_VERSION="2.17.00"

mkdir -p "${TOOLS_DIR}"
cd "${TOOLS_DIR}"

if [[ -d juicer/.git ]]; then
  echo "  → juicer/ already cloned, checking out pinned commit"
  (cd juicer && git fetch --depth=1 origin "${JUICER_COMMIT}" 2>/dev/null || git fetch --unshallow; git checkout "${JUICER_COMMIT}")
else
  echo "  ↓ cloning aidenlab/juicer"
  git clone https://github.com/aidenlab/juicer.git
  (cd juicer && git checkout "${JUICER_COMMIT}")
fi

cd juicer
mkdir -p scripts/common
cp CPU/*.* scripts/common/
cp CPU/common/* scripts/common/
chmod +x scripts/common/*.sh

if [[ -f scripts/common/juicer_tools.jar ]]; then
  echo "  → juicer_tools.jar already present, skipping"
else
  echo "  ↓ juicer_tools_${JUICER_TOOLS_VERSION}.jar"
  wget -q "https://github.com/aidenlab/Juicebox/releases/download/v${JUICER_TOOLS_VERSION}/juicer_tools_${JUICER_TOOLS_VERSION}.jar"
  mv "juicer_tools_${JUICER_TOOLS_VERSION}.jar" scripts/common/juicer_tools.jar
fi

echo ""
echo "Juicer installed at ${TOOLS_DIR}/juicer (commit ${JUICER_COMMIT})."
echo "Point config/config.yaml's juicer.juicer_dir at ${TOOLS_DIR}/juicer if you used a custom TOOLS_DIR."
echo "Next: scripts/02_generate_restriction_sites.sh"
