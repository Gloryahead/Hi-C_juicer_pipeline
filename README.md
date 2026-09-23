# NGS101 Hi-C (Juicer) Pipeline

A fully reproducible, modular pipeline covering the [NGS101 Hi-C
tutorial](https://ngs101.com/how-to-analyze-hi-c-data-for-absolute-beginners-from-raw-reads-to-3d-genome-organization-with-juicer/) —
raw paired-end reads through TADs, chromatin loops, and A/B compartments.
Runs on an HPC with SLURM via Snakemake (a Nextflow skeleton is also
included — see below), every tool version-pinned and containerized with
Apptainer (Singularity) for bit-for-bit reproducibility.

**Start here if you're new to this repo:**
[docs/KEY_CONCEPTS.md](docs/KEY_CONCEPTS.md) — the non-obvious decisions
(why Juicer over HiC-Pro/nf-core, why the core rule is monolithic, version
pinning, the HiCCUPS GPU caveat, and what the QC benchmarks actually mean)
in about 15 minutes.

---

## What this pipeline does

| Stage | Tutorial step | Tools |
|-------|---------------|-------|
| QC + trimming | Steps 5-7 | FastQC · Trim Galore |
| Core Juicer pipeline | Step 8 | BWA · Juicer (`juicer.sh`) — alignment, chimeric-read splitting, PCR-dup removal, `.hic` generation |
| Validation + QC gate | Step 8 (+ this repo) | `juicer_tools validate` · `workflow/scripts/hic_qc_check.py` (valid-pairs / inter-chromosomal-% sanity check) |
| TAD calling | Step 9 | `juicer_tools arrowhead` (KR-normalized, 25kb) + TAD-count QC gate |
| Contact matrices | Step 10 | `juicer_tools dump` (observed, KR, 25kb) |
| A/B compartments | Step 11 | `juicer_tools eigenvector` + `pearsons` (KR, 100kb) |
| Chromatin loops | Step 12 | `juicer_tools hiccups` (optional — GPU-shaped workload, see `config.yaml`) |
| Visualization | Step 13 | Juicebox (interactive desktop/web app — not automated; see below) |

The tutorial names exactly one worked example, **SRR1658570** — verified via
ENA/SRA to be GEO `GSM1551550` ("HIC001", Rao et al. 2014 *Cell*, GM12878 in
situ Hi-C, MboI). It's a full ~202M-read-pair replicate, not a toy dataset —
see [docs/KEY_CONCEPTS.md](docs/KEY_CONCEPTS.md) §9 before assuming a quick
run. `config/samples.tsv` ships with one row (`hic1`) wired to it; add more
rows for your own samples or biological replicates.

---

## Repository layout

```
.
├── config/
│   ├── config.yaml            # ALL pipeline parameters — edit this first
│   └── samples.tsv             # sample metadata (sample, sra_accession, restriction_enzyme, fastq_r1/r2)
├── containers/
│   ├── hic_juicer.def          # Apptainer: fastqc/trim-galore/bwa/samtools + Juicer + juicer_tools.jar (pinned)
│   └── build_all_sifs.sh
├── environments/                # numbered micromamba YAMLs
├── scripts/                     # setup/download scripts (00-03, numbered) + run_pipeline.sh
├── slurm/                       # SLURM submission scripts (setup, then pipeline)
├── data/accessions/             # accession provenance (SRR1658570 -> GEO/paper, verified)
├── runs/manifest.tsv            # run tracking (fill in yourself after a run)
├── docs/KEY_CONCEPTS.md         # the 15% that gives you 85% of the understanding
└── workflow/
    ├── Snakefile                 # master DAG (primary, maintained implementation)
    ├── profiles/slurm/           # SLURM executor config
    ├── rules/                    # 3 Snakemake rule files
    ├── scripts/hic_qc_check.py   # Hi-C QC-gate parser
    └── nextflow/                 # Nextflow skeleton (QC/trim + Juicer core only — see its main.nf header)
```

---

## Step-by-step setup

### Step 0 — Clone and configure

```bash
git clone <this-repo-url> Hi-C_juicer_pipeline
cd Hi-C_juicer_pipeline
```

Open `config/config.yaml` and check:
- `genome.restriction_enzyme` / `genome.restriction_sites` — must match your
  actual Hi-C library prep (see [docs/KEY_CONCEPTS.md](docs/KEY_CONCEPTS.md) §4)
- `loops.run` / `loops.use_gpu` — HiCCUPS is slow on CPU; see §7
- `juicer.threads`, and the `juicer_pipeline` rule's resources in
  `workflow/rules/02_juicer_core.smk` if your data is much bigger/smaller
  than the tutorial's ~200M-read-pair example

### Step 1 — Install micromamba (if not already on HPC)

```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

### Step 2 — Build environments

```bash
bash environments/setup_all_envs.sh
```

Creates `snakemake_env` and `hic_juicer_env`.

### Step 3 — Download references, install Juicer, generate restriction sites

```bash
REF_DIR=ref bash scripts/00_download_refs.sh                          # hg38 FASTA + BWA index (~1 GB, ~1h incl. index)
TOOLS_DIR=tools bash scripts/01_setup_juicer.sh                       # clones Juicer @ pinned commit + juicer_tools.jar
REF_DIR=ref TOOLS_DIR=tools ENZYME=MboI \
  bash scripts/02_generate_restriction_sites.sh                       # hg38_MboI.txt
```

Or submit `slurm/01_setup.slurm` to do all of the above (plus env creation)
in one SLURM job.

### Step 4 — Get the test data

```bash
DATA_DIR=data/fastq bash scripts/03_download_test_data.sh
```

Or submit `slurm/02_download_test_data.slurm` (recommended — ~202M read
pairs / ~40GB, not login-node-appropriate; see
[docs/KEY_CONCEPTS.md](docs/KEY_CONCEPTS.md) §9).

Downloads SRR1658570 as paired FASTQs. Then edit `config/samples.tsv`'s
`hic1` row: `fastq_r1`/`fastq_r2` → these two paths. Add more rows for your
own samples.

### Step 5 — Run the pipeline

```bash
# Dry run first — shows what will run without executing
bash scripts/run_pipeline.sh --dry-run

# Full run
bash scripts/run_pipeline.sh
```

Or submit `slurm/03_pipeline.slurm` (recommended — a long-running controller
job that submits every rule as its own child SLURM job automatically).

Run just through the core `.hic` file (skip TADs/loops/compartments):

```bash
bash scripts/run_pipeline.sh --target results/hic1/aligned/inter_30.hic
```

### Step 6 — (Optional) Build the Apptainer container

```bash
bash containers/build_all_sifs.sh
```

Then flip `use_singularity: true` in `config.yaml` and re-run. The container
bakes in the same pinned Juicer commit + `juicer_tools.jar` version as
`scripts/01_setup_juicer.sh` (see
[docs/KEY_CONCEPTS.md](docs/KEY_CONCEPTS.md) §6), so results are identical
either way.

### Step 7 — Visualize in Juicebox (manual — not part of the pipeline)

Download [Juicebox](https://github.com/aidenlab/Juicebox/releases) or use
the [web version](https://aidenlab.org/juicebox/). Load, per sample:

| File | Role |
|------|------|
| `results/{sample}/aligned/inter_30.hic` | primary contact map |
| `results/{sample}/aligned/tads_genome_wide/25000_blocks.bedpe` | TAD annotations |
| `results/{sample}/aligned/loops/postprocessed_pixels_5000.bedpe` | loop annotations (if `loops.run: true`) |
| `results/{sample}/aligned/compartments/{chrom}_eigenvector_100000.txt` | compartment track |

Start at 1MB resolution → zoom to 100kb (compartments) → 25kb (TAD/loop detail).

---

## Known limitations and manual steps

| Item | Issue | Workaround |
|------|-------|------------|
| Only one sample ships pre-wired | Tutorial names a single worked example (SRR1658570) | Add rows to `config/samples.tsv`; biological replicates are needed for IDR-style loop reproducibility filtering, which this tutorial doesn't demonstrate |
| HiCCUPS loop calling is GPU-shaped | `juicer_tools hiccups`'s CPU fallback is very slow genome-wide | `loops.use_gpu: true` + point `workflow/profiles/slurm` / the rule's `slurm_partition` at a real GPU partition; or `loops.run: false` to skip |
| Juicer has no recent version tags | `aidenlab/juicer` last tagged Nov 2021; real fixes are on `main` | Pinned to an exact, dated commit SHA in both the container and `scripts/01_setup_juicer.sh` — see `docs/KEY_CONCEPTS.md` §6 |
| Restriction-site script differs from the tutorial's own | Tutorial hand-rolls a cut-site finder | This repo uses Juicer's own `misc/generate_site_positions.py` instead — the implementation the rest of Juicer is validated against |
| QC benchmarks are heuristics | `config.yaml`'s `qc_benchmarks` are the tutorial author's stated ranges, not a citation | QC rules report `WARN`, never fail the DAG — treat as "go look," not a hard gate |
| Nextflow skeleton is partial | Only QC/trim + Juicer core are ported | Snakemake is the complete, maintained implementation; extend the Nextflow modules the same way if you adopt it as primary |

---

## Configuration quick reference

All parameters live in `config/config.yaml`. Key toggles:

```yaml
genome:
  restriction_enzyme: "MboI"   # must match the actual wet-lab library prep
loops:
  run: true                    # HiCCUPS loop calling
  use_gpu: false                # flip once a GPU partition is configured
use_singularity: false          # flip after containers/build_all_sifs.sh
```

---

## Reproducibility

- Every environment pins exact tool versions (`environments/*.yml`)
- Juicer itself is pinned to an exact commit SHA (no recent upstream tags —
  see `docs/KEY_CONCEPTS.md` §6), identically in the Apptainer `.def` and
  the host-side setup script
- `config/config.yaml` records every parameter used in a given run
- Snakemake logs every shell command with timestamps in `logs/`
- `--use-singularity` mode uses an immutable SIF container for bit-for-bit
  reproducibility
- `snakemake --report` generates an HTML provenance report after each run
- `runs/manifest.tsv` — append a row per run you want a durable record of
  (config snapshot, git commit, notes)
