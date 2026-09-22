# Key Concepts — Hi-C (Juicer) Pipeline

A ~15-minute read covering the non-obvious decisions baked into this repo.
Source tutorial: [ngs101.com — "How to Analyze Hi-C Data for Absolute
Beginners: From Raw Reads to 3D Genome Organization with
Juicer"](https://ngs101.com/how-to-analyze-hi-c-data-for-absolute-beginners-from-raw-reads-to-3d-genome-organization-with-juicer/).

---

## 1. Why Juicer, not HiC-Pro / nf-core/hic

`nf-core/hic` exists and is actively maintained, but it wraps **HiC-Pro +
cooler**: a different aligner/dedup/normalization chain with its own stats
format and its own TAD/compartment tooling. The tutorial this repo implements
teaches the **Aiden-lab Juicer** toolchain specifically — `juicer.sh` for
alignment through `.hic` generation, and `juicer_tools.jar` for everything
downstream (Arrowhead TADs, HiCCUPS loops, eigenvector/Pearson compartments).
These aren't interchangeable: outputs, normalization methods, and QC metrics
differ. So this repo is a **custom Snakemake/Nextflow assembly around
containerized Juicer**, not an nf-core pipeline with different tools swapped
in — that would silently change what the pipeline actually computes.

## 2. Why `juicer.sh` is one monolithic rule, not decomposed

`juicer.sh` manages its own working-directory layout
(`topDir/{fastq,splits,aligned}/...`) across alignment, chimeric-read
splitting, sorting, PCR-duplicate removal, and `.hic` matrix building. It
isn't designed to resume from arbitrary intermediate files the way GATK's
BQSR/HaplotypeCaller stages are — so Snakemake wraps it as one rule
(`workflow/rules/02_juicer_core.smk::juicer_pipeline`), the same reasoning
that keeps HiC-Pro as a single containerized unit in pipelines that use it.
You still get DAG-level resume: if `inter_30.hic` already exists, the whole
rule is skipped.

## 3. `inter_30.hic` vs `inter.hic`

`juicer.sh` produces both. `inter_30.hic` keeps only alignments with
**MAPQ ≥ 30** (confidently unique); `inter.hic` includes everything
(MAPQ ≥ 0), which is noisier but higher-depth. Every downstream rule in this
pipeline reads `inter_30.hic` — this matches the tutorial and is the
standard choice for TAD/loop calling, where mapping ambiguity directly
produces false contacts.

## 4. Restriction enzyme is a fact about the data, not a tunable parameter

`config.genome.restriction_enzyme` (default `MboI`) must match whatever
enzyme was actually used in the wet-lab Hi-C library prep — it is not
something to sweep for "better resolution." Using the wrong enzyme here
silently corrupts fragment-level deduplication and normalization, it won't
error out. MboI/DpnII (4-cutters, `GATC`) give kb-scale resolution at higher
sequencing cost; HindIII (6-cutter, `AAGCTT`) is coarser but needs less
depth.

## 5. Restriction-site file: upstream `generate_site_positions.py`, not a
   custom reimplementation

The tutorial includes its own short Python script for finding cut sites.
This repo instead calls Juicer's own `misc/generate_site_positions.py`
(`scripts/02_generate_restriction_sites.sh`) — the field-tested implementation
the rest of Juicer's fragment-aware steps are validated against, rather than
a second, independently-written cut-site finder that could disagree with it
at edge cases (chromosome-end fragments, ambiguous IUPAC bases, etc.).

## 6. Version pinning: Juicer has no recent tags

`aidenlab/juicer` tags infrequently (last tag: `encode_tag_disk1`, Nov 2021;
real fixes have landed on `main` since). A commit SHA is therefore the only
reproducible reference available — this repo pins
`177eb610397e4207fc56db1df169b2d08d06d43a` (verified current as of
2025-08-31) in **both** `containers/hic_juicer.def` (baked into the SIF at
build time) and `scripts/01_setup_juicer.sh` (the `--use-conda` host-side
path), deliberately kept identical. `juicer_tools_2.17.00.jar` is a real
tagged release asset (`aidenlab/Juicebox` v2.17.00) and is pinned the same
way in both places.

## 7. HiCCUPS needs a GPU to be fast

`juicer_tools hiccups` is GPU-accelerated by design; the CPU fallback works
but can take many hours to days on a genome-wide human map at 5-10kb
resolution. `config.yaml`'s `loops.use_gpu` toggles the SLURM partition and
resource request in `workflow/rules/03_downstream_analysis.smk::hiccups_loops`
— set it to `true` and point it at your HPC's actual GPU partition name
before relying on loop calls for anything beyond a quick sanity check.
`loops.run: false` skips loop calling entirely if no GPU is available.

## 8. QC benchmarks are the tutorial author's heuristics, not a citation

`config.yaml`'s `qc_benchmarks` block (inter-chromosomal 5-40%, valid pairs
10-50M/100-500M, TAD count 2,000-5,000) comes from the tutorial's own stated
expected ranges for a human genome. They're reasonable sanity checks — real
signals of a failed ligation step or over/under-sequencing — but they are
not a peer-reviewed statistical threshold. `workflow/scripts/hic_qc_check.py`
and the `tad_qc_summary` rule report `WARN`, never fail the pipeline: treat a
WARN as "go look at this," not as a hard rejection.

## 9. `SRR1658570` is a full replicate, not a toy dataset

The tutorial's one worked example (verified against ENA/SRA:
`data/accessions/01_hic_test_data.txt`) is GEO `GSM1551550` ("HIC001") from
Rao et al. 2014 (*Cell*), GM12878 in situ Hi-C, MboI digestion — **~202M read
pairs**. That's within this pipeline's own "high-resolution" QC band. Budget
the resources `workflow/rules/02_juicer_core.smk` requests (128GB /
48h for `juicer_pipeline`) — this is not a quick few-minute demo run.

## 10. Juicebox is not part of the pipeline, by design

The tutorial's final step (loading `inter_30.hic` + TAD/loop/compartment
tracks into Juicebox) is an interactive desktop or web application
(`aidenlab.org/juicebox`) for visual exploration — there is nothing to batch
or containerize. README.md documents which output files to load and in what
order (1MB → zoom to 100kb for compartments → 25kb for loop detail); this
pipeline's job ends at producing those files.

## 11. Trimming thresholds differ from the sibling WGS/RNA-seq repos

`trimming.min_length: 20` here (vs. 50 in `WGS_variant_pipeline`). Hi-C reads
spanning a ligation junction are legitimately short/chimeric — aggressive
length filtering tuned for clean genomic reads would discard real proximity-
ligation signal, not just noise.
