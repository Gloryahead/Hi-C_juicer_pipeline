"""
Rule: the core Juicer pipeline (tutorial step 8) — two-step BWA alignment
(normal + chimeric-read splitting), PCR-duplicate removal, fragment-level
filtering, and .hic contact-matrix generation.

This is deliberately ONE Snakemake rule wrapping juicer.sh end-to-end rather
than decomposed into per-stage rules (align / dedup / merge / build-matrix).
juicer.sh's internal stages share a tightly-coupled working-directory layout
(topDir/{fastq,splits,aligned}/...) that Juicer itself manages and that isn't
designed to be resumed from arbitrary intermediate files — the same reasoning
that keeps HiC-Pro as a single containerized unit elsewhere in this pipeline
family (see docs/KEY_CONCEPTS.md). Snakemake still gets DAG-level resume: if
inter_30.hic already exists, this whole rule is skipped.

Environment: hic_juicer_env (01_hic_juicer.yml)
"""

_JUICER = config["juicer"]
_GENOME = config["genome"]

# ── Core Juicer pipeline: FASTQ -> inter_30.hic ─────────────────────────
rule juicer_pipeline:
    input:
        r1  = f"{OUTDIR}/{{sample}}/fastq/{{sample}}_R1.fastq",
        r2  = f"{OUTDIR}/{{sample}}/fastq/{{sample}}_R2.fastq",
        ref = _GENOME["fasta"],
        chrom_sizes = _GENOME["chrom_sizes"],
        restriction_sites = _GENOME["restriction_sites"],
    output:
        hic       = f"{OUTDIR}/{{sample}}/aligned/inter_30.hic",
        stats     = f"{OUTDIR}/{{sample}}/aligned/inter_30.txt",
        merged    = f"{OUTDIR}/{{sample}}/aligned/merged_nodups.txt",
    log: f"{LOGDIR}/juicer/{{sample}}.log"
    threads: _JUICER["threads"]
    resources:
        # Real-depth human Hi-C (100-500M read pairs) routinely needs
        # 64-128GB RAM for the sort/dedup stage and can run 12-48h. The
        # tutorial's own SRR1658570 (~202M read pairs, see
        # data/accessions/01_hic_test_data.txt) is squarely in that range —
        # this is not a quick toy run, size the job accordingly.
        mem_mb = 131072, runtime = 2880, slurm_partition = "standard",
    params:
        activate  = mamba_activate("hic_juicer_env"),
        juicer_sh = f"{_JUICER['juicer_dir']}/scripts/common/juicer.sh",
        juicer_d  = _JUICER["juicer_dir"],
        topdir    = f"{OUTDIR}/{{sample}}",
        genome_id = _GENOME["assembly"],
        enzyme    = get_enzyme,
    shell:
        """
        set -eo pipefail
        {params.activate}
        bash {params.juicer_sh} \
            -D {params.juicer_d} \
            -d {params.topdir} \
            -g {params.genome_id} \
            -s {params.enzyme} \
            -p {input.chrom_sizes} \
            -y {input.restriction_sites} \
            -z {input.ref} \
            -t {threads} \
            2>{log}
        """


# ── Validate the .hic file structure (tutorial's own explicit check) ────
rule juicer_validate:
    input: hic = f"{OUTDIR}/{{sample}}/aligned/inter_30.hic"
    output: log = f"{OUTDIR}/{{sample}}/aligned/inter_30_validate.log"
    log: f"{LOGDIR}/juicer_validate/{{sample}}.log"
    threads: 1
    resources: mem_mb=32000, runtime=60, slurm_partition="standard"
    params:
        activate = mamba_activate("hic_juicer_env"),
        jar      = f"{_JUICER['juicer_dir']}/scripts/common/juicer_tools.jar",
    shell:
        """
        set -eo pipefail
        {params.activate}
        java -Xmx28g -jar {params.jar} validate {input.hic} > {output.log} 2>{log}
        """


# ── Assay-specific QC gate: parse Juicer's own inter_30.txt statistics and
# flag values outside the tutorial's stated expected ranges. This is a real
# pipeline step producing a checkable file, not a manual post-hoc read of
# the stats — see docs/KEY_CONCEPTS.md for why these are heuristics, not a
# hard pass/fail bar. ────────────────────────────────────────────────────
rule hic_qc_summary:
    input:
        stats = f"{OUTDIR}/{{sample}}/aligned/inter_30.txt",
    output:
        summary = f"{OUTDIR}/qc/{{sample}}_hic_qc_summary.txt",
    log: f"{LOGDIR}/hic_qc_summary/{{sample}}.log"
    threads: 1
    resources: mem_mb=2000, runtime=15
    params:
        activate = mamba_activate("hic_juicer_env"),
        inter_min = config["qc_benchmarks"]["inter_chrom_pct_min"],
        inter_max = config["qc_benchmarks"]["inter_chrom_pct_max"],
        pairs_min = config["qc_benchmarks"]["valid_pairs_min"],
        pairs_hi  = config["qc_benchmarks"]["valid_pairs_high_res"],
    shell:
        """
        set -eo pipefail
        {params.activate}
        mkdir -p {OUTDIR}/qc
        python3 workflow/scripts/hic_qc_check.py \
            --stats {input.stats} \
            --sample {wildcards.sample} \
            --inter-min {params.inter_min} --inter-max {params.inter_max} \
            --pairs-min {params.pairs_min} --pairs-high-res {params.pairs_hi} \
            --out {output.summary} \
            2>{log}
        """
