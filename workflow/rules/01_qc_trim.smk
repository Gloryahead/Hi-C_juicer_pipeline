"""
Rules: raw-read QC, adapter trimming, and Juicer-convention FASTQ staging.
Tutorial steps 5-7 ("Download & Quality Control Raw Data" through
"Prepare Input for Juicer").
Environment: hic_juicer_env (01_hic_juicer.yml)
"""

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

# ── FastQC on raw reads ───────────────────────────────────────────────
rule fastqc_raw:
    input:
        r1 = get_r1,
        r2 = get_r2,
    output:
        zip_r1 = f"{OUTDIR}/qc/fastqc/{{sample}}_R1_fastqc.zip",
        zip_r2 = f"{OUTDIR}/qc/fastqc/{{sample}}_R2_fastqc.zip",
    log:   f"{LOGDIR}/fastqc/{{sample}}.log"
    threads: 2
    resources:
        mem_mb = 4000, runtime = 30, slurm_partition = "standard",
    params:
        activate = mamba_activate("hic_juicer_env"),
        outdir   = f"{OUTDIR}/qc/fastqc",
    shell:
        """
        set -eo pipefail
        {params.activate}
        mkdir -p {params.outdir}
        fastqc --threads {threads} --outdir {params.outdir} {input.r1} {input.r2} 2>{log}
        base_r1=$(basename {input.r1} .fastq.gz); base_r2=$(basename {input.r2} .fastq.gz)
        mv {params.outdir}/${{base_r1}}_fastqc.zip {output.zip_r1}
        mv {params.outdir}/${{base_r2}}_fastqc.zip {output.zip_r2}
        """


# ── Trim Galore (adapter trim + embedded FastQC) ─────────────────────
# Tutorial values (q20, length>=20) are deliberately more permissive than the
# sibling WGS pipeline's (length>=50): Hi-C ligation-junction reads are
# expected to be chimeric, so filtering short reads discards real signal
# rather than noise.
rule trim_galore:
    input:
        r1 = get_r1,
        r2 = get_r2,
    output:
        r1     = temp(f"{OUTDIR}/trimmed/{{sample}}_R1_val_1.fq.gz"),
        r2     = temp(f"{OUTDIR}/trimmed/{{sample}}_R2_val_2.fq.gz"),
        report = f"{OUTDIR}/qc/trimming/{{sample}}_trimming_report.txt",
    log:   f"{LOGDIR}/trim_galore/{{sample}}.log"
    threads: config["trimming"]["cores"]
    resources:
        mem_mb = 8000, runtime = 120, slurm_partition = "standard",
    params:
        activate   = mamba_activate("hic_juicer_env"),
        quality    = config["trimming"]["quality"],
        min_length = config["trimming"]["min_length"],
        outdir     = f"{OUTDIR}/trimmed",
    shell:
        """
        set -eo pipefail
        {params.activate}
        trim_galore \
            --quality {params.quality} \
            --length {params.min_length} \
            --stringency 3 \
            --fastqc --paired --cores {threads} \
            --output_dir {params.outdir} \
            {input.r1} {input.r2} \
            2>{log}
        mv {params.outdir}/{wildcards.sample}*_val_1.fq.gz {output.r1} 2>>{log} || true
        mv {params.outdir}/{wildcards.sample}*_val_2.fq.gz {output.r2} 2>>{log} || true
        mv {params.outdir}/{wildcards.sample}*_trimming_report.txt {output.report} 2>>{log} || true
        """


# ── Stage trimmed reads into Juicer's expected <topDir>/fastq/ layout ──
# juicer.sh scans <topDir>/fastq/ for *_R1*.fastq(.gz) / *_R2*.fastq(.gz)
# pairs. The tutorial gunzips before handing off to Juicer; juicer.sh does
# accept gzipped input directly in recent versions, but this repo follows
# the tutorial's own tested path rather than an unverified shortcut.
rule prep_juicer_fastq:
    input:
        r1 = f"{OUTDIR}/trimmed/{{sample}}_R1_val_1.fq.gz",
        r2 = f"{OUTDIR}/trimmed/{{sample}}_R2_val_2.fq.gz",
    output:
        r1 = temp(f"{OUTDIR}/{{sample}}/fastq/{{sample}}_R1.fastq"),
        r2 = temp(f"{OUTDIR}/{{sample}}/fastq/{{sample}}_R2.fastq"),
    log: f"{LOGDIR}/prep_juicer_fastq/{{sample}}.log"
    threads: 1
    resources: mem_mb=4000, runtime=60, slurm_partition="standard"
    shell:
        """
        set -eo pipefail
        mkdir -p {OUTDIR}/{wildcards.sample}/fastq
        gunzip -c {input.r1} > {output.r1} 2>{log}
        gunzip -c {input.r2} > {output.r2} 2>>{log}
        """


# ── MultiQC (aggregate raw + trimming QC; Juicer's own alignment stats are
# summarized separately by 03_downstream_analysis.smk's hic_qc_summary rule —
# MultiQC has no Juicer module) ────────────────────────────────────────
rule multiqc:
    input:
        expand(f"{OUTDIR}/qc/fastqc/{{sample}}_R1_fastqc.zip", sample=ALL_SAMPLES),
        expand(f"{OUTDIR}/qc/trimming/{{sample}}_trimming_report.txt", sample=ALL_SAMPLES),
    output: f"{OUTDIR}/qc/multiqc_report.html"
    log:    f"{LOGDIR}/multiqc/multiqc.log"
    resources: mem_mb=8000, runtime=30
    params:
        activate = mamba_activate("hic_juicer_env"),
        indir    = f"{OUTDIR}/qc",
        outdir   = f"{OUTDIR}/qc",
    shell:
        """
        set -eo pipefail
        {params.activate}
        multiqc {params.indir} --outdir {params.outdir} --filename multiqc_report.html --force 2>{log}
        """
