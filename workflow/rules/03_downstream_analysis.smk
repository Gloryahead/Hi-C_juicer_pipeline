"""
Rules: downstream 3D-genome analysis on inter_30.hic (tutorial steps 9-12) —
TAD calling (Arrowhead), contact-matrix extraction, A/B compartment analysis
(eigenvector + Pearson correlation), and chromatin loop calling (HiCCUPS).
Juicebox itself (tutorial step 13) is an interactive desktop/web viewer, not
a batch step — it is documented in docs/KEY_CONCEPTS.md and README.md, not
automated here.
Environment: hic_juicer_env (01_hic_juicer.yml)
"""

_TADS   = config["tads"]
_MATS   = config["matrices"]
_COMP   = config["compartments"]
_LOOPS  = config["loops"]
_JAR    = f"{config['juicer']['juicer_dir']}/scripts/common/juicer_tools.jar"

# ── TAD calling (Arrowhead) ─────────────────────────────────────────────
rule arrowhead_tads:
    input: hic = f"{OUTDIR}/{{sample}}/aligned/inter_30.hic"
    output:
        bedpe = f"{OUTDIR}/{{sample}}/aligned/tads_genome_wide/{_TADS['resolution']}_blocks.bedpe",
        bed   = f"{OUTDIR}/{{sample}}/aligned/tads_genome_wide/{_TADS['resolution']}_blocks.bed",
    log: f"{LOGDIR}/arrowhead/{{sample}}.log"
    threads: 4
    resources: mem_mb=32000, runtime=240, slurm_partition="standard"
    params:
        activate = mamba_activate("hic_juicer_env"),
        chroms   = ",".join(config["genome"]["chroms"]),
        min_res  = _TADS["min_resolution"],
        res      = _TADS["resolution"],
        norm     = _TADS["norm"],
        outdir   = f"{OUTDIR}/{{sample}}/aligned/tads_genome_wide",
    shell:
        """
        set -eo pipefail
        {params.activate}
        mkdir -p {params.outdir}
        java -Xmx28g -jar {_JAR} arrowhead \
            -c {params.chroms} \
            -m {params.min_res} \
            -r {params.res} \
            -k {params.norm} \
            {input.hic} {params.outdir} \
            2>{log}
        # bedpe -> BED (chrom, start, end of the first anchor; matches the
        # tutorial's own "converted to BED format" step for genome-browser use)
        cut -f1-3 {output.bedpe} > {output.bed}
        """


# ── Assay-specific QC: TAD count against the tutorial's own stated
# expected range for a human genome (2,000-5,000) ───────────────────────
rule tad_qc_summary:
    input: bedpe = f"{OUTDIR}/{{sample}}/aligned/tads_genome_wide/{_TADS['resolution']}_blocks.bedpe"
    output: summary = f"{OUTDIR}/qc/{{sample}}_tad_qc_summary.txt"
    threads: 1
    resources: mem_mb=1000, runtime=10
    params:
        tad_min = config["qc_benchmarks"]["tad_count_min"],
        tad_max = config["qc_benchmarks"]["tad_count_max"],
    shell:
        """
        n=$(tail -n +2 {input.bedpe} | wc -l)
        status="PASS"
        if [ "$n" -lt {params.tad_min} ] || [ "$n" -gt {params.tad_max} ]; then
          status="WARN (outside {params.tad_min}-{params.tad_max} expected range)"
        fi
        {{
          echo "TAD count: $n"
          echo "STATUS: $status"
        }} > {output.summary}
        """


# ── Contact-matrix extraction (dump) ────────────────────────────────────
rule dump_contacts:
    input: hic = f"{OUTDIR}/{{sample}}/aligned/inter_30.hic"
    output: matrix = f"{OUTDIR}/{{sample}}/aligned/matrices/{{chrom}}_contacts_{_MATS['dump_resolution']}.txt"
    log: f"{LOGDIR}/dump_contacts/{{sample}}_{{chrom}}.log"
    threads: 2
    resources: mem_mb=16000, runtime=60, slurm_partition="standard"
    params:
        activate = mamba_activate("hic_juicer_env"),
        res      = _MATS["dump_resolution"],
        outdir   = f"{OUTDIR}/{{sample}}/aligned/matrices",
    shell:
        """
        set -eo pipefail
        {params.activate}
        mkdir -p {params.outdir}
        java -Xmx14g -jar {_JAR} dump \
            observed KR {input.hic} \
            {wildcards.chrom} {wildcards.chrom} \
            BP {params.res} \
            {output.matrix} \
            2>{log}
        """


# ── A/B compartment analysis: eigenvector decomposition + Pearson matrix ─
rule eigenvector_compartments:
    input: hic = f"{OUTDIR}/{{sample}}/aligned/inter_30.hic"
    output: f"{OUTDIR}/{{sample}}/aligned/compartments/{{chrom}}_eigenvector_{_COMP['resolution']}.txt"
    log: f"{LOGDIR}/eigenvector/{{sample}}_{{chrom}}.log"
    threads: 2
    resources: mem_mb=16000, runtime=60, slurm_partition="standard"
    params:
        activate = mamba_activate("hic_juicer_env"),
        res      = _COMP["resolution"],
        norm     = _COMP["norm"],
        outdir   = f"{OUTDIR}/{{sample}}/aligned/compartments",
    shell:
        """
        set -eo pipefail
        {params.activate}
        mkdir -p {params.outdir}
        java -Xmx14g -jar {_JAR} eigenvector \
            {params.norm} {input.hic} {wildcards.chrom} BP {params.res} \
            {output} \
            2>{log}
        """


rule pearsons_matrix:
    input: hic = f"{OUTDIR}/{{sample}}/aligned/inter_30.hic"
    output: f"{OUTDIR}/{{sample}}/aligned/compartments/{{chrom}}_pearsons_{_COMP['resolution']}.txt"
    log: f"{LOGDIR}/pearsons/{{sample}}_{{chrom}}.log"
    threads: 2
    resources: mem_mb=16000, runtime=60, slurm_partition="standard"
    params:
        activate = mamba_activate("hic_juicer_env"),
        res      = _COMP["resolution"],
        norm     = _COMP["norm"],
        outdir   = f"{OUTDIR}/{{sample}}/aligned/compartments",
    shell:
        """
        set -eo pipefail
        {params.activate}
        mkdir -p {params.outdir}
        java -Xmx14g -jar {_JAR} pearsons \
            {params.norm} {input.hic} {wildcards.chrom} BP {params.res} \
            {output} \
            2>{log}
        """


# ── Chromatin loop calling (HiCCUPS) — optional, see config.yaml ────────
# HiCCUPS is GPU-accelerated by design; juicer_tools' CPU fallback works but
# can take many hours-to-days on a genome-wide human map at 5-10kb
# resolution. loops.use_gpu flips the SLURM partition/resources below —
# point it at an actual GPU partition name for your HPC before enabling.
rule hiccups_loops:
    input: hic = f"{OUTDIR}/{{sample}}/aligned/inter_30.hic"
    output:
        bedpe = f"{OUTDIR}/{{sample}}/aligned/loops/postprocessed_pixels_{_LOOPS['resolutions'][0]}.bedpe",
    log: f"{LOGDIR}/hiccups/{{sample}}.log"
    threads: 8
    resources:
        mem_mb = 65536,
        runtime = 2880 if not _LOOPS["use_gpu"] else 360,
        slurm_partition = "gpu" if _LOOPS["use_gpu"] else "standard",
        slurm_extra = "--gres=gpu:1" if _LOOPS["use_gpu"] else "",
    params:
        activate = mamba_activate("hic_juicer_env"),
        m        = 512,
        res      = ",".join(str(r) for r in _LOOPS["resolutions"]),
        fdr      = ",".join(str(x) for x in _LOOPS["fdr"]),
        peak     = ",".join(str(x) for x in _LOOPS["peak_width"]),
        window   = ",".join(str(x) for x in _LOOPS["window"]),
        radius   = ",".join(str(x) for x in _LOOPS["cluster_radius"]),
        outdir   = f"{OUTDIR}/{{sample}}/aligned/loops",
    shell:
        """
        set -eo pipefail
        {params.activate}
        mkdir -p {params.outdir}
        java -Xmx60g -jar {_JAR} hiccups \
            -m {params.m} \
            -r {params.res} \
            -f {params.fdr} \
            -p {params.peak} \
            -i {params.window} \
            -d {params.radius} \
            {input.hic} {params.outdir} \
            2>{log}
        """
