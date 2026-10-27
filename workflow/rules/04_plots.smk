"""
Rule: publication-style figures (contact maps, compartments, TAD/loop summaries)
drawn from the text outputs of 03_downstream_analysis.smk by
workflow/scripts/plot_hic.py. Draws only; recomputes nothing.
Environment: hic_juicer_env (01_hic_juicer.yml: needs matplotlib, numpy, pandas)
"""

_PLOTS = config["plots"]
_PCHR  = _PLOTS["chrom"]
_PNAMES = ["fig1_contact_maps", "fig2_compartments", "fig3_tads_loops"]

# Fine-resolution matrix for the zoom panel. juicer_tools 2.17.00's
# `dump ... chr1:x1:x2 chr1:x1:x2` region syntax silently writes an empty file
# for every normalization and bin size (checked on this map 2026-10-22), so the
# whole chromosome is dumped to a temporary file and the window is cut out.
rule dump_zoom:
    input: hic = f"{OUTDIR}/{{sample}}/aligned/inter_30.hic"
    output: f"{OUTDIR}/{{sample}}/aligned/matrices/{_PCHR}_zoom_{_PLOTS['zoom_start']}_{_PLOTS['zoom_end']}_{_PLOTS['zoom_res']}.txt"
    log: f"{LOGDIR}/plots/{{sample}}_dump_zoom.log"
    threads: 2
    resources: mem_mb=32000, runtime=120, slurm_partition="standard"
    params:
        activate = mamba_activate("hic_juicer_env"),
        jar      = f"{config['juicer']['juicer_dir']}/scripts/common/juicer_tools.jar",
        chrom    = _PCHR,
        res      = _PLOTS["zoom_res"],
        norm     = config["matrices"]["norm"],
        zstart   = _PLOTS["zoom_start"],
        zend     = _PLOTS["zoom_end"],
        tmp      = f"{OUTDIR}/{{sample}}/aligned/matrices/tmp_{_PCHR}_{_PLOTS['zoom_res']}.txt",
    shell:
        """
        set -eo pipefail
        {params.activate}
        java -Xmx28g -jar {params.jar} dump \
            observed {params.norm} {input.hic} \
            {params.chrom} {params.chrom} \
            BP {params.res} \
            {params.tmp} \
            > {log} 2>&1
        # juicer_tools can exit 0 after a normalization error and leave an
        # empty file (SCALE may be missing at 5 kb); fail the job instead.
        test -s {params.tmp}
        awk -v s={params.zstart} -v e={params.zend} \
            '$1>=s && $1<e && $2>=s && $2<e' {params.tmp} > {output}
        rm -f {params.tmp}
        test -s {output}
        """


rule plot_figures:
    input:
        # Snakemake does not track the contents of a script called from a
        # shell block; listing it here reruns the plots whenever it changes.
        script  = "workflow/scripts/plot_hic.py",
        zoom    = rules.dump_zoom.output,
        matrix  =f"{OUTDIR}/{{sample}}/aligned/matrices/{_PCHR}_contacts_{config['matrices']['dump_resolution']}.txt",
        pearson = f"{OUTDIR}/{{sample}}/aligned/compartments/{_PCHR}_pearsons_{config['compartments']['resolution']}.txt",
        eigen   = f"{OUTDIR}/{{sample}}/aligned/compartments/{_PCHR}_eigenvector_{config['compartments']['eigenvector_resolution']}.txt",
        tads    = f"{OUTDIR}/{{sample}}/aligned/tads_genome_wide/{config['tads']['resolution']}_blocks.bedpe",
        # merged_loops.bedpe combines every HiCCUPS resolution (5 kb + 10 kb);
        # hiccups_loops only declares the first resolution's file, but this one
        # is written by the same run.
        loops   = (f"{OUTDIR}/{{sample}}/aligned/loops/merged_loops.bedpe"
                   if config["loops"]["run"] else []),
    output:
        expand(f"{OUTDIR}/{{{{sample}}}}/plots/{{name}}.{{ext}}", name=_PNAMES, ext=["png", "pdf"]),
    log: f"{LOGDIR}/plots/{{sample}}.log"
    threads: 1
    # The chr1 matrix is held dense (about 10,000 x 10,000 float32 = 400 MB,
    # plus a log10 copy); 16 GB is generous.
    resources: mem_mb=16000, runtime=30, slurm_partition="standard"
    params:
        activate = mamba_activate("hic_juicer_env"),
        chrom    = _PCHR,
        sizes    = config["genome"]["chrom_sizes"],
        m_res    = config["matrices"]["dump_resolution"],
        z_res    = _PLOTS["zoom_res"],
        p_res    = config["compartments"]["resolution"],
        e_res    = config["compartments"]["eigenvector_resolution"],
        zstart   = _PLOTS["zoom_start"],
        zend     = _PLOTS["zoom_end"],
        loops    = lambda wc, input: f"--loops {input.loops}" if config["loops"]["run"] else "",
        outdir   = f"{OUTDIR}/{{sample}}/plots",
    shell:
        """
        set -eo pipefail
        {params.activate}
        python3 workflow/scripts/plot_hic.py \
            --sample {wildcards.sample} --chrom {params.chrom} \
            --chrom-sizes {params.sizes} \
            --matrix {input.matrix} --matrix-res {params.m_res} \
            --zoom-matrix {input.zoom} --zoom-res {params.z_res} \
            --pearsons {input.pearson} --pearsons-res {params.p_res} \
            --eigenvector {input.eigen} --eigenvector-res {params.e_res} \
            --tads {input.tads} {params.loops} \
            --zoom-start {params.zstart} --zoom-end {params.zend} \
            --outdir {params.outdir} \
            > {log} 2>&1
        """
