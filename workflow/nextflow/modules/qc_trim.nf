process FASTQC_RAW {
    tag "${sample}"
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'

    input:
        tuple val(sample), path(r1), path(r2)
    output:
        tuple val(sample), path("*_fastqc.zip"), emit: zips

    script:
    """
    fastqc --threads ${task.cpus} ${r1} ${r2}
    """
}

process TRIM_GALORE {
    tag "${sample}"
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: '*trimming_report.txt'

    input:
        tuple val(sample), path(r1), path(r2)
    output:
        tuple val(sample), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed
        path "*trimming_report.txt"

    script:
    // Tutorial values (q20, length>=20) — deliberately permissive vs a typical
    // WGS trim (length>=50): Hi-C ligation-junction reads are chimeric by
    // design, so aggressive length filtering would discard real signal.
    """
    trim_galore \
        --quality ${params.trim_quality} \
        --length ${params.trim_min_length} \
        --stringency 3 \
        --fastqc --paired --cores ${task.cpus} \
        ${r1} ${r2}
    """
}
