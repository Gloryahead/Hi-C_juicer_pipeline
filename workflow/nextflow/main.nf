#!/usr/bin/env nextflow
/*
 * Hi-C (Juicer) pipeline — Nextflow skeleton.
 * Covers QC/trim -> Juicer core (FASTQ -> inter_30.hic) only. The Snakemake
 * implementation (workflow/Snakefile) is this project's complete, maintained
 * pipeline, including TADs/loops/compartments/QC gates — see that file and
 * docs/KEY_CONCEPTS.md. Run with: nextflow run main.nf -profile slurm
 */
nextflow.enable.dsl = 2

include { FASTQC_RAW; TRIM_GALORE } from './modules/qc_trim.nf'
include { JUICER_CORE }             from './modules/juicer_core.nf'

workflow {
    samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .filter { it.fastq_r1 && it.fastq_r2 }
        .map { row -> tuple(row.sample, file(row.fastq_r1), file(row.fastq_r2)) }

    FASTQC_RAW(samples_ch)
    trimmed_ch = TRIM_GALORE(samples_ch).trimmed
    JUICER_CORE(trimmed_ch)
}
