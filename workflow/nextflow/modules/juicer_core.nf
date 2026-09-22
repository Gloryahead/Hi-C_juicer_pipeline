process JUICER_CORE {
    // One process wrapping juicer.sh end-to-end, same rationale as the
    // Snakemake rule of the same purpose (workflow/rules/02_juicer_core.smk):
    // juicer.sh's internal stages share a working-directory layout it
    // manages itself and isn't designed to resume from arbitrary
    // intermediate files, so it isn't decomposed further here.
    tag "${sample}"
    publishDir "${params.outdir}/${sample}/aligned", mode: 'copy'

    input:
        tuple val(sample), path(r1_fastq), path(r2_fastq)
    output:
        tuple val(sample), path("inter_30.hic"), path("inter_30.txt"), path("merged_nodups.txt"), emit: hic

    script:
    """
    mkdir -p fastq
    ln -s ../${r1_fastq} fastq/${sample}_R1.fastq
    ln -s ../${r2_fastq} fastq/${sample}_R2.fastq

    bash ${params.juicer_dir}/scripts/common/juicer.sh \
        -D ${params.juicer_dir} \
        -d \$(pwd) \
        -g ${params.genome_id} \
        -s ${params.restriction_enzyme} \
        -p ${params.chrom_sizes} \
        -y ${params.restriction_sites} \
        -z ${params.genome_fasta} \
        -t ${task.cpus}

    mv aligned/inter_30.hic aligned/inter_30.txt aligned/merged_nodups.txt .
    """
}
