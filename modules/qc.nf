process MOSDEPTH {
    tag "$meta_id"
    label 'process_low'
    publishDir "${params.outdir}/qc/mosdepth", mode: 'copy'

    input:
    tuple val(meta_id), path(bam), path(bai)

    output:
    tuple val(meta_id), path("${meta_id}*"), emit: results

    script:
    """
    mosdepth -n -x "${meta_id}" ${bam}
    """
}

process FLAGSTAT {
    tag "$meta_id"
    label 'process_low'
    publishDir "${params.outdir}/qc/flagstat", mode: 'copy'

    input:
    tuple val(meta_id), path(bam), path(bai)

    output:
    tuple val(meta_id), path("${meta_id}.flagstat"), emit: results

    script:
    """
    samtools flagstat ${bam} > "${meta_id}.flagstat"
    """
}

process NANOPLOT {
    tag "$meta_id"
    label 'process_medium'
    publishDir "${params.outdir}/qc/nanoplot", mode: 'copy'

    input:
    tuple val(meta_id), path(bam), path(bai)

    output:
    tuple val(meta_id), path("${meta_id}_*"), emit: results

    script:
    """
    NanoPlot --bam ${bam} \\
             --outdir . \\
             --maxlength 4000 \\
             --no_static \\
             -p "${meta_id}_"
    """
}
