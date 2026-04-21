process MOSDEPTH {
    tag "$meta_id"
    label 'process_high'
    publishDir "${params.outdir}/${params.batch}/qc/mosdepth", mode: 'copy'

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
    publishDir "${params.outdir}/${params.batch}/qc/flagstat", mode: 'copy'

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
    label 'process_high'
    publishDir "${params.outdir}/${params.batch}/qc/nanoplot", mode: 'copy'

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
             -p "${meta_id}_" \\
             -t 4
    """
}

process MULTIQC {
    label 'process_low'
    // Save the final report to the main results folder
    publishDir "${params.outdir}/${params.batch}/multiqc", mode: 'copy'

    input:
    // Takes a collected list of all QC files from all samples
    path qc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data"       , emit: data

    script:
    """
    # MultiQC will scan the current Nextflow work directory 
    # where all the qc_files have been symlinked
    multiqc .
    """
}

process EVALUATE_QC {
    label 'process_low'
    
    debug true 
    
    publishDir "${params.outdir}/${params.batch}/evaluation", mode: 'copy'

    input:
    path multiqc_data_dir
    path mosdepth_files 

    output:
    path "qc_summary.csv", emit: report
    
    script:
    """
    # Run the Python evaluator, outputting to CSV
    evaluate_qc.py ${multiqc_data_dir}/multiqc_general_stats.txt . qc_summary.csv
    """
}
