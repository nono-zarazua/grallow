process RUN_QC {
    tag "$meta_id"
    label 'process_medium'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    // Takes the output directly from REHEADER_BAM.out.bam
    tuple val(meta_id), path(bam), path(bai)

    output:
    tuple val(meta_id), path("mosdepth/*"), emit: mosdepth_results
    tuple val(meta_id), path("flagstat/*.flagstat"), emit: flagstat_results
    tuple val(meta_id), path("nanoplot/*"), emit: nanoplot_results

    script:
    """
    # Create output directories to keep results organized
    mkdir -p mosdepth flagstat nanoplot

    # 1. Mosdepth (Coverage stats)
    mosdepth -n -x "mosdepth/${meta_id}" ${bam}
    
    # 2. Samtools Flagstat (Alignment/mapping rates)
    samtools flagstat ${bam} > "flagstat/${meta_id}.flagstat"

    # 3. NanoPlot (Read length and quality distributions)
    NanoPlot --bam ${bam} \\
             --outdir nanoplot \\
             --maxlength 4000 \\
             --no_static \\
             -p "${meta_id}_"
    """
}
