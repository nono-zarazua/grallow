process REHEADER_BAM {
    tag "$meta_id"
    label 'process_medium'
    // Moves results from the 'work' directory to your project folder 
    publishDir "${params.outdir}/${params.barch}/reheadered", mode: 'copy'

    input:
    tuple val(meta_id), path(bam)

    output:
    tuple val(meta_id), path("${meta_id}.fixed.bam"), path("${meta_id}.fixed.bam.bai"), emit: bam_bai

    script:
    """
    # Check if header already contains 'chr' [cite: 2]
    if samtools view -H ${bam} | grep '^@SQ' | head -n 1 | grep -q 'SN:chr'; then
        echo "Headers already standardized for ${meta_id}. Linking original."
        # Create a copy with the expected output name [cite: 2]
        cp ${bam} ${meta_id}.fixed.bam
        samtools index ${meta_id}.fixed.bam
    else
        echo "Fixing headers for ${meta_id}..."
        # Reheader logic using your sed patterns [cite: 3]
        samtools reheader <( \\
            samtools view -H ${bam} | \\
            sed 's/SN:\\([0-9XY]\\)/SN:chr\\1/g' | \\
            sed 's/SN:MT/SN:chrM/g' \\
        ) ${bam} > ${meta_id}.fixed.bam
        samtools index ${meta_id}.fixed.bam
    fi
    """
}
