process REHEADER_BAM {
    tag "$meta_id"
    label 'process_low'

    input:
    tuple val(meta_id), path(bam)

    output:
    tuple val(meta_id), path("${meta_id}.fixed.bam"), emit: bam
    tuple val(meta_id), path("${meta_id}.fixed.bam.bai"), emit: bai

    script:
    """
    # Check if the header already contains 'chr'
    if samtools view -H ${bam} | grep '^@SQ' | head -n 1 | grep -q 'SN:chr'; then
        echo "Headers already standardized for ${meta_id}. Linking original."
        # Use cp -rs (symbolic link) or copy to ensure the output filename matches the output block
        cp ${bam} ${meta_id}.fixed.bam
        samtools index ${meta_id}.fixed.bam
    else
        echo "Fixing headers for ${meta_id}..."
        samtools reheader <( \
            samtools view -H ${bam} | \
            sed 's/SN:\\([0-9XY]\\)/SN:chr\\1/g' | \
            sed 's/SN:MT/SN:chrM/g' \
        ) ${bam} > ${meta_id}.fixed.bam
        samtools index ${meta_id}.fixed.bam
    fi
    """
}
