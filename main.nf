nextflow.enable.dsl=2

// Parameters
params.input_bams = "data/raw_bams/*.bam"
params.outdir     = "results"

// Import the module
include { REHEADER_BAM } from './modules/reheader'

workflow {
    // 1. Create a channel of [SampleID, BamPath]
    // This uses the filename (minus .bam) as the ID
    ch_bams = Channel.fromPath(params.input_bams)
        .map { file -> tuple(file.simpleName, file) }

    // 2. Run the reheader process
    REHEADER_BAM(ch_bams)
    
    // 3. The output of REHEADER_BAM.out.bam is now ready for QUILT2
}
