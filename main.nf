nextflow.enable.dsl=2

// Point this to your trial CSV
params.input_csv = "trial_samples.csv"
params.outdir    = "results"

include { REHEADER_BAM } from './modules/reheader'

workflow {
    // 1. Parse the CSV and create the [ID, path] tuple
    ch_bams = Channel.fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.bam_path)) }

    // 2. Feed into your reheader process
    REHEADER_BAM(ch_bams)
}
