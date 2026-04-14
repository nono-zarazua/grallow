nextflow.enable.dsl=2

params.input_csv = "trial_samples.csv"
params.outdir    = "results"

include { REHEADER_BAM }             from './modules/reheader'
include { MOSDEPTH; FLAGSTAT; NANOPLOT } from './modules/qc'

workflow {
    ch_bams = Channel.fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.bam_path)) }

    REHEADER_BAM(ch_bams)

    // Feed the output of REHEADER into the three QC tools
    // Because we use the same output channel, they run in parallel
    MOSDEPTH(REHEADER_BAM.out.bam)
    FLAGSTAT(REHEADER_BAM.out.bam)
    NANOPLOT(REHEADER_BAM.out.bam)
}
