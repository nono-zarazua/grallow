nextflow.enable.dsl=2

params.input_csv = "trial_samples.csv"
params.outdir    = "results"

include { REHEADER_BAM }                 from './modules/reheader'
include { MOSDEPTH; FLAGSTAT; NANOPLOT; MULTIQC } from './modules/qc'

workflow {
    ch_bams = Channel.fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.bam_path)) }

    REHEADER_BAM(ch_bams)

    // Run paraller QC
    MOSDEPTH(REHEADER_BAM.out.bam_bai)
    FLAGSTAT(REHEADER_BAM.out.bam_bai)
    NANOPLOT(REHEADER_BAM.out.bam_bai)

    // Extract just the files from the tuples (it[1] is the path)
    ch_mosdepth_files = MOSDEPTH.out.results.map { it[1] }
    ch_flagstat_files = FLAGSTAT.out.results.map { it[1] }
    ch_nanoplot_files = NANOPLOT.out.results.map { it[1] }

    // Mix them all together, collect into a list, and run MultiQC
    ch_all_qc = ch_mosdepth_files
        .mix(ch_flagstat_files, ch_nanoplot_files)
        .collect()

    MULTIQC(ch_all_qc)

    // We collect the mosdepth files into a list to pass to the evaluator
    ch_mosdepth_for_eval = ch_mosdepth_files.collect()

    // Pass the MultiQC data AND the Mosdepth files into our new process
    EVALUATE_QC(MULTIQC.out.data, ch_mosdepth_for_eval)
}
