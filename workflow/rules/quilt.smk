
rule quilt_prepare_regular:
    input:
        vcf=rules.subset_refpanel_by_chunkid.output.vcf,
    output:
        os.path.join(
            OUTDIR_QUILT1_REF,
            "refsize{size}",
            "QUILT_prepared_reference.{chrom}.chunk_{chunkid}.RData",
        ),
    params:
        time=config["time"],
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        nGen=config["quilt1"]["nGen"],
        buffer=config["quilt1"]["buffer"],
        lowram=config["quilt1"]["lowram"],
        impute_rare_common=config["quilt1"]["impute_rare_common"],
        rare_af_threshold=config["quilt1"]["rare_af_threshold"],
        start=get_quilt_chunk_region_start,
        end=get_quilt_chunk_region_end,
        gmap=if_use_quilt_map_in_refpanel,
    log:
        os.path.join(
            OUTDIR_QUILT1_REF,
            "refsize{size}",
            "QUILT_prepared_reference.{chrom}.chunk_{chunkid}.RData.llog",
        ),
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        (
        if [ -s {params.gmap} ];then \
        {params.time} -v QUILT_prepare_reference.R \
            --genetic_map_file='{params.gmap}' \
            --reference_vcf_file={input.vcf} \
            --chr={wildcards.chrom} \
            --regionStart={params.start} \
            --regionEnd={params.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --use_hapMatcherR={params.lowram} \
            --use_mspbwt=FALSE \
            --impute_rare_common={params.impute_rare_common} \
            --rare_af_threshold={params.rare_af_threshold} \
            --outputdir={params.outdir} \
            --output_file={output} \
        ; else \
        {params.time} -v QUILT_prepare_reference.R \
            --reference_vcf_file={input.vcf} \
            --chr={wildcards.chrom} \
            --regionStart={params.start} \
            --regionEnd={params.end} \
            --buffer={params.buffer} \
            --use_hapMatcherR={params.lowram} \
            --nGen={params.nGen} \
            --use_mspbwt=FALSE \
            --impute_rare_common={params.impute_rare_common} \
            --rare_af_threshold={params.rare_af_threshold} \
            --outputdir={params.outdir} \
            --output_file={output} \
        ; fi
        ) &> {log}
        """


rule quilt_run_regular:
    input:
        vcf=rules.subset_refpanel_by_chunkid.output.vcf,
        bams=rules.bamlist.output,
        rdata=rules.quilt_prepare_regular.output,
    output:
        os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.regular.{chrom}.chunk_{chunkid}.vcf.gz",
        ),
    params:
        time=config["time"],
        nGen=config["quilt1"]["nGen"],
        buffer=config["quilt1"]["buffer"],
        start=get_quilt_chunk_region_start,
        end=get_quilt_chunk_region_end,
        Ksubset=config["quilt1"]["Ksubset"],
        Knew=config["quilt1"]["Knew"],
        nGibbsSamples=config["quilt1"]["nGibbsSamples"],
        lowram=config["quilt1"]["lowram"],
        impute_rare_common=config["quilt1"]["impute_rare_common"],
        rare_af_threshold=config["quilt1"]["rare_af_threshold"],
        n_seek_its=config["quilt1"]["n_seek_its"],
        block_gibbs=config["quilt1"]["small_ref_panel_block_gibbs_iterations"],
        gibbs_iters=config["quilt1"]["small_ref_panel_gibbs_iterations"],
    log:
        os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.regular.{chrom}.chunk_{chunkid}.vcf.gz.llog",
        ),
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        {params.time} -v QUILT.R \
            --reference_vcf_file={input.vcf} \
            --prepared_reference_filename={input.rdata} \
            --bamlist={input.bams} \
            --chr={wildcards.chrom} \
            --regionStart={params.start} \
            --regionEnd={params.end} \
            --buffer={params.buffer} \
            --nGen={params.nGen} \
            --use_mspbwt=FALSE \
            --Ksubset={params.Ksubset} \
            --Knew={params.Knew} \
            --nGibbsSamples={params.nGibbsSamples} \
            --use_hapMatcherR={params.lowram} \
            --impute_rare_common={params.impute_rare_common} \
            --rare_af_threshold={params.rare_af_threshold} \
            --n_seek_its={params.n_seek_its} \
            --small_ref_panel_block_gibbs_iterations='{params.block_gibbs}' \
            --small_ref_panel_gibbs_iterations={params.gibbs_iters} \
            --output_filename={output} &> {log}
        """


rule quilt_ligate_regular:
    input:
        get_quilt_regular_outputs,
    output:
        vcf=os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.regular.{chrom}.vcf.gz",
        ),
        sample=os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.regular.{chrom}.bcf.sample",
        ),
        tmp=temp(
            os.path.join(
                OUTDIR_QUILT1,
                "refsize{size}",
                "{chrom}",
                "quilt.down{depth}x.regular.{chrom}.bcf",
            )
        ),
        lst=temp(
            os.path.join(
                OUTDIR_QUILT1,
                "refsize{size}",
                "{chrom}",
                "quilt.down{depth}x.regular.{chrom}.vcf.list",
            )
        ),
    log:
        os.path.join(
            OUTDIR_QUILT1,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.regular.{chrom}.vcf.gz.llog",
        ),
    params:
        N="quilt_ligate_regular",
        extra=config["extra_buffer_in_quilt"],
        sample=config["samples"],
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        ( \
        if [ {params.extra} -gt 0 ];then \
           echo {input} | tr ' ' '\n' > {output.lst} && \
           bcftools concat --ligate --file-list {output.lst} --threads 4 -o {output.tmp} \
        ; else \
           echo {input} | tr ' ' '\n' > {output.lst} && \
           bcftools concat --file-list {output.lst} --threads 4 -o {output.tmp} \
        ; fi 
        awk 'NR>1 {{ print $1 }}' {params.sample} > {output.sample} && \
        bcftools reheader -s {output.sample} -o {output.vcf} {output.tmp} && \
        bcftools index -f {output.vcf}
        ) &> {log}
        """


rule quilt_prepare_mspbwt:
    input:
        vcf=rules.subset_refpanel_by_chunkid.output.vcf,
    output:
        os.path.join(
            OUTDIR_QUILT2_REF,
            "refsize{size}",
            "QUILT_prepared_reference.{chrom}.chunk_{chunkid}.RData",
        ),
    params:
        time=config["time"],
        N="quilt_prepare_mspbwt",
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        nGen=config["quilt2"]["nGen"],
        buffer=config["quilt2"]["buffer"],
        start=get_quilt_chunk_region_start,
        end=get_quilt_chunk_region_end,
        gmap=if_use_quilt_map_in_refpanel,
        lowram=config["quilt2"]["lowram"],
        impute_rare_common=config["quilt2"]["impute_rare_common"],
        rare_af_threshold=config["quilt2"]["rare_af_threshold"],
        nindices=config["quilt2"]["mspbwt-nindices"],
    log:
        os.path.join(
            OUTDIR_QUILT2_REF,
            "refsize{size}",
            "QUILT_prepared_reference.{chrom}.chunk_{chunkid}.RData.llog",
        ),
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        (
        # We define a helper to handle failure: if R fails, create an empty file
        run_quilt_prep() {{
            if ! "$@"; then
                echo "WARNING: QUILT_prepare_reference.R failed. Assuming empty region/centromere."
                echo "Creating empty placeholder: {output}"
                rm -f {output} && touch {output}
            fi
        }}

        if [ -s {params.gmap} ]; then
            run_quilt_prep {params.time} -v QUILT_prepare_reference.R \
                --genetic_map_file='{params.gmap}' \
                --reference_vcf_file={input.vcf} \
                --chr={wildcards.chrom} \
                --regionStart={params.start} \
                --regionEnd={params.end} \
                --use_hapMatcherR={params.lowram} \
                --buffer={params.buffer} \
                --nGen={params.nGen} \
                --use_mspbwt=TRUE \
                --impute_rare_common={params.impute_rare_common} \
                --rare_af_threshold={params.rare_af_threshold} \
                --mspbwt_nindices={params.nindices} \
                --outputdir={params.outdir} \
                --output_file={output}
        else
            run_quilt_prep {params.time} -v QUILT_prepare_reference.R \
                --reference_vcf_file={input.vcf} \
                --chr={wildcards.chrom} \
                --regionStart={params.start} \
                --regionEnd={params.end} \
                --buffer={params.buffer} \
                --use_hapMatcherR={params.lowram} \
                --nGen={params.nGen} \
                --use_mspbwt=TRUE \
                --rare_af_threshold={params.rare_af_threshold} \
                --impute_rare_common={params.impute_rare_common} \
                --mspbwt_nindices={params.nindices} \
                --outputdir={params.outdir} \
                --output_file={output}
        fi
        ) &> {log}
        """


rule quilt_run_mspbwt:
    input:
        vcf=rules.subset_refpanel_by_chunkid.output.vcf,
        bams=rules.bamlist.output,
        rdata=rules.quilt_prepare_mspbwt.output,
        sex="config/samples-sex.tsv"
    output:
        os.path.join(
            OUTDIR_QUILT2,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.mspbwt.{chrom}.chunk_{chunkid}.vcf.gz",
        ),
    params:
        raw_samples="config/samples.tsv",
        time=config["time"],
        N="quilt_run_mspbwt",
        nGen=config["quilt2"]["nGen"],
        buffer=config["quilt2"]["buffer"],
        start=get_quilt_chunk_region_start,
        end=get_quilt_chunk_region_end,
        Ksubset=config["quilt2"]["Ksubset"],
        Knew=config["quilt2"]["Knew"],
        nGibbsSamples=config["quilt2"]["nGibbsSamples"],
        n_seek_its=config["quilt2"]["n_seek_its"],
        lowram=config["quilt2"]["lowram"],
        rare_af_threshold=config["quilt2"]["rare_af_threshold"],
        impute_rare_common=config["quilt2"]["impute_rare_common"],
        block_gibbs=config["quilt2"]["small_ref_panel_block_gibbs_iterations"],
        gibbs_iters=config["quilt2"]["small_ref_panel_gibbs_iterations"],
        mspbwtM=config["quilt2"]["mspbwtM"],
        mspbwtL=config["quilt2"]["mspbwtL"],
    log:
        os.path.join(
            OUTDIR_QUILT2,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.mspbwt.{chrom}.chunk_{chunkid}.vcf.gz.llog",
        ),
    conda:
        "../envs/quilt.yaml"
    threads: 1
    shell:
        """
        # 1. GENERATE CLEAN SAMPLE LIST
        # Extract 1st column, skip header (NR>1), save to temp file
        awk 'NR>1 {{print $1}}' {params.raw_samples} > {output}.sample_names

        # Check if the RData file exists and has size > 0
        if [ -s {input.rdata} ]; then
            {params.time} -v QUILT.R \
                --reference_vcf_file={input.vcf} \
                --prepared_reference_filename={input.rdata} \
                --bamlist={input.bams} \
                --use_hapMatcherR={params.lowram} \
                --impute_rare_common={params.impute_rare_common} \
                --chr={wildcards.chrom} \
                --regionStart={params.start} \
                --regionEnd={params.end} \
                --buffer={params.buffer} \
                --nGen={params.nGen} \
                --use_mspbwt=TRUE \
                --mspbwtM={params.mspbwtM} \
                --mspbwtL={params.mspbwtL} \
                --Ksubset={params.Ksubset} \
                --Knew={params.Knew} \
                --nGibbsSamples={params.nGibbsSamples} \
                --n_seek_its={params.n_seek_its} \
                --rare_af_threshold={params.rare_af_threshold} \
                --small_ref_panel_block_gibbs_iterations='{params.block_gibbs}'\
                --small_ref_panel_gibbs_iterations={params.gibbs_iters} \
                --output_filename={output} &> {log}
        else
            echo "Skipping QUILT.R because input RData is empty (likely centromere/gap)." &> {log}
            # Create a valid empty VCF (header only) so ligation doesn't fail
            # We use the input reference chunk to grab a valid header
            bcftools view -h {input.vcf} | bgzip -c > {output} 2>> {log}
        fi
        """        



rule quilt_ligate_mspbwt:
    input:
        get_quilt_mspbwt_outputs,
    output:
        vcf=os.path.join(
            OUTDIR_QUILT2,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.mspbwt.{chrom}.vcf.gz",
        ),
        sample=temp(
            os.path.join(
                OUTDIR_QUILT2,
                "refsize{size}",
                "{chrom}",
                "quilt.down{depth}x.mspbwt.{chrom}.vcf.sample",
            )
        ),
        tmp=temp(
            os.path.join(
                OUTDIR_QUILT2,
                "refsize{size}",
                "{chrom}",
                "quilt.down{depth}x.mspbwt.{chrom}.tmp.vcf.gz",
            )
        ),
        lst=temp(
            os.path.join(
                OUTDIR_QUILT2,
                "refsize{size}",
                "{chrom}",
                "quilt.down{depth}x.mspbwt.{chrom}.vcf.list",
            )
        ),
    log:
        os.path.join(
            OUTDIR_QUILT2,
            "refsize{size}",
            "{chrom}",
            "quilt.down{depth}x.mspbwt.{chrom}.vcf.gz.llog",
        ),
    params:
        N="quilt_ligate_mspbwt",
        extra=config["extra_buffer_in_quilt"],
        sample=config["samples"],
    conda:
        "../envs/quilt.yaml"
    shell:
       """
        set -euo pipefail
        (
        echo "Generating verified input list for {wildcards.chrom}..."

        > {output.lst}
        
        for f in {input}; do
            # NEGATIVE FILTER (Robust Version):
            # We search the raw file for 'HG00096' (found in Ref Panel headers).
            # -m 1: Stop reading after the first match (fast).
            # If zgrep finds it (exit code 0), it's a Ghost Chunk -> SKIP.
            
            if zgrep -m 1 -q "HG00096" "$f"; then
                echo "WARNING: Skipping $f (Found HG00096 - treating as Ghost Chunk/Gap)."
            else
                echo "$f" >> {output.lst}
            fi
        done

        # Safety Check: Did we find any valid chunks?
        if [ ! -s {output.lst} ]; then
            echo "ERROR: No valid chunks found for {wildcards.chrom}!"
            exit 1
        fi

        # --- 2. CONCATENATION ---
        if [ {params.extra} -gt 0 ]; then
            echo "Attempting strict ligation..."
            if ! bcftools concat \
                --ligate \
                --file-list {output.lst} \
                --threads 4 \
                -o {output.tmp} \
                -O z; then

                echo ">> Ligation failed. Switching to ROBUST MODE (--allow-overlaps)..."
                rm -f {output.tmp}
                bcftools concat \
                    --file-list {output.lst} \
                    --allow-overlaps \
                    --threads 4 \
                    -o {output.tmp} \
                    -O z
            fi
        else
            bcftools concat \
                --file-list {output.lst} \
                --threads 4 \
                -o {output.tmp} \
                -O z
        fi

        # --- 3. REHEADER & INDEX ---
        awk 'NR>1 {{ print $1 }}' {params.sample} > {output.sample}
        if [ ! -s {output.sample} ]; then
             echo "Sample1" > {output.sample}
        fi

        bcftools reheader \
            -s {output.sample} \
            -o {output.vcf} \
            {output.tmp}

        bcftools index -t -f {output.vcf}


        echo "Job finished successfully."
        ) &>> {log}
        """

rule quilt_concat_genome:
    input:
        vcfs=get_quilt_genome_inputs
    output:
        vcf=os.path.join(OUTDIR_QUILT2, "refsize{size}", "quilt.down{depth}x.mspbwt.genome.vcf.gz"),
        tbi=os.path.join(OUTDIR_QUILT2, "refsize{size}", "quilt.down{depth}x.mspbwt.genome.vcf.gz.tbi")
    log:
        os.path.join(OUTDIR_QUILT2, "refsize{size}", "quilt.down{depth}x.mspbwt.genome.vcf.gz.log")
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        echo "Sorting chromosome inputs to ensure numeric order..." > {log}
        
        # Sort the input files naturally (chr2 before chr10)
        SORTED_VCFS=$(echo "{input.vcfs}" | tr ' ' '\\n' | sort -V | tr '\\n' ' ')
        
        echo "Concatenating..." >> {log}
        bcftools concat \
            --threads 8 \
            -O z \
            -o {output.vcf} \
            $SORTED_VCFS \
            2>> {log}

        echo "Indexing..." >> {log}
        bcftools index -t {output.vcf} 2>> {log}
        """


rule quilt_split_by_sample:
    input:
        vcf=os.path.join(OUTDIR_QUILT2, "refsize{size}", f"quilt.down{config['downsample'][0]}x.mspbwt.genome.vcf.gz"),
        tbi=os.path.join(OUTDIR_QUILT2, "refsize{size}", f"quilt.down{config['downsample'][0]}x.mspbwt.genome.vcf.gz.tbi")
    output:
        vcf=os.path.join(OUTDIR_QUILT2, "refsize{size}", "split_files", "{sample}.vcf.gz"),
        tbi=os.path.join(OUTDIR_QUILT2, "refsize{size}", "split_files", "{sample}.vcf.gz.tbi")
    log:
        os.path.join(OUTDIR_QUILT2, "refsize{size}", "split_files", "{sample}.split.log")
    conda:
        "../envs/quilt.yaml"
    shell:
        """
        echo "Extracting sample {wildcards.sample} from combined genome..." > {log}
        
        # Extract only the target sample, keep all INFO tags, and compress
        bcftools view \
            -s {wildcards.sample} \
            --threads 4 \
            -O z \
            -o {output.vcf} \
            {input.vcf} 2>> {log}

        echo "Indexing the sample VCF..." >> {log}
        bcftools index -t {output.vcf} 2>> {log}
        
        echo "Done!" >> {log}
        """
