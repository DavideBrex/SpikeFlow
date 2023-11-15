rule fastqc: 
    input: 
        reads=get_reads
    output: 
        html="{}results/QC/fastqc/{{id}}_fastqc.html".format(outdir),
    log:
        "{}results/logs/QC/fastqc/{{id}}_fastqc.log".format(outdir),
    conda: 
        "../envs/qc.yaml"
    params: 
        threads=4,
        outdir="{}results/QC/fastqc".format(outdir),
    shell:
        """
        fastqc -o {params.outdir} -t {params.threads} --extract {input} 2> {log}
        NAME="$(basename {input} .fastq.gz )"
        NAME+="_fastqc"
        mv {params.outdir}/$NAME/fastqc_report.html {output.html}
        rm -r {params.outdir}/$NAME
        """


rule plotFingerprint:
    input: 
        "{}results/bam/{{id}}.clean.bam".format(outdir),
    output: 
        qualMetrics="{}results/QC/fingerPrint/{{id}}.qualityMetrics.tsv".format(outdir),
        raw_counts="{}results/QC/fingerPrint/{{id}}.rawcounts.tsv".format(outdir),
        plot="{}results/QC/fingerPrint/{{id}}.plot.pdf".format(outdir),
    log:
        "{}results/logs/QC/fingerPrint/{{id}}.fingerPrint.log".format(outdir),
    conda: 
        "../envs/qc.yaml"
    threads:
        5,
    shell:
        """
        plotFingerprint -b {input} \
        -p {threads} \
        --outQualityMetrics {output.qualMetrics} \
        --outRawCounts {output.raw_counts} \
        --plotFile {output.plot} \
        2> {log}
        """


rule phantom_peak_qual:
    input: 
        "{}results/bam/{{id}}.clean.bam".format(outdir),
    output:
        spp="{}results/QC/phantompeakqual/{{id}}.spp.out".format(outdir), 
    log:
        "{}results/logs/QC/phantompeakqual/{{id}}.spp.log".format(outdir),
    conda: 
        "../envs/Renv.yaml"
    threads:
        5
    params:
        out_dir="{}results/QC/phantompeakqual/".format(outdir), 
    shell:
        """
        Rscript --vanilla workflow/scripts/run_spp_nodups.R \
        -c={input} -savp -rf -p={threads} -odir={params.out_dir}  -out={output} -tmpdir={params.out_dir}  \
        > {log} 2>&1
        """


rule create_qc_table_spike:
    input:
        logFile=expand("{}results/logs/spike/{{id}}.removeSpikeDups".format(outdir), id=set(idSamples)),
        normFactor=expand("{}results/logs/spike/{{id}}.normFactor".format(outdir), id=set(idSamples)),
    output:
        tab="{}results/QC/Spikein_qc_FULL.tsv".format(outdir),
        tab_multiqc="{}results/QC/Spike-in_Reads_mqc.tsv".format(outdir),
    log:
        "{}results/logs/QC/Spike_in_qc.log".format(outdir),
    run:
        import pandas as pd
        dict_allsamp = {}
        for f in input.logFile:
            idName = f.split("/")[-1].split(".")[0]
            with open(f, "r") as file:
                info_sample = file.read().strip().split("\n")
                sample_nr = int(info_sample[1].split(":")[-1])
                spike_nr = int(info_sample[2].split(":")[-1])
                dict_allsamp[idName] =  [(sample_nr + spike_nr), sample_nr, spike_nr, ((spike_nr / (sample_nr + spike_nr)) * 100)]

        df_info = pd.DataFrame.from_dict(dict_allsamp, orient='index')
        df_info.columns=["Total mapped Reads", "Tot. Sample Reads","Tot. spikeIn Reads","Percentage spikeIn"]
        df_info.to_csv(output.tab, float_format='%.3f')
        df_info.iloc[:,1:3].to_csv(output.tab_multiqc, float_format='%.2f')



rule multiqc: 
    input: 
        fastqc=expand("{}results/QC/fastqc/{{id}}_fastqc.html".format(outdir), id=set(idSamples)),
        spp=expand("{}results/QC/phantompeakqual/{{id}}.spp.out".format(outdir), id=set(idSamples)), 
        plot=expand("{}results/QC/fingerPrint/{{id}}.plot.pdf".format(outdir), id=set(idSamples)),
        tab="{}results/QC/Spike-in_Reads_mqc.tsv".format(outdir),
    output:
        multiqc="{}results/QC/multiqc/multiqc_report.html".format(outdir),
    log:
        "{}results/logs/QC/multiqc/multiqc.log".format(outdir),
    conda: 
        "../envs/qc.yaml"
    params:
        outdir="{}results".format(outdir),
    shell:
        """ 
        multiqc {params.outdir} \
        --outdir {params.outdir}/QC/multiqc \
        --ignore-samples *spike \
        --force \
        2> {log}
        """


