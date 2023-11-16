rule fastqc:
    input:
        reads=get_reads,
    output:
        html="{}results/QC/fastqc/{{id}}_fastqc.html".format(outdir),
    log:
        "{}results/logs/QC/fastqc/{{id}}_fastqc.log".format(outdir),
    conda:
        "../envs/qc.yaml"
    params:
        threads=4,
        out_dir=lambda w, output: os.path.dirname(output.html),
    shell:
        """
        fastqc -o {params.out_dir} -t {params.threads} --extract {input} 2> {log}
        NAME="$(basename {input} .fastq.gz )"
        NAME+="_fastqc"
        mv {params.out_dir}/$NAME/fastqc_report.html {output.html}
        rm -r {params.out_dir}/$NAME
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
    threads: 5
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
    threads: 5
    params:
        out_dir=lambda w, output: os.path.dirname(output.spp),
    shell:
        """
        Rscript --vanilla workflow/scripts/run_spp_nodups.R \
        -c={input} -savp -rf -p={threads} -odir={params.out_dir}  -out={output} -tmpdir={params.out_dir}  \
        > {log} 2>&1
        """


rule create_qc_table_spike:
    input:
        logFile=expand(
            "{}results/logs/spike/{{id}}.removeSpikeDups".format(outdir),
            id=set(idSamples),
        ),
        normFactor=expand(
            "{}results/logs/spike/{{id}}.normFactor".format(outdir), id=set(idSamples)
        ),
    output:
        tab="{}results/QC/Spikein_qc_FULL.tsv".format(outdir),
        tab_multiqc="{}results/QC/Spike-in_Reads_mqc.tsv".format(outdir),
    log:
        "{}results/logs/QC/Spike_in_qc.log".format(outdir),
    script:
        "../scripts/createSpikeQCtab.py"


rule multiqc:
    input:
        fastqc=expand(
            "{}results/QC/fastqc/{{id}}_fastqc.html".format(outdir), id=set(idSamples)
        ),
        spp=expand(
            "{}results/QC/phantompeakqual/{{id}}.spp.out".format(outdir),
            id=set(idSamples),
        ),
        plot=expand(
            "{}results/QC/fingerPrint/{{id}}.plot.pdf".format(outdir),
            id=set(idSamples),
        ),
        tab="{}results/QC/Spike-in_Reads_mqc.tsv".format(outdir),
    output:
        multiqc="{}results/QC/multiqc/multiqc_report.html".format(outdir),
    log:
        "{}results/logs/QC/multiqc/multiqc.log".format(outdir),
    conda:
        "../envs/qc.yaml"
    params:
        out_dir=lambda w, output: os.path.dirname(output.multiqc),
        whereTofind=lambda w, input: Path(input.tab).parents[1],
    shell:
        """ 
        multiqc {params.whereTofind} \
        --outdir {params.out_dir} \
        --ignore-samples *spike \
        --force \
        2> {log}
        """
