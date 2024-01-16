rule fastqc:
    input:
        reads=get_reads,
    output:
        html="{}results/QC/fastqc/{{id}}_fastqc.html".format(outdir),
        zip="{}results/QC/fastqc/{{id}}_fastqc.zip".format(outdir),
    log:
        "{}results/logs/QC/fastqc/{{id}}_fastqc.log".format(outdir),
    conda:
        "../envs/qc.yaml"
    params:
        threads=config["threads"]["qc"],
        out_dir=lambda w, output: os.path.dirname(output.html),
    shell:
        """
        zcat {input.reads[0]} | fastqc stdin:{wildcards.id} -o {params.out_dir} &> {log}
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
    threads: config["threads"]["qc"]
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
    threads: config["threads"]["qc"]
    params:
        out_dir=lambda w, output: os.path.dirname(output.spp),
    script:
        "../scripts/run_spp_nodups.R"


rule create_qc_table_spike:
    input:
        logFile=expand(
            "{}results/logs/spike/{{id}}.removeSpikeDups".format(outdir),
            id=set(idSamples),
        ),
    output:
        tab="{}results/QC/Spikein_qc_FULL.tsv".format(outdir),
        tab_multiqc="{}results/QC/Spike-in_Reads_mqc.tsv".format(outdir),
    log:
        "{}results/logs/QC/Spike_in_qc.log".format(outdir),
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/createSpikeQCtab.py"


rule create_qc_table_epic2:
    input:
        logFile=expand(
            "{}results/logs/peakCalling/epic2/{{sample}}_summary.txt".format(outdir),
            sample=broadSamples,
        ),
    output:
        tab="{}results/QC/epic2_peaks_mqc.tsv".format(outdir),
    log:
        "{}results/logs/QC/epic2_qc.log".format(outdir),
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/createEpic2QCtab.py"


rule create_qc_table_edd:
    input:
        logFile=expand(
            "{}results/logs/peakCalling/edd/{{sample}}_summary.txt".format(outdir),
            sample=veryBroadSamples,
        ),
    output:
        tab="{}results/QC/edd_peaks_mqc.tsv".format(outdir),
    log:
        "{}results/logs/QC/edd_qc.log".format(outdir),
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/createEddQCtab.py"


rule create_qc_table_macs2:
    input:
        logFile=expand(
            "{}results/peakCalling/macs2_ref/{{sample}}_peaks.narrowPeak".format(
                outdir
            ),
            sample=narrowSamples,
        ),
    output:
        tab="{}results/QC/macs2_peaks_mqc.tsv".format(outdir),
    log:
        "{}results/logs/QC/macs2_qc.log".format(outdir),
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/createMacs2QCtab.py"


rule create_qc_table_peakAnnot:
    input:
        logFile=expand(
            "{}results/logs/peakCalling/peakAnnot/{{sample}}_annotInfo.txt".format(
                outdir
            ),
            sample=narrowSamples + broadSamples,
        ),
    output:
        tab="{}results/QC/peaks_annotation_mqc.tsv".format(outdir),
    log:
        "{}results/logs/QC/peaks_annotation_qc.log".format(outdir),
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/createPeakAnnotQCtab.py"


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
        tab_spike="{}results/QC/Spike-in_Reads_mqc.tsv".format(outdir),
        tab_macs2=rules.create_qc_table_macs2.output.tab,
        tab_epic2=rules.create_qc_table_epic2.output.tab,
        tab_edd=rules.create_qc_table_edd.output.tab,
        tab_peakAnnot=rules.create_qc_table_peakAnnot.output.tab,
    output:
        multiqc="{}results/QC/multiqc/multiqc_report.html".format(outdir),
    log:
        "{}results/logs/QC/multiqc/multiqc.log".format(outdir),
    conda:
        "../envs/qc.yaml"
    params:
        out_dir=lambda w, output: os.path.dirname(output.multiqc),
        whereTofind=lambda w, input: Path(input.tab_spike).parents[1],
    shell:
        """ 
        multiqc {params.whereTofind} \
        --outdir {params.out_dir} \
        --exclude macs2 \
        --ignore-samples *spike \
        --force \
        2> {log}
        """


ruleorder: epic2_callBroadPeaks > multiqc
