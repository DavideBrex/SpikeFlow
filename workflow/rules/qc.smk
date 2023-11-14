rule fastqc: 
    input: 
        reads=get_reads
    output: 
        html = "{}results/QC/fastqc/{{id}}_fastqc.html".format(outdir),
    log:
        "{}results/logs/QC/fastqc/{{id}}_fastqc.log".format(outdir),
    conda: 
        "../envs/qc.yaml"
    params: 
        threads = 4,
        outdir = "{}results/QC/fastqc".format(outdir),
    shell:
        """
        fastqc -o {params.outdir} -t {params.threads} --extract {input} 2> {log}
        NAME="$(basename {input} .fastq.gz )"
        NAME+="_fastqc"
        mv {params.outdir}/$NAME/fastqc_report.html {output.html}
        rm -r {params.outdir}/$NAME
        """

rule multiqc: 
    input: 
        expand("{}results/QC/fastqc/{{id}}_fastqc.html".format(outdir), id=set(idSamples)),
    output:
        multiqc = "{}results/QC/multiqc/multiqc_report.html".format(outdir),
    log:
        "{}results/logs/QC/multiqc/multiqc.log".format(outdir),
    conda: 
        "../envs/qc.yaml"
    params:
        outdir = "{}results".format(outdir),
    shell:
        """ 
        multiqc {params.outdir} \
        --outdir {params.outdir}/QC/multiqc \
        --ignore-samples *spike \
        2> {log}
        """


rule plotFingerprint:
    input: 
        "{}results/bam/{{id}}.clean.bam".format(outdir),
    output: 
        qualMetrics = "{}results/QC/fingerPrint/{{id}}.qualityMetrics.tsv".format(outdir),
        raw_counts = "{}results/QC/fingerPrint/{{id}}.rawcounts.tsv".format(outdir),
        plot = "{}results/QC/fingerPrint/{{id}}.plot.pdf".format(outdir),
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
        spp = "{}results/QC/phantompeakqual/{{id}}.spp.out".format(outdir), 
    log:
        "{}results/logs/QC/phantompeakqual/{{id}}.spp.log".format(outdir),
    conda: 
        "../envs/Renv.yaml"
    threads:
        5
    params:
        out_dir = "{}results/QC/phantompeakqual/".format(outdir), 
    shell:
        """
        Rscript --vanilla workflow/scripts/run_spp_nodups.R \
        -c={input} -savp -rf -p={threads} -odir={params.out_dir}  -out={output} -tmpdir={params.out_dir}  \
        > {log} 2>&1
        """
