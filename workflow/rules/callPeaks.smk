rule macs2_callNarrowPeak:
    input:
        treatment="{}results/bam/{{sample}}.clean.bam".format(outdir),
        control=lambda w: "{}results/bam/{}.clean.bam".format(
            outdir, sample_to_input[w.sample]
        ),
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext(
            "{}results/peakCalling/macs2_ref/{{sample}}".format(outdir),
            "_peaks.xls",
            "_peaks.narrowPeak",
        ),
    log:
        "{}results/logs/peakCalling/macs2/{{sample}}_callpeak.log".format(outdir),
    params:
        lambda w: "--format BAM --nomodel"
        if is_single_end(w.sample)
        else "--format BAMPE"
        + " --gsize "
        + str(config["params"]["deeptools"]["effective_genome_length"])
        + " --pvalue "
        + config["params"]["peakCalling"]["macs2"]["pvalue"]
        + " --keep-dup all",
    benchmark:
        "{}results/.benchmarks/{{sample}}.macs2.benchmark.txt".format(outdir)
    wrapper:
        "v2.6.0/bio/macs2/callpeak"


rule epic2_callBroadPeaks:
    input:
        treatment="{}results/bam/{{sample}}.clean.bam".format(outdir),
        control=lambda w: "{}results/bam/{}.clean.bam".format(
            outdir, sample_to_input[w.sample]
        ),
    output:
        broadPeaks="{}results/peakCalling/epic2/{{sample}}_broadPeaks.bed".format(
            outdir
        ),
    log:
        "{}results/logs/peakCalling/epic2/{{sample}}_callpeak.log".format(outdir),
    params:
        fdr=config["params"]["peakCalling"]["epic2"]["fdr"],
        genome=config["params"]["peakCalling"]["genome"],
    conda:
        "../envs/various.yaml"
    shell:
        """
        epic2 -t {input.treatment} \
              -c {input.control} \
              --guess-bampe \
              -fdr {params.fdr} \
              --genome {params.genome} \
              --output {output.broadPeaks} \
              2> {log}
        """


rule edd_callVeryBroadPeaks:
    input:
        treatment="{}results/bam/{{sample}}.clean.bam".format(outdir),
        control=lambda w: "{}results/bam/{}.clean.bam".format(
            outdir, sample_to_input[w.sample]
        ),
    output:
        veryBroadPeaks="{}results/peakCalling/edd/{{sample}}/{{sample}}_peaks.bed".format(
            outdir
        ),
    log:
        "{}results/logs/peakCalling/edd/{{sample}}_callpeak.log".format(outdir),
    params:
        chrom_size=config["params"]["peakCalling"]["chrom_sizes"],
        fdr=config["params"]["peakCalling"]["edd"]["fdr"],
        blacklist=config["resources"]["ref"]["blacklist"],
        output_dir="{}results/peakCalling/edd/{{sample}}".format(outdir),
    conda:
        "../envs/edd.yaml"
    threads: 8
    shell:
        """
        edd --nprocs {threads} \
            --fdr {params.fdr} \
            {params.chrom_size} \
            {params.blacklist} \
            {input.treatment} {input.control} \
            {params.output_dir} \
            2> {log}
        mv {params.output_dir}/log.txt {params.output_dir}/{wildcards.sample}_runlog.txt
        mv {params.output_dir}/{wildcards.sample}_runlog.txt $(dirname {log})
        """
