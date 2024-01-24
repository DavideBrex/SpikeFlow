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
        + str(config["params"]["peakCalling"]["macs2"]["pvalue"])
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
        summary="{}results/logs/peakCalling/epic2/{{sample}}_summary.txt".format(outdir),
        bedfile="{}results/logs/peakCalling/epic2/{{sample}}.bed".format(outdir),
    log:
        "{}results/logs/peakCalling/epic2/{{sample}}_callpeak.log".format(outdir),
    params:
        fdr=config["params"]["peakCalling"]["epic2"]["fdr"],
        egf=config["params"]["peakCalling"]["epic2"]["egf"],
        chrom_size=config["params"]["peakCalling"]["chrom_sizes"],
    benchmark:
        "{}results/.benchmarks/{{sample}}.epic2.benchmark.txt".format(outdir)
    conda:
        "../envs/various.yaml"
    shell:
        """
        epic2 -t {input.treatment} \
              -c {input.control} \
              --guess-bampe \
              -fdr {params.fdr} \
              --chromsizes {params.chrom_size} \
              --effective-genome-fraction {params.egf} \
              --output {output.broadPeaks} \
              &> {log}
        
        # we create a simpler bed output file with only the chromosome, start, end and peak score
        awk '{{print $1, $2, $3, $6, int($5), $6, $7, $4, $9}}' {output.broadPeaks} > {output.bedfile}
        # Count the number of peaks and create a summary file
        echo "Total Peaks: $(wc -l < {output.broadPeaks})" > {output.summary}
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
        summary="{}results/logs/peakCalling/edd/{{sample}}_summary.txt".format(outdir),
    log:
        "{}results/logs/peakCalling/edd/{{sample}}_callpeak.log".format(outdir),
    params:
        chrom_size=config["params"]["peakCalling"]["chrom_sizes"],
        fdr=config["params"]["peakCalling"]["edd"]["fdr"],
        blacklist=config["resources"]["ref"]["blacklist"],
        extraConfig=config["params"]["peakCalling"]["edd"]["extraParameters"],
        output_dir="{}results/peakCalling/edd/{{sample}}".format(outdir),
    benchmark:
        "{}results/.benchmarks/{{sample}}.edd.benchmark.txt".format(outdir)
    conda:
        "../envs/edd.yaml"
    threads: config["threads"]["edd"]
    shell:
        """
        edd --nprocs {threads} \
            --fdr {params.fdr} \
            --config-file {params.extraConfig} \
            {params.chrom_size} \
            {params.blacklist} \
            {input.treatment} {input.control} \
            {params.output_dir} &> {log} || touch {output.veryBroadPeaks} 
            
        mv {params.output_dir}/log.txt {params.output_dir}/{wildcards.sample}_runlog.txt
        mv {params.output_dir}/{wildcards.sample}_runlog.txt $(dirname {log})
        # Count the number of peaks and create a summary file
        echo "Total Peaks: $(wc -l < {output.veryBroadPeaks})" > {output.summary}
        """
