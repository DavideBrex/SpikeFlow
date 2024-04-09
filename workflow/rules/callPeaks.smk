rule macs2_callNarrowPeak:
    input:
        treatment="{}results/bam/{{sample}}_ref.sorted.bam".format(outdir),
        control=lambda w: "{}results/bam/{}_ref.sorted.bam".format(
            outdir, sample_to_input[w.sample]
        ),
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext(
            "{}results/peakCalling/macs2/{{sample}}".format(outdir),
            "_peaks.xls",
            "_peaks.narrowPeak",
        ),
    log:
        "{}results/logs/peakCalling/macs2/{{sample}}_callpeak.log".format(outdir),
    params:
        isPairedEnd=lambda w: "--format BAM --nomodel"
        if is_single_end(w.sample)
        else "--format BAMPE",
        gsize=config["params"]["deeptools"]["effective_genome_length"],
        qval=config["params"]["peakCalling"]["macs2"]["qvalue"],
        otherParams=config["params"]["peakCalling"]["macs2"]["extraOptions"],
        output_dir="{}results/peakCalling/macs2/".format(outdir),
    benchmark:
        "{}results/.benchmarks/{{sample}}.macs2.benchmark.txt".format(outdir)
    conda:
        "../envs/various.yaml"
    shell:
        """
        macs2 callpeak \
            --treatment {input.treatment} \
            --control {input.control} \
            --gsize {params.gsize} \
            --outdir {params.output_dir} \
            --name {wildcards.sample} \
            --qvalue {params.qval} \
            {params.isPairedEnd} \
            {params.otherParams} &> {log}
        """


rule macs2_callNormPeaks:
    input:
        treatment="{}results/bam/{{sample}}_ref.sorted.bam".format(outdir),
        control=lambda w: "{}results/bam/{}_ref.sorted.bam".format(
            outdir, sample_to_input[w.sample]
        ),
        logFile="{}results/logs/spike/{{sample}}.normFactor".format(outdir),
        logFileInput=lambda wildcards: "{}results/logs/spike/{}.normFactor".format(
            outdir, sample_to_input[wildcards.sample]
        )
        if not pd.isna(sample_to_input[wildcards.sample])
        else "{}results/logs/spike/{{sample}}.normFactor".format(outdir),
    output:
        multiext(
            "{}results/peakCallingNorm/{{sample}}".format(outdir),
            ".control.pileup.max.bdg",
            ".control.pileup.max.bigWig",
            ".treat.pileup.bdg",
            ".treat.pileup.SpikeIn_scaled.bdg",
            ".treat.pileup.SpikeIn_scaled.bigWig",
    ),
    log:
        "{}results/logs/peakCallingNorm/{{sample}}_callNormPeak.log".format(outdir),
    params:
        scaleFactors=lambda wildcards: spiker_normalization_factor(wildcards),
        output_dir="{}results/peakCallingNorm/{{sample}}".format(outdir),
        peak_type=lambda wildcards: check_peak_type(wildcards),
        gsize=config["params"]["deeptools"]["effective_genome_length"],
        qvalue=config["params"]["peakCalling"]["macs2"]["qvalue"],
    benchmark:
        "{}results/.benchmarks/{{sample}}.spiker.benchmark.txt".format(outdir)
    conda:
        "../envs/various.yaml"
    shell:
        """
        spiker.py -t {input.treatment} \
            -c {input.control} \
            --spikeIn {params.scaleFactors} \
            {params.peak_type} \
            --q-peak {params.qvalue} \
            --genome-size {params.gsize} \
            --bw  \
            --cleanup \
            -o {params.output_dir} &> {log}
        """



rule epic2_callBroadPeaks:
    input:
        treatment="{}results/bam/{{sample}}_ref.sorted.bam".format(outdir),
        control=lambda w: "{}results/bam/{}_ref.sorted.bam".format(
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
        treatment="{}results/bam/{{sample}}_ref.sorted.bam".format(outdir),
        control=lambda w: "{}results/bam/{}_ref.sorted.bam".format(
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
