rule macs2_callNarrowPeak:
    input:
        treatment="{}results/bam/{{sample}}_ref.sorted.bam".format(outdir),
        control=lambda w: "{}results/bam/{}_ref.sorted.bam".format(
            outdir, sample_to_input[w.sample]
        ),
    output:
        multiext(
            "{}results/peakCalling/macs2/{{sample}}".format(outdir),
            "_peaks.xls",
            "_peaks.narrowPeak",
        ),
    log:
        "{}results/logs/peakCalling/macs2/{{sample}}_callpeak.log".format(outdir),
    params:
        isPairedEnd=lambda w: (
            "--format BAM --nomodel" if is_single_end(w.sample) else "--format BAMPE"
        ),
        gsize=config["params"]["deeptools"]["effective_genome_length"],
        qval=config["params"]["peakCalling"]["macs2"]["qvalue"],
        otherParams=config["params"]["peakCalling"]["macs2"]["extraOptions"],
        output_dir=lambda w, output: os.path.dirname(output[0]),
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


rule macs2_callNormPeaks_narrow:
    input:
        treatment="{}results/bam/{{sample}}_ref.sorted.bam".format(outdir),
        treatmentIndex="{}results/bam/{{sample}}_ref.sorted.bam.bai".format(outdir),
        control=lambda w: "{}results/bam/{}_ref.sorted.bam".format(
            outdir, sample_to_input[w.sample]
        ),
        controlIndex=lambda w: "{}results/bam/{}_ref.sorted.bam.bai".format(
            outdir, sample_to_input[w.sample]
        ),
        logFile="{}results/logs/spike/{{sample}}.normFactor".format(outdir),
        logFileInput=lambda wildcards: (
            "{}results/logs/spike/{}.normFactor".format(
                outdir, sample_to_input[wildcards.sample]
            )
            if not pd.isna(sample_to_input[wildcards.sample])
            else "{}results/logs/spike/{{sample}}.normFactor".format(outdir)
        ),
    output:
        multiext(
            "{}results/peakCallingNorm/{{sample}}".format(outdir),
            ".treat.pileup.SpikeIn_scaled.bigWig",
            "_narrowPeaks.narrowPeak",
        ),
    log:
        "{}results/logs/peakCallingNorm/{{sample}}_callNormPeak_narrow.log".format(
            outdir
        ),
    params:
        scaleFactors=lambda wildcards: spiker_normalization_factor(wildcards),
        output_prefix=lambda w, output: output[0].split(os.extsep)[0],
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
            --q-peak {params.qvalue} \
            --genome-size {params.gsize} \
            --bw  \
            --cleanup \
            -o {params.output_prefix} &> {log}

        #we modify the output peak names in the narrow (first case) or broad (else) files to match the naming convention for normal peak calling
        sed -i -e 's_results/peakCallingNorm/__' -e 's/\\(Peak\\)\\([0-9]\\)/\\1_\\2/g' {params.output_prefix}.narrowPeak 2>> {log}
        #change file names 
        mv {params.output_prefix}.narrowPeak {params.output_prefix}_narrowPeaks.narrowPeak  2>> {log}
        #we remove other files that are not needed
        rm {params.output_prefix}.control.*   {params.output_prefix}.treat.pileup.bdg \
        {params.output_prefix}.treat.pileup.SpikeIn_scaled.bdg {params.output_prefix}.treat.predict_d.r 2>> {log}
        """


rule macs2_callNormPeaks_broad:
    input:
        treatment="{}results/bam/{{sample}}_ref.sorted.bam".format(outdir),
        control=lambda w: "{}results/bam/{}_ref.sorted.bam".format(
            outdir, sample_to_input[w.sample]
        ),
        logFile="{}results/logs/spike/{{sample}}.normFactor".format(outdir),
        logFileInput=lambda wildcards: (
            "{}results/logs/spike/{}.normFactor".format(
                outdir, sample_to_input[wildcards.sample]
            )
            if not pd.isna(sample_to_input[wildcards.sample])
            else "{}results/logs/spike/{{sample}}.normFactor".format(outdir)
        ),
    output:
        multiext(
            "{}results/peakCallingNorm/{{sample}}".format(outdir),
            ".treat.pileup.SpikeIn_scaled.bigWig",
            "_broadPeaks.broadPeak",
        ),
    log:
        "{}results/logs/peakCallingNorm/{{sample}}_callNormPeak_narrow.log".format(
            outdir
        ),
    params:
        scaleFactors=lambda wildcards: spiker_normalization_factor(wildcards),
        output_prefix=lambda w, output: output[0].split(os.extsep)[0],
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
            --broad \
            --q-peak {params.qvalue} \
            --genome-size {params.gsize} \
            --bw  \
            --cleanup \
            -o {params.output_prefix} &> {log}
        #we modify the output peak names in the narrow (first case) or broad (else) files to match the naming convention for normal peak calling
        sed -i -e 's_results/peakCallingNorm/__' -e 's/\\(Region\\)\\([0-9]\\)/\\1_\\2/g' {params.output_prefix}.broadPeak 2>> {log}
        #change file name
        mv {params.output_prefix}.broadPeak {params.output_prefix}_broadPeaks.broadPeak 2>> {log}
        #we remove other files that are not needed
        rm {params.output_prefix}.control.*   {params.output_prefix}.treat.pileup.bdg \
        {params.output_prefix}.treat.pileup.SpikeIn_scaled.bdg {params.output_prefix}.treat.predict_d.r 2>> {log}
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
