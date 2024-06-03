rule calculate_norm_factors:
    input:
        logFiles=expand(
            "{}results/logs/spike/{{id}}.removeSpikeDups".format(outdir),
            id=sample_to_input.keys(),
        ),
    output:
        expand(
            "{}results/logs/spike/{{id}}.normFactor".format(outdir),
            id=sample_to_input.keys(),
        ),
    params:
        typeOfNorm=config["normalization_type"],
        sampleToInput=sample_to_input,
        antibody_dict=antibody_dict,
        outdir=outdir,
    conda:
        "../envs/various.yaml"
    log:
        "{}results/logs/normFactors/calcuteNormFactors.log".format(outdir),
    script:
        "../scripts/computeNormFactors.py"


rule bam2bigwig_general:
    input:
        bam="{}results/bam/{{id}}_ref.sorted.bam".format(outdir),
        bamIndex="{}results/bam/{{id}}_ref.sorted.bam.bai".format(outdir),
        logFile="{}results/logs/spike/{{id}}.normFactor".format(outdir),
        blacklist=config["resources"]["ref"]["blacklist"],
    output:
        out="{}results/bigWigs/{{id}}.bw".format(outdir),
    params:
        extra=lambda w, input: normalization_factor(w, input.logFile),
        effective_genome_size=config["params"]["deeptools"]["effective_genome_length"],
    message:
        "Generating bigwig file for {input.bam} using bamCoverage"
    conda:
        "../envs/qc.yaml"
    log:
        "{}results/logs/bam2bigwig/{{id}}.log".format(outdir),
    threads: config["threads"]["generateBigWig"]
    benchmark:
        "{}results/.benchmarks/{{id}}.bigwigs.benchmark.txt".format(outdir)
    shell:
        """
        bamCoverage --blackListFileName {input.blacklist} \
                    {params.extra} \
                    --numberOfProcessors {threads} \
                    --effectiveGenomeSize {params.effective_genome_size} \
                    --bam {input.bam} \
                    --outFileName {output.out} \
                    --outFileFormat bigwig \
                    > {log} 2>&1
        """
