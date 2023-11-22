rule bam2bigwig_general:
    input:
        bam="{}results/bam/{{id}}.clean.bam".format(outdir),
        logFile="{}results/logs/spike/{{id}}.removeSpikeDups".format(outdir),
        logFileInput=lambda wildcards: "{}results/logs/spike/{}.removeSpikeDups".format(
            outdir, sample_to_input[wildcards.id]
        )
        if not pd.isna(sample_to_input[wildcards.id])
        else "{}results/logs/spike/{{id}}.removeSpikeDups".format(outdir),
        blacklist=config["resources"]["ref"]["blacklist"],
    output:
        out="{}results/bigWigs/{{id}}.bw".format(outdir),
    params:
        extra=lambda wildcards: normalization_factor(wildcards),
        effective_genome_size=config["params"]["deeptools"]["effective_genome_length"],
    log:
        "{}results/logs/bam2bigwig/{{id}}.log".format(outdir),
    threads: 5
    benchmark:
        "{}results/.benchmarks/{{id}}.bigwigs.benchmark.txt".format(outdir)
    wrapper:
        "v2.6.0/bio/deeptools/bamcoverage"


# ruleorder: clean_spike > bam2bigwig_general
