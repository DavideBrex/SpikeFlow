rule bam2bigwig_general:
    input:
        bam="results/bam/{id}.clean.bam",
        logFile="results/logs/spike/{id}.removeSpikeDups",
        logFileInput=lambda wildcards: "results/logs/spike/{}.removeSpikeDups".format(
            sample_to_input[wildcards.id]
        )
        if not pd.isna(sample_to_input[wildcards.id])
        else "results/logs/spike/{id}.removeSpikeDups",
        blacklist=config["resources"]["ref"]["blacklist"],
    output:
        out="results/bigWigs/{id}.bw",
    params:
        extra=lambda wildcards: normalization_factor(wildcards),
        effective_genome_size=config["params"]["deeptools"]["effective_genome_length"],
    log:
        "results/logs/bam2bigwig/{id}.log",
    threads: 5
    wrapper:
        "v2.6.0/bio/deeptools/bamcoverage"


# ruleorder: clean_spike > bam2bigwig_general
