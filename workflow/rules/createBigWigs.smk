rule bam2bigwig:
    input: 
        bamfile=lambda wildcards: "results/bam/{}.bam".format(wildcards.id)
        if not is_spike(wildcards.id),
        blacklist=config["resources"]["ref"]["blacklist"],
    output:
        out='/results/bigWigs/{sample}.bw',
    params: 
        extra='--scaleFactor RPKM'
    log: 'results/logs/bam2bigwig/{sample}.log'
    threads: 5
    wrapper:
        "v2.6.0/bio/deeptools/bamcoverage"


rule bam2bigwig_spikein:
    input: 
        bamfile=lambda wildcards: "results/bam/{}.bam".format(wildcards.id)
        if not is_spike(wildcards.id),
        blacklist=config["resources"]["ref"]["blacklist"],
    output:
        out='/results/bigWigs/{sample}_spike.bw',
    params: 
        extra = '--scaleFactor '
    log: 'results/logs/bam2bigwig/{sample}.log'
    threads: 5
    wrapper: "v2.6.0/bio/deeptools/bamcoverage"