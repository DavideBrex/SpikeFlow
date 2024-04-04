# if index file for reference genome does not exist, create it


rule return_genome_path:
    output:
        genome=directory("resources/reference_genome/index/"),
    log:
        "{}results/logs/ref/return_genome_path.log".format(outdir),
    params:
        genome_path=config["resources"]["ref"]["index"],
    shell:
        """
        if [ -n "{params.genome_path}" ]; then
            mkdir resources/reference_genome/index/
            prefix=$(basename {params.genome_path})
            directory=$(dirname {params.genome_path})
            ln -s $directory/*  {output.genome}
        fi
        """


rule get_reference_genome:
    output:
        faFile="resources/reference_genome/{assembly}.fa",
    log:
        "{}results/logs/ref/get_reference_{{assembly}}.log".format(outdir),
    params:
        provider="ucsc",  # optional, defaults to ucsc. Choose from ucsc, ensembl, and ncbi
        assembly=config["resources"]["ref"]["assembly"],
    cache: "omit-software"  # mark as eligible for between workflow caching
    priority: 10
    conda:
        "../envs/various.yaml"
    script:
        "../scripts/download_genome.py"


rule create_bowtie_index:
    input:
        refGenome=expand(
            "resources/reference_genome/{assembly}.fa",
            assembly=config["resources"]["ref"]["assembly"],
        ),
        spikeGenome=expand(
            "resources/spike_genome/{assemblySpike}.fa",
            assemblySpike=config["resources"]["ref_spike"]["spike_assembly"],
        ),
    output:
        genome=directory("resources/reference_genome/index/"),
    log:
        "{}results/logs/ref/indexing_reference.log".format(outdir),
    message:
        "Creating bowtie index"
    threads: config["threads"]["bowtie2"]
    params:
        genome_path=config["resources"]["ref"]["index"],
    cache: True
    shell:
        """
        # Add a condition in the shell script to determine if commands should run
        if [ -z "{params.genome_path}" ]; then
            mkdir resources/reference_genome/index/

            # add flag to Exogenous genome chromosome names to distinguish from endogenous genome
            awk '/^>/ {{sub(/^>/, ">EXO_")}} 1' {input.spikeGenome} > resources/spike_genome/exo.tmp.fa

            # merge endogenous and exogenous genomes
            cat {input.refGenome} resources/spike_genome/exo.tmp.fa > resources/reference_genome/mergedGenome.fa

            rm resources/spike_genome/exo.tmp.fa

            #build index
            bowtie2-build --threads {threads} resources/reference_genome/mergedGenome.fa {output}/index_ref >>{log} 2>&1
        fi
        """


# Rule priority
if config["resources"]["ref"]["index"]:

    ruleorder: return_genome_path > create_bowtie_index

else:

    ruleorder: create_bowtie_index > return_genome_path


# SPIKE IN RULES



rule get_spike_genome:
    output:
        faFile="resources/spike_genome/{assemblySpike}.fa",
    log:
        "{}results/logs/ref/get_spike_{{assemblySpike}}.log".format(outdir),
    params:
        provider="ucsc",  # optional, defaults to ucsc. Choose from ucsc, ensembl, and ncbi
        assembly=config["resources"]["ref_spike"]["spike_assembly"],
    cache: "omit-software"  # mark as eligible for between workflow caching
    priority: 10
    conda:
        "../envs/various.yaml"
    script:
        "../scripts/download_genome.py"

