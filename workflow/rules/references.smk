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
        faFile="resources/reference_genome/{assembly}/{assembly}.fa",
    log:
        "{}results/logs/ref/get_reference_{{assembly}}.log".format(outdir),
    params:
        provider="ucsc",  # optional, defaults to ucsc. Choose from ucsc, ensembl, and ncbi
        assemblyRef=config["resources"]["ref"]["assembly"],
    cache: "omit-software"  # mark as eligible for between workflow caching
    conda:
        "../envs/genomepy.yaml"
    shell:
        """
        genomepy plugin enable blacklist
        genomepy install {params.assemblyRef} -g resources/reference_genome \
        --provider {params.provider} >>{log} 2>&1
        """


rule create_bowtie_index_reference:
    input:
        expand(
            "resources/reference_genome/{assembly}/{assembly}.fa",
            assembly=config["resources"]["ref"]["assembly"],
        ),
    output:
        genome=directory("resources/reference_genome/index/"),
    log:
        "{}results/logs/ref/indexing_reference.log".format(outdir),
    message:
        "Creating bowtie index"
    conda:
        "../envs/bowtie.yaml"
    threads: config["threads"]["bowtie"]
    params:
        genome_path=config["resources"]["ref"]["index"],
    cache: True
    shell:
        """
        # Add a condition in the shell script to determine if commands should run
        if [ -z "{params.genome_path}" ]; then
            mkdir resources/reference_genome/index/
            bowtie-build --threads {threads} {input} {output}/index_ref >>{log} 2>&1
        fi
        """


# Rule priority
if config["resources"]["ref"]["index"]:

    ruleorder: return_genome_path > create_bowtie_index_reference

else:

    ruleorder: create_bowtie_index_reference > return_genome_path


# SPIKE IN RULES


rule return_spike_path:
    output:
        genome=directory("resources/spike_genome/index/"),
    log:
        "{}results/logs/ref/return_spike_path.log".format(outdir),
    params:
        genome_path=config["resources"]["ref_spike"]["index_spike"],
    shell:
        """
        if [ -n "{params.genome_path}" ]; then
            mkdir resources/spike_genome/index/
            prefix=$(basename {params.genome_path})
            directory=$(dirname {params.genome_path})
            ln -s $directory/*  {output.genome}
        fi
        """


rule get_spike_genome:
    output:
        faFile="resources/spike_genome/{assemblySpike}/{assemblySpike}.fa",
    log:
        "{}results/logs/ref/get_spike_{{assemblySpike}}.log".format(outdir),
    params:
        provider="ucsc",  # optional, defaults to ucsc. Choose from ucsc, ensembl, and ncbi
        spike_assembly=config["resources"]["ref_spike"]["spike_assembly"],
    cache: "omit-software"  # mark as eligible for between workflow caching
    conda:
        "../envs/genomepy.yaml"
    shell:
        """
        genomepy install {params.spike_assembly} -g resources/spike_genome \
        --provider {params.provider} >>{log} 2>&1
        """


rule create_bowtie_index_spike:
    input:
        expand(
            "resources/spike_genome/{assemblySpike}/{assemblySpike}.fa",
            assemblySpike=config["resources"]["ref_spike"]["spike_assembly"],
        ),
    output:
        genome=directory("resources/spike_genome/index/"),
    log:
        "{}results/logs/ref/indexing_spike.log".format(outdir),
    message:
        "Creating bowtie index for spike-in genome"
    conda:
        "../envs/bowtie.yaml"
    threads: config["threads"]["bowtie_spike"]
    params:
        genome_path=config["resources"]["ref_spike"]["index_spike"],
    cache: True
    shell:
        """
        # Add a condition in the shell script to determine if commands should run
        if [ -z "{params.genome_path}" ]; then
            mkdir resources/spike_genome/index/
            bowtie-build --threads {threads} {input} {output}/index_spike >>{log} 2>&1
        fi
        """


# Rule priority
if config["resources"]["ref_spike"]["index_spike"]:

    ruleorder: return_spike_path > create_bowtie_index_spike

else:

    ruleorder: create_bowtie_index_spike > return_spike_path
