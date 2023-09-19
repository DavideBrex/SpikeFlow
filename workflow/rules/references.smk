

#if index file for reference genome does not exist, create it

rule return_genome_path:
    output:
        genome = directory("resources/reference_genome/genome/")
    log:
        "results/logs/ref/return_genome_path.log"
    params:
        genome_path=config["resources"]["ref"]["index"]
    shell:
        """
        if [ -n "{params.genome_path}" ]; then
            mkdir resources/reference_genome/genome/
            prefix=$(basename {params.genome_path})
            directory=$(dirname {params.genome_path})
            ln -s $directory/*  {output.genome}
        fi
        """


rule get_reference_genome:
    output:
        "resources/reference_genome/genome.fasta"
    log:
        "results/logs/ref/get_reference_genome.log"
    params:
        species=config["resources"]["ref"]["species"],
        datatype="dna",
        build=config["resources"]["ref"]["build"],
        release=config["resources"]["ref"]["release"]
    cache: True
    wrapper:
        "v2.6.0/bio/reference/ensembl-sequence"


rule create_bowtie_index_reference:
    input:
        "resources/reference_genome/genome.fasta"
    output:
        genome=directory("resources/reference_genome/genome/")
    log:
        "results/logs/ref/indexing_reference.log"
    message:
        "Creating bowtie index"
    conda:
        "../envs/bowtie.yaml"
    threads: 10
    params:
        genome_path=config["resources"]["ref"]["index"]
    cache: True
    shell: 
        """
        # Add a condition in the shell script to determine if commands should run
        if [ -z "{params.genome_path}" ]; then
            bowtie-build --threads {threads} {input} {output}/genome
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
        genome=directory("resources/spike_genome/genome/")
    log:
        "results/logs/ref/return_spike_path.log"
    params:
        genome_path=config["resources"]["ref_spike"]["index_spike"]
    shell:
        """
        if [ -n "{params.genome_path}" ]; then
            mkdir resources/spike_genome/genome/
            prefix=$(basename {params.genome_path})
            directory=$(dirname {params.genome_path})
            ln -s $directory/*  {output.genome}
        fi
        """

rule get_spike_genome:
    output:
        "resources/spike_genome/genome.fasta"
    log:
        "results/logs/ref/get_spike_genome.log"
    params:
        species=config["resources"]["ref_spike"]["species"],
        datatype="dna",
        build=config["resources"]["ref_spike"]["build"],
        release=config["resources"]["ref_spike"]["release"],
    cache: True
    wrapper:
        "v2.6.0/bio/reference/ensembl-sequence"

rule create_bowtie_index_spike:
    input:
        "resources/spike_genome/genome.fasta"
    output:
        genome=directory("resources/spike_genome/genome/")
    log:
        "results/logs/ref/indexing_spike.log"
    message:
        "Creating bowtie index for spike-in genome"
    conda:
        "../envs/bowtie.yaml"
    threads: 10
    params:
        genome_path=config["resources"]["ref_spike"]["index_spike"]
    cache: True
    shell: 
        """
        # Add a condition in the shell script to determine if commands should run
        if [ -z "{params.genome_path}" ]; then
            bowtie-build --threads {threads} {input} {output}/genome
        fi
        """
        

# Rule priority
if config["resources"]["ref_spike"]["index_spike"]:
    ruleorder: return_spike_path > create_bowtie_index_spike
else:
    ruleorder: create_bowtie_index_spike > return_spike_path
