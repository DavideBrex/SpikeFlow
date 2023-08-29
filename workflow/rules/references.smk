


rule get_reference_genome:
    output:
        "resources/genome.fasta"
    log:
        "logs/ref/get_reference_genome.log"
    params:
        species=config["resources"]["ref"]["species"],
        datatype="dna",
        build=config["resources"]["ref"]["build"],
        release=config["resources"]["ref"]["release"]
    cache: True
    wrapper:
       "v2.6.0/bio/reference/ensembl-sequence"


rule create_bowtie_index:
    input:
        "resources/genome.fasta"
    output:
        multiext(
            "resources/reference_genome/genome",
            ".1.ebwt",
            ".2.ebwt",
            ".3.ebwt",
            ".4.ebwt",
            ".rev.1.ebwt",
            ".rev.2.ebwt",
        ),
    log:
        "logs/ref/indexing_reference.log"
    message:
        "Creating bowtie index"
    conda:
        "../envs/bowtie.yaml"
    threads: 10
    priority:2
    shell:
        """
        bowtie-build --threads {threads} {input} genome
        """



rule get_spike_genome:
    output:
        "resources/ref/genome_spike.fasta"
    log:
        "logs/ref/get_spike_genome.log"
    params:
        species=config["resources"]["ref_spike"]["species"],
        datatype="dna",
        build=config["resources"]["ref_spike"]["build"],
        release=config["resources"]["ref_spike"]["release"],
    cache: True
    wrapper:
        "v2.6.0/bio/reference/ensembl-sequence"

