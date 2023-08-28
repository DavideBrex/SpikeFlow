


rule get_reference_genome:
    output:
        "resources/ref/genome.fasta"
    log:
        "logs/ref/get_reference_genome.log"
    params:
        species=config["resources"]["ref"]["species"],
        datatype="dna",
        build=config["resources"]["ref"]["build"],
        release=config["resources"]["ref"]["release"],
        chromosome=config["resources"]["ref"]["chromosome"]
    cache: True
    wrapper:
       "v2.6.0/bio/reference/ensembl-sequence"

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
        chromosome=config["resources"]["ref_spike"]["chromosome"]
    cache: True
    wrapper:
        "v2.6.0/bio/reference/ensembl-sequence"



rule 