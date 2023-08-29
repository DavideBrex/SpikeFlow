



if config["aligner"] == "bowtie":


    rule align_bowtie:
        input:
            get_reads
        output:
            bam   = temp("results/bam/{id}.bam"),
            index = temp("results/bam/{id}.bam.tmp.bai")
        threads:
            8
        params:
            index  	     = config["resources"]["ref"]["index"],
            bowtie 	     = config["params"]["bowtie"]["global"],
            samtools_mem = config["params"]["samtools"]["memory"],
            input  	     = (
                            lambda wildcards, input: "{}".format(*input)
                            if len(input) == 1
                            else config["params"]["bowtie"]["pe"] + " -1 {} -2 {}".format(*input)
                            )
        message:
            "Aligning {input} with parameters {params.bowtie}"
        conda:
            "../envs/bowtie.yaml"
        log:
            align   = "results/00log/alignments/{id}.log",
            rm_dups = "results/00log/alignments/rm_dup/{id}.log",
        benchmark:
            "results/.benchmarks/{id}.align.benchmark.txt"
        shell:
            """
            bowtie -p {threads} {params.bowtie} {params.index} {params.input} 2> {log.align} \
            | samblaster --removeDups 2> {log.rm_dups} \
            | samtools view -Sb -F 4 - \
            | samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
            samtools index {output.bam}
            """