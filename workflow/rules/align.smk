



if config["aligner"] == "bowtie":


    rule align_bowtie:
        input:
            reads=get_reads,
            idx= "resources/reference_genome/genome.1.ebwt"
            
        output:
            bam   = temp("results/bam/{id}.bam"),
            index = temp("results/bam/{id}.bam.tmp.bai")
        threads:
            8
        params:
            bowtie 	     = config["params"]["bowtie"]["global"],
            samtools_mem = config["params"]["samtools"]["memory"],
            inputsel  	     = (
                            lambda wildcards, input: input.reads
                            if len(input.reads) == 1
                            else config["params"]["bowtie"]["pe"] + " -1 {0} -2 {1}".format(*input.reads)
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
            bowtie -p {threads} {params.bowtie} -x {input.idx} {params.inputsel} 2> {log.align} \
            | samblaster --removeDups 2> {log.rm_dups} \
            | samtools view -Sb -F 4 - \
            | samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
            samtools index {output.bam}
            """

