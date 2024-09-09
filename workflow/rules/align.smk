if config["aligner"] == "bowtie2":

    rule align_bowtie2:
        input:
            reads=get_reads,
            idx="resources/reference_genome/index/",
        output:
            bam=temp("{}results/bam/{{id}}.tmp.bam".format(outdir)),
            index=temp("{}results/bam/{{id}}.tmp.bam.bai".format(outdir)),
        threads: config["threads"]["bowtie2"]
        params:
            index=(
                config["resources"]["ref"]["index"]
                if config["resources"]["ref"]["index"] != ""
                else "resources/reference_genome/index/index_ref"
            ),
            bowtie2=config["params"]["bowtie2"]["global"],
            samtools_mem=config["params"]["samtools"]["memory"],
            inputsel=(
                lambda wildcards, input: (
                    " -U {0} ".format(input.reads)
                    if len(input.reads) == 1
                    else config["params"]["bowtie2"]["pe"]
                    + " -1 {0} -2 {1}".format(*input.reads)
                )
            ),
            forSamblaster=(
                lambda wildcards, input: (
                    " --ignoreUnmated ".format(input.reads)
                    if len(input.reads) == 1
                    else ""
                )
            ),
        message:
            "Aligning {input} with parameters {params.bowtie2}"
        conda:
            "../envs/bowtie2.yaml"
        log:
            align="{}results/logs/alignments/{{id}}.log".format(outdir),
            rm_dups="{}results/logs/alignments/rm_dup/{{id}}.log".format(outdir),
        benchmark:
            "{}results/.benchmarks/{{id}}.align.benchmark.txt".format(outdir)
        shell:
            """
            bowtie2 -p {threads} {params.bowtie2} -x {params.index} {params.inputsel} 2> {log.align} \
            | samblaster {params.forSamblaster} --removeDups 2> {log.rm_dups} \
            | samtools view -Sb -F 4 - \
            | samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
            samtools index {output.bam}
            """


rule split_bam:
    input:
        sample_mixed="{}results/bam/{{id}}.tmp.bam".format(outdir),
        sample_mixed_index="{}results/bam/{{id}}.tmp.bam.bai".format(outdir),
    output:
        sample_ref="{}results/bam/{{id}}_ref.sorted.bam".format(outdir),
        sample_spike="{}results/bam/{{id}}_spike.sorted.bam".format(outdir),
        sample_ref_bai="{}results/bam/{{id}}_ref.sorted.bam.bai".format(outdir),
        sample_spike_bai="{}results/bam/{{id}}_spike.sorted.bam.bai".format(outdir),
    threads: config["threads"]["samtools"]
    params:
        outprefix=lambda w, input: input.sample_mixed.split(os.extsep)[0],
        map_qual=config["params"]["bowtie2"]["map_quality"],
    conda:
        "../envs/bowtie2.yaml"
    log:
        readsInfo="{}results/logs/spike/{{id}}.removeSpikeDups".format(outdir),
        otherLog="{}results/logs/spike/{{id}}.split_bam.log".format(outdir),
    benchmark:
        "{}results/.benchmarks/{{id}}.split_bam.benchmark.txt".format(outdir)
    script:
        "../scripts/split_bam.py"
