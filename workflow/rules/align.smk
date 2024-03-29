if config["aligner"] == "bowtie":

    rule align_bowtie:
        input:
            reads=get_reads,
            idx="resources/reference_genome/index/",
        output:
            bam=temp("{}results/bam/{{id}}.tmp.bam".format(outdir)),
            index=temp("{}results/bam/{{id}}.tmp.bam.bai".format(outdir)),
        threads: config["threads"]["bowtie"]
        params:
            index=config["resources"]["ref"]["index"]
            if config["resources"]["ref"]["index"] != ""
            else "resources/reference_genome/index/index_ref",
            bowtie=config["params"]["bowtie"]["global"],
            samtools_mem=config["params"]["samtools"]["memory"],
            inputsel=(
                lambda wildcards, input: input.reads
                if len(input.reads) == 1
                else config["params"]["bowtie"]["pe"]
                + " -1 {0} -2 {1}".format(*input.reads)
            ),
        message:
            "Aligning {input} with parameters {params.bowtie}"
        conda:
            "../envs/bowtie.yaml"
        log:
            align="{}results/logs/alignments/{{id}}.log".format(outdir),
            rm_dups="{}results/logs/alignments/rm_dup/{{id}}.log".format(outdir),
        benchmark:
            "{}results/.benchmarks/{{id}}.align.benchmark.txt".format(outdir)
        shell:
            """
            bowtie -p {threads} {params.bowtie} -x {params.index} {params.inputsel} 2> {log.align} \
            | samblaster --removeDups 2> {log.rm_dups} \
            | samtools view -Sb -F 4 - \
            | samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
            samtools index {output.bam}
            """




rule clean_spike:
    input:
        sample_ref="{}results/bam/{{id}}.tmp.bam".format(outdir),
        sample_spike="{}results/bam_spike/{{id}}_spike.tmp.bam".format(outdir),
        sample_ref_index="{}results/bam/{{id}}.tmp.bam.bai".format(outdir),
        sample_spike_index="{}results/bam_spike/{{id}}_spike.tmp.bam.bai".format(outdir),
    output:
        sample_ref="{}results/bam/{{id}}.clean.bam".format(outdir),
        sample_spike="{}results/bam_spike/{{id}}_spike.clean.bam".format(outdir),
    conda:
        "../envs/various.yaml"
    log:
        "{}results/logs/spike/{{id}}.removeSpikeDups".format(outdir),
    benchmark:
        "{}results/.benchmarks/{{id}}.clean_spike.benchmark.txt".format(outdir)
    script:
        "../scripts/remove_spikeDups.py"



