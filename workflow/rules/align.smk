if config["aligner"] == "bowtie":

    rule align_bowtie:
        input:
            reads=get_reads,
            idx="resources/reference_genome/index/",
        output:
            bam=temp("{}results/bam/{{id}}.tmp.bam".format(outdir)),
            index=temp("{}results/bam/{{id}}.tmp.bam.bai".format(outdir)),
        threads: 8
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
            align="{}results/logs/alignments_bowtie/{{id}}.log".format(outdir),
            rm_dups="{}results/logs/alignments_bowtie/rm_dup/{{id}}.log".format(outdir),
        benchmark:
            "{}results/.benchmarks/{{id}}.bowtie_align.benchmark.txt".format(outdir)
        shell:
            """
            bowtie -p {threads} {params.bowtie} -x {params.index} {params.inputsel} 2> {log.align} \
            | samblaster --removeDups 2> {log.rm_dups} \
            | samtools view -Sb -F 4 - \
            | samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
            samtools index {output.bam}
            """

    # SPIKE alignment
    rule align_bowtie_spike:
        input:
            reads=get_reads,
            idx="resources/spike_genome/index/",
        output:
            bam=temp("{}results/bam_spike/{{id}}_spike.tmp.bam".format(outdir)),
            index=temp("{}results/bam_spike/{{id}}_spike.tmp.bam.bai".format(outdir)),
        threads: 8
        params:
            index=config["resources"]["ref_spike"]["index_spike"]
            if config["resources"]["ref_spike"]["index_spike"] != ""
            else "resources/spike_genome/index/index_spike",
            bowtie=config["params"]["bowtie"]["global"],
            samtools_mem=config["params"]["samtools"]["memory"],
            inputsel=(
                lambda wildcards, input: input.reads
                if len(input.reads) == 1
                else config["params"]["bowtie"]["pe"]
                + " -1 {0} -2 {1}".format(*input.reads)
            ),
        message:
            "SPIKE-IN - Aligning {input} with parameters {params.bowtie}"
        conda:
            "../envs/bowtie.yaml"
        log:
            align="{}results/logs/alignments_bowtie/spike/{{id}}_spike.log".format(outdir),
            rm_dups="{}results/logs/alignments_bowtie/spike/rm_dup/{{id}}_spike.log".format(
                outdir
            ),
        benchmark:
            "{}results/.benchmarks/{{id}}_spike.bowtie_align.benchmark.txt".format(outdir)
        shell:
            """
            bowtie -p {threads} {params.bowtie} -x {params.index} {params.inputsel} 2> {log.align} \
            | samblaster --removeDups 2> {log.rm_dups} \
            | samtools view -Sb -F 4 - \
            | samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
            samtools index {output.bam}
            """

if config["aligner"] == "chromap":

    rule create_chromap_index_reference:
        input:
            ref=rules.get_reference_genome.output,
        output:
            index="resources/reference_genome/index_chromap/chromap_index_ref",
        conda:
            "../envs/various.yaml"
        log:
            "{}results/logs/ref/chromap_create_index_ref.log".format(outdir),
        benchmark:
            "{}results/.benchmarks/chromap_create_index_ref.benchmark.txt".format(outdir)
        shell:
            """
            chromap -i -r {input.ref} -o {output.index} 2> {log}
            """

    rule align_chromap:
        input:
            reads=get_reads,
            idx="resources/reference_genome/index_chromap/chromap_index_ref",
            ref=rules.get_reference_genome.output,
        output:
            bam=temp("{}results/bam/{{id}}.tmp.bam".format(outdir)),
            index=temp("{}results/bam/{{id}}.tmp.bam.bai".format(outdir)),
        threads: 8
        params:
            samtools_mem=config["params"]["samtools"]["memory"],
            inputsel=(
                lambda wildcards, input: "-1 "+str(input.reads)
                if len(input.reads) == 1
                else " -1 {0} -2 {1} -l 2000 ".format(*input.reads)
            ),
        conda:
            "../envs/various.yaml"
        log:
            align="{}results/logs/alignments_chromap/{{id}}.log".format(outdir),
        benchmark:
            "{}results/.benchmarks/{{id}}.chromap_align.benchmark.txt".format(outdir)
        shell:
            """
            chromap -x {input.idx} -r {input.ref} {params.inputsel} \
            --remove-pcr-duplicates --SAM  -t {threads} -o /dev/stdout 2> {log.align} \
            | samtools view -Sb -F 4 - \
            | samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
            samtools index {output.bam}
            """


    rule create_chromap_index_spike:
        input:
            ref=rules.get_spike_genome.output,
        output:
            index="resources/spike_genome/index_chromap/chromap_index_spike",
        conda:
            "../envs/various.yaml"
        log:
            "{}results/logs/ref/chromap_create_index_spike.log".format(outdir),
        benchmark:
            "{}results/.benchmarks/chromap_create_index_spike.benchmark.txt".format(outdir)
        shell:
            """
            chromap -i -r {input.ref} -o {output.index} 2> {log}
            """

    # SPIKE alignment
    rule align_chromaps_spike:
        input:
            reads=get_reads,
            idx="resources/spike_genome/index_chromap/chromap_index_spike",
            ref=rules.get_spike_genome.output,
        output:
            bam=temp("{}results/bam_spike/{{id}}_spike.tmp.bam".format(outdir)),
            index=temp("{}results/bam_spike/{{id}}_spike.tmp.bam.bai".format(outdir)),
        threads: 8
        params:
            samtools_mem=config["params"]["samtools"]["memory"],
            inputsel=(
                lambda wildcards, input: "-1 "+str(input.reads)
                if len(input.reads) == 1
                else " -1 {0} -2 {1} -l 2000 ".format(*input.reads)
            ),
        conda:
            "../envs/various.yaml"
        log:
            align="{}results/logs/alignments_chromap/spike/{{id}}_spike.log".format(outdir),
        benchmark:
            "{}results/.benchmarks/{{id}}_spike.chromap_align.benchmark.txt".format(outdir)
        shell:
            """
            chromap -x {input.idx} -r {input.ref} {params.inputsel} \
            --remove-pcr-duplicates --SAM  -t {threads} -o /dev/stdout 2> {log.align} \
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
