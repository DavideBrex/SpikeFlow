

ruleorder: merge_lanes_pe > merge_lanes_se


rule merge_lanes_se:
    input:
        sample=get_fastq,
    output:
        temp("{}results/fastq/{{id}}.fastq.gz".format(outdir)),
    log:
        "{}results/logs/pipe-fastqs/{{id}}.log".format(outdir),
    threads: 1
    shell:
        "cat {input.sample} > {output} 2> {log}"


rule merge_lanes_pe:
    input:
        unpack(get_fastq),
    output:
        fw=temp("{}results/fastq/{{id}}_1.fastq.gz".format(outdir)),
        rv=temp("{}results/fastq/{{id}}_2.fastq.gz".format(outdir)),
    log:
        "{}results/logs/pipe-fastqs/{{id}}.log".format(outdir),
    threads: 1
    shell:
        """
        cat {input.fw} >  {output.fw}
        cat {input.rv} >  {output.rv}
        """


if config["trimming"]:

    ruleorder: fastp_pe > fastp_se

    rule fastp_se:
        input:
            sample=get_fastq_trimming,
        output:
            trimmed=temp("{}results/trimmed/{{id}}.fastq.gz".format(outdir)),
            html="{}report/trimmed/{{id}}.html".format(outdir),
            json="{}report/trimmed/{{id}}.json".format(outdir),
        log:
            "{}results/logs/fastp/{{id}}.log".format(outdir),
        message:
            "Trimming {input.sample} with fastp-se."
        conda:
            "../envs/various.yaml"
        params:
            extra=config["params"]["fastp-se"],
        threads: config["threads"]["fastp"]
        shell:
            """
            if ! command -v fastp &> /dev/null
            then
                echo "fastp could not be found"
                exit 1
            fi
            
            fastp --thread {threads} \
                -i {input.sample} \
                -o {output.trimmed} \
                {params.extra} \
                --html {output.html} \
                --json {output.json} \
                2> {log}
            """

    rule fastp_pe:
        input:
            sample=get_fastq_trimming,
        output:
            trimmed=temp(
                [
                    "{}results/trimmed/{{id}}_1.fastq.gz".format(outdir),
                    "{}results/trimmed/{{id}}_2.fastq.gz".format(outdir),
                ]
            ),
            html="{}report/trimmed/{{id}}.html".format(outdir),
            json="{}report/trimmed/{{id}}.json".format(outdir),
        log:
            "{}results/logs/fastp/{{id}}.log".format(outdir),
        message:
            "Trimming {input.sample} with fastp-pe."
        conda:
            "../envs/various.yaml"
        params:
            extra=config["params"]["fastp-pe"],
        threads: config["threads"]["fastp"]
        shell:
            """
            if ! command -v fastp &> /dev/null
            then
                echo "fastp could not be found"
                exit 1
            fi
            fastp --thread {threads} \
                -i {input.sample[0]} \
                -I {input.sample[1]} \
                -o {output.trimmed[0]} \
                -O {output.trimmed[1]} \
                --html {output.html} \
                --json {output.json} \
                {params.extra} \
                2> {log}
            """
