

ruleorder: merge_lanes_pe > merge_lanes_se


rule merge_lanes_se:
    input:
        sample=get_fastq
    output:
        temp("results/fastq/{id}.fastq.gz"),
    log:
        "results/logs/pipe-fastqs/{id}.log",
    threads: 1
    shell:
        "cat {input.sample} > {output} 2> {log}"


rule merge_lanes_pe:
    input:
        unpack(get_fastq)
    output:
        fw=temp("results/fastq/{id}_1.fastq.gz"),
        rv=temp("results/fastq/{id}_2.fastq.gz")
    log:
        "results/logs/pipe-fastqs/{id}.log",
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
            sample=get_fastq_trimming
        output:
            trimmed ="results/trimmed/{id}.fastq.gz",
            html    ="report/trimmed/{id}.html",
            json    ="report/trimmed/{id}.json"
        log:
            "results/logs/fastp/{id}.log"
        params:
            adapters=config["params"]["fastp-se"],
            extra=""
        threads: 3
        wrapper:
            "v2.6.0/bio/fastp"



    rule fastp_pe:
        input:
            sample=get_fastq_trimming
        output:
            trimmed =["results/trimmed/{id}_1.fastq.gz", "results/trimmed/{id}_2.fastq.gz"],
            html    ="report/trimmed/{id}.html",
            json    ="report/trimmed/{id}.json"
        log:
            "results/logs/fastp/{id}.log"
        params:
            adapters=config["params"]["fastp-pe"],
            extra   =""
        threads: 3
        wrapper:
            "v2.6.0/bio/fastp"
        