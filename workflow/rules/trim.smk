




rule fastp_se:
    input:
        get_fastq
    output:
        trimmed="results/trimmed/{id}.fastq.gz",
        # html="report/trimmed/{id}.html",
        # json="report/trimmed/{id}.json"
    log:
        "logs/fastp/{id}.log"
    params:
        adapters=config["params"]["fastp-se"],
        extra=""
    threads: 3
    wrapper:
        "v2.6.0/bio/fastp"



rule fastp_pe:
    input:
        get_fastq
    output:
        trimmed=["results/trimmed/{id}.1.fastq.gz", "results/trimmed/{id}.2.fastq.gz"],
        # html="report/trimmed/{id}.html",
        # json="report/trimmed/{id}.json"
    log:
        "logs/fastp/{id}.log"
    params:
        adapters=config["params"]["fastp-pe"],
        extra=""
    threads: 3
    wrapper:
        "v2.6.0/bio/fastp"
    