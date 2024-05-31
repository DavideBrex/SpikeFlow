
rule consensus_peaks:
    input:
        narrowCalled=expand(
            "{}results/peakCallingNorm/{{sample}}_narrowPeaks.narrowPeak".format(
                outdir
            ),
            sample=narrowSamples,
        )
        if config["diffPeakAnalysis"]["useSpikeinCalledPeaks"]
        else expand(
            "{}results/peakCalling/macs2/{{sample}}_peaks.narrowPeak".format(outdir),
            sample=narrowSamples,
        ),
        broadCalled=expand(
            "{}results/peakCallingNorm/{{sample}}_broadPeaks.broadPeak".format(outdir),
            sample=broadSamples,
        )
        if config["diffPeakAnalysis"]["useSpikeinCalledPeaks"]
        else expand(
            "{}results/peakCalling/epic2/{{sample}}_broadPeaks.bed".format(outdir),
            sample=broadSamples,
        ),
    output:
        "{}results/peakCallingNorm/mergedPeaks/{{antibody}}_consensusPeaks.bed".format(
            outdir
        )
        if config["diffPeakAnalysis"]["useSpikeinCalledPeaks"]
        else "{}results/peakCalling/mergedPeaks/{{antibody}}_consensusPeaks.bed".format(
            outdir
        ),
    params:
        min_num_reps=config["diffPeakAnalysis"]["minNumSamples"],
        antibody=lambda w: config["diffPeakAnalysis"]["contrasts"][w.antibody],
        sampleNamesToUse=lambda w: antibody_dict[w.antibody],
    log:
        "{}results/logs/peakCallingNorm/mergedPeaks/{{antibody}}_consensusPeaks.log".format(
            outdir
        )
        if config["diffPeakAnalysis"]["useSpikeinCalledPeaks"]
        else "{}results/logs/peakCalling/mergedPeaks/{{antibody}}_consensusPeaks.log".format(
            outdir
        ),
    threads: 1
    benchmark:
        "{}results/.benchmarks/{{antibody}}_consensusPeaks.benchmark.txt".format(outdir)
    conda:
        "../envs/various.yaml"
    script:
        "../scripts/consensusPeaks.py"


rule count_reads_on_peaks:
    input:
        bamFiles=get_bams_by_antibody,
        consensus_peaks="{outdir}results/peakCallingNorm/mergedPeaks/{{antibody}}_consensusPeaks.bed".format(
            outdir=outdir
        )
        if config["diffPeakAnalysis"]["useSpikeinCalledPeaks"]
        else "{outdir}results/peakCalling/mergedPeaks/{{antibody}}_consensusPeaks.bed".format(
            outdir=outdir
        ),
    output:
        output_tsv="{outdir}results/peakCallingNorm/mergedPeaks/{{antibody}}_readsOnConsensusPeaks.tsv".format(
            outdir=outdir
        )
        if config["diffPeakAnalysis"]["useSpikeinCalledPeaks"]
        else "{outdir}results/peakCalling/mergedPeaks/{{antibody}}_readsOnConsensusPeaks.tsv".format(
            outdir=outdir
        ),
    params:
        joined_bams=lambda w, input: ",".join(input.bamFiles),
        map_qual=config["params"]["bowtie2"]["map_quality"],
    log:
        "{}results/logs/peakCallingNorm/mergedPeaks/{{antibody}}_countReadsOnPeaks.log".format(
            outdir
        )
        if config["diffPeakAnalysis"]["useSpikeinCalledPeaks"]
        else "{}results/logs/peakCalling/mergedPeaks/{{antibody}}_countReadsOnPeaks.log".format(
            outdir
        ),
    benchmark:
        "{}results/.benchmarks/{{antibody}}_countReadsOnPeaks.benchmark.txt".format(
            outdir
        )
    conda:
        "../envs/various.yaml"
    shell:
        """
        python3 workflow/scripts/frag_count.py -b {input.consensus_peaks} \
            -i {params.joined_bams} \
            -o {output.output_tsv} \
            --mapq {params.map_qual} &> {log}
        """


rule peakAnnot_singleRep:
    input:
        peaks=get_singleRep_peaks,
    output:
        annotations="{outdir}results/peakCalling/peakAnnot/{{sample}}_annot.txt".format(
            outdir=outdir
        ),
        promoBed="{outdir}results/peakCalling/peakAnnot/{{sample}}_promoPeaks.bed".format(
            outdir=outdir
        ),
        distalBed="{outdir}results/peakCalling/peakAnnot/{{sample}}_distalPeaks.bed".format(
            outdir=outdir
        ),
        annotInfo="{outdir}results/logs/peakCalling/peakAnnot/{{sample}}_annotInfo.txt".format(
            outdir=outdir
        ),
    params:
        before=config["params"]["peaksAnnotation"]["promoter"]["upstream"],
        after=config["params"]["peaksAnnotation"]["promoter"]["downstream"],
        genome=config["resources"]["ref"]["assembly"],
    log:
        "{outdir}results/logs/peakCalling/peakAnnot/{{sample}}_singleRep.log".format(
            outdir=outdir
        ),
    benchmark:
        "{}results/.benchmarks/{{sample}}.singleRep.annotatePeaks.benchmark.txt".format(
            outdir
        )
    conda:
        "../envs/Renv.yaml"
    script:
        "../scripts/PeakAnnot.R"


rule peakAnnot_singleRep_normPeaks:
    input:
        peaks=get_singleRep_peaksnorm,
    output:
        annotations="{outdir}results/peakCallingNorm/peakAnnot/{{sample}}_annot.txt".format(
            outdir=outdir
        ),
        promoBed="{outdir}results/peakCallingNorm/peakAnnot/{{sample}}_promoPeaks.bed".format(
            outdir=outdir
        ),
        distalBed="{outdir}results/peakCallingNorm/peakAnnot/{{sample}}_distalPeaks.bed".format(
            outdir=outdir
        ),
        annotInfo="{outdir}results/logs/peakCallingNorm/peakAnnot/{{sample}}_annotInfo.txt".format(
            outdir=outdir
        ),
    params:
        before=config["params"]["peaksAnnotation"]["promoter"]["upstream"],
        after=config["params"]["peaksAnnotation"]["promoter"]["downstream"],
        genome=config["resources"]["ref"]["assembly"],
    log:
        "{outdir}results/logs/peakCallingNorm/peakAnnot/{{sample}}_singleRep.log".format(
            outdir=outdir
        ),
    benchmark:
        "{}results/.benchmarks/{{sample}}.singleRep.annotatePeaks.benchmark.txt".format(
            outdir
        )
    conda:
        "../envs/Renv.yaml"
    script:
        "../scripts/PeakAnnot.R"
