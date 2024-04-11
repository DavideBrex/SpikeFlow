
rule merge_rep_peaks:
    input:
        get_replicate_peaks,
    output:
        out="{outdir}results/peakCalling/mergedPeaks/{{unique_rep}}_merged_optimal.bed".format(
            outdir=outdir
        ),
    params:
        outputPrefix=lambda w, output: output.out[:-12],
        size=config["params"]["peakMerge"]["chipr"]["size"],
    log:
        "{outdir}results/logs/peakCalling/mergedPeaks/{{unique_rep}}.log".format(
            outdir=outdir
        ),
    threads: 1
    benchmark:
        "{}results/.benchmarks/{{unique_rep}}.merge_peaks.benchmark.txt".format(outdir)
    conda:
        "../envs/various.yaml"
    shell:
        """
        chipr -i {input} -o {params.outputPrefix} -m 2 -s {params.size} &> {log}
        mv {params.outputPrefix}_log.txt $(dirname {log})
        mv {params.outputPrefix}_all.bed $(dirname {log})
        """

rule consensus_peaks:
    input:
        expand(
            "{}results/peakCalling/macs2/{{sample}}_peaks.narrowPeak".format(
                outdir
            ),
            sample=narrowSamples)
    output:
        "{outdir}results/peakCalling/mergedPeaks/consensusPeaks.bed".format(outdir=outdir),
    params:
        min_num_reps=2,
    log:
        "{outdir}results/logs/peakCalling/mergedPeaks/consensusPeaks.log".format(outdir=outdir),
    threads: 1
    benchmark:
        "{}results/.benchmarks/consensusPeaks.benchmark.txt".format(outdir)
    conda:
        "../envs/various.yaml"
    script:
        "../scripts/consensusPeaks.py"


rule count_reads_on_peaks:
    input:
        bamFiles=expand(
            "{}results/bam/{{sample}}_ref.sorted.bam".format(
                outdir
            ),
            sample=narrowSamples
        ),
        consensus_peaks="{outdir}results/peakCalling/mergedPeaks/consensusPeaks.bed".format(outdir=outdir),
    output:
        output_tsv="{outdir}results/peakCalling/mergedPeaks/readsOnConsensusPeaks.tsv".format(outdir=outdir),
    params:
        joined_bams=lambda w, input: ",".join(input.bamFiles),
        map_qual=config["params"]["bowtie2"]["map_quality"],
    log:
        "{}results/logs/peakCalling/mergedPeaks/countReadsOnPeaks.log".format(outdir),
    benchmark:
        "{}results/.benchmarks/countReadsOnPeaks.benchmark.txt".format(outdir)
    conda:
        "../envs/various.yaml"
    shell:
        """
        python3 workflow/scripts/frag_count.py -b {input.consensus_peaks} \
            -i {params.joined_bams} \
            -o {output.output_tsv} \
            --mapq {params.map_qual} &> {log}
        """



rule merge_rep_peaks_edd:
    input:
        get_replicate_peaks_edd,
    output:
        "{outdir}results/peakCalling/mergedPeaks/{{unique_rep}}_merged_intersected.bed".format(
            outdir=outdir
        ),
    log:
        "{outdir}results/logs/peakCalling/mergedPeaks/{{unique_rep}}.log".format(
            outdir=outdir
        ),
    threads: 1
    benchmark:
        "{}results/.benchmarks/{{unique_rep}}.merge_peaks.benchmark.txt".format(outdir)
    conda:
        "../envs/various.yaml"
    shell:
        """
        bedtools multiinter -i {input} | awk '$4 >= 2' > {output} 2> {log}
        """


rule peakAnnot_mergedReps:
    input:
        peaks="{outdir}results/peakCalling/mergedPeaks/{{unique_rep}}_merged_optimal.bed".format(
            outdir=outdir
        ),
    output:
        annotations="{outdir}results/peakCalling/peakAnnot/{{unique_rep}}_annot.txt".format(
            outdir=outdir
        ),
        promoBed="{outdir}results/peakCalling/peakAnnot/{{unique_rep}}_promoPeaks.bed".format(
            outdir=outdir
        ),
        distalBed="{outdir}results/peakCalling/peakAnnot/{{unique_rep}}_distalPeaks.bed".format(
            outdir=outdir
        ),
        annotInfo="{outdir}results/logs/peakCalling/peakAnnot/{{unique_rep}}_annotInfo.txt".format(
            outdir=outdir
        ),
    params:
        before=config["params"]["peaksAnnotation"]["promoter"]["upstream"],
        after=config["params"]["peaksAnnotation"]["promoter"]["downstream"],
        genome=config["resources"]["ref"]["assembly"],
    log:
        "{outdir}results/logs/peakCalling/peakAnnot/{{unique_rep}}_mergedReps.log".format(
            outdir=outdir
        ),
    benchmark:
        "{}results/.benchmarks/{{unique_rep}}.mergedReps.annotatePeaks.benchmark.txt".format(
            outdir
        )
    conda:
        "../envs/Renv.yaml"
    script:
        "../scripts/PeakAnnot.R"


rule peakAnnot_singleRep:
    input:
        peaks=get_singelRep_peaks,
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


ruleorder: peakAnnot_singleRep > peakAnnot_mergedReps
