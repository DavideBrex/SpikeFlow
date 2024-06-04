
rule differential_peaks:
    input:
        rawReadsOnPeaks="{outdir}results/peakCallingNorm/mergedPeaks/{{antibody}}_readsOnConsensusPeaks.tsv".format(
            outdir=outdir
        )
        if config["diffPeakAnalysis"]["useSpikeinCalledPeaks"]
        else "{outdir}results/peakCalling/mergedPeaks/{{antibody}}_readsOnConsensusPeaks.tsv".format(
            outdir=outdir
        ),
        logFile=get_normFactor_by_antibody,
    output:
        diffTab="{outdir}results/differentialAnalysis/NormalisedPeaks/{{antibody}}/{{antibody}}_{{contrast}}_diffPeaks.tsv".format(
            outdir=outdir
        )
        if config["diffPeakAnalysis"]["useSpikeinCalledPeaks"]
        else "{outdir}results/differentialAnalysis/{{antibody}}/{{antibody}}_{{contrast}}_diffPeaks.tsv".format(
            outdir=outdir
        ),
    params:
        padjCutoff=config["diffPeakAnalysis"]["padjust"],
        log2FCcutoff=config["diffPeakAnalysis"]["log2FCcutoff"],
        normMethod=config["normalization_type"],
        outdir=lambda w, output: os.path.dirname(output["diffTab"]) + "/",
    log:
        "{}results/logs/DifferentialAnalysis/{{antibody}}_{{contrast}}_DifferentialAnalysis.log".format(
            outdir
        ),
    conda:
        "../envs/Renv.yaml"
    script:
        "../scripts/DifferentialPeaksAnalysis.R"
