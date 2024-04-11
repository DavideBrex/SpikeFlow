
rule differential_peaks:
    input:
        rawReadsOnPeaks="{outdir}results/peakCalling/mergedPeaks/readsOnConsensusPeaks.tsv".format(
            outdir=outdir
        ),
        logFile=expand("{}results/logs/spike/{{sample}}.normFactor".format(outdir),
            sample=narrowSamples,
        ),
    output:
        diffTab="{outdir}results/differentialAnalysis/diffPeaks.tsv".format(outdir=outdir),
        volcanoPlot="{outdir}results/differentialAnalysis/volcanoPlot.pdf".format(outdir=outdir),
        pcaPlot="{outdir}results/differentialAnalysis/pcaPlot.pdf".format(outdir=outdir),
    params:
        contrast=config["diffPeakAnalysis"]["contrast"],
        padjCutoff=config["diffPeakAnalysis"]["padjust"],
        log2FCcutoff=config["diffPeakAnalysis"]["log2FCcutoff"],
        normMethod=config["normalization_type"],
    log:
        "{}results/logs/DifferentialAnalysis/DifferentialAnalysis.log".format(outdir),
    conda:
        "../envs/Renv.yaml"
    script:
        "../scripts/DifferentialPeaksAnalysis.R"