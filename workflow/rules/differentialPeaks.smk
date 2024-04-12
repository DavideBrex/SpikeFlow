
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
    params:
        contrast=config["diffPeakAnalysis"]["contrast"],
        padjCutoff=config["diffPeakAnalysis"]["padjust"],
        log2FCcutoff=config["diffPeakAnalysis"]["log2FCcutoff"],
        normMethod=config["normalization_type"],
        outdir='{outdir}results/differentialAnalysis/'.format(outdir=outdir),
    log:
        "{}results/logs/DifferentialAnalysis/DifferentialAnalysis.log".format(outdir),
    conda:
        "../envs/Renv.yaml"
    script:
        "../scripts/DifferentialPeaksAnalysis.R"