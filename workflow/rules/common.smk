import os
from snakemake.utils import validate
import pandas as pd
import numpy as np

# We start by checking the input files (samples_sheet and config.yaml) to ensure that their format is correct

samples_sheet = (
    pd.read_csv(
        config["samples_sheet"],
        dtype={"replicate": "Int64", "control_replicate": "Int64"},
        sep=",",
    )
    .set_index(["sample", "replicate"], drop=False)
    .sort_index()
)

# -------------------- validation of config and sample sheet --------------------#

validate(samples_sheet, schema="../schemas/sampleSheet.schema.yaml")

validate(config, schema="../schemas/config.schema.yaml")

# print(samples_sheet)

# -------------------- global variables defintion --------------------#

# let's get the samples that need to be merged due to presence of multiple lanes
duplicated_indices = samples_sheet.index[
    samples_sheet.index.duplicated(keep=False)
].unique()
multiLanes_samp = ["{}-rep{}".format(a, b) for a, b in duplicated_indices]

# create a dictionary of sample-input match
idSamples = samples_sheet["sample"].str.cat(
    samples_sheet["replicate"].astype(str), sep="-rep"
)
inputSamples = samples_sheet["control"].str.cat(
    samples_sheet["control_replicate"].astype(str), sep="-rep"
)

sample_to_input = dict(zip(idSamples, inputSamples))

# define narrow, broad and very-broad samples
narrowSamples = [
    "{}-rep{}".format(sample, rep)
    for sample, rep in samples_sheet[samples_sheet["peak_type"] == "narrow"]
    .index.unique()
    .tolist()
]
broadSamples = [
    "{}-rep{}".format(sample, rep)
    for sample, rep in samples_sheet[samples_sheet["peak_type"] == "broad"]
    .index.unique()
    .tolist()
]
veryBroadSamples = [
    "{}-rep{}".format(sample, rep)
    for sample, rep in samples_sheet[samples_sheet["peak_type"] == "very-broad"]
    .index.unique()
    .tolist()
]


def retrieve_index(id):
    """Function to retrieve sample and replicate from the id"""
    samp, rep = id.split("-rep")
    return (samp, int(rep))


# if there are replicates, create a dictionary to assign them to a sampleID for the different peak types
reps_dict_narrow = {}
reps_dict_broad = {}
reps_dict_verybroad = {}

if narrowSamples:
    for sample in narrowSamples:
        if retrieve_index(sample)[0] in reps_dict_narrow:
            reps_dict_narrow[retrieve_index(sample)[0]].append(sample)
        else:
            reps_dict_narrow[retrieve_index(sample)[0]] = [sample]
    reps_dict_narrow = {
        key: value for key, value in reps_dict_narrow.items() if len(value) > 1
    }  # we keep only those with more than 1 rep
if broadSamples:
    for sample in broadSamples:
        if retrieve_index(sample)[0] in reps_dict_broad:
            reps_dict_broad[retrieve_index(sample)[0]].append(sample)
        else:
            reps_dict_broad[retrieve_index(sample)[0]] = [sample]
    reps_dict_broad = {
        key: value for key, value in reps_dict_broad.items() if len(value) > 1
    }  # we keep only those with more than 1 rep
if veryBroadSamples:
    for sample in veryBroadSamples:
        if retrieve_index(sample)[0] in reps_dict_verybroad:
            reps_dict_verybroad[retrieve_index(sample)[0]].append(sample)
        else:
            reps_dict_verybroad[retrieve_index(sample)[0]] = [sample]
    reps_dict_verybroad = {
        key: value for key, value in reps_dict_verybroad.items() if len(value) > 1
    }  # we keep only those with more than 1 rep

# since the consensus peaks are divided by antibody, we need to create a dictionary with the samples for each antibody
antibody_dict = {}  
for antibody in samples_sheet["antibody"].dropna().unique():
    antibody_dict[antibody] = [
    "{}-rep{}".format(sample, rep)
    for sample, rep in samples_sheet[samples_sheet["antibody"] == antibody]
    .index.unique()
    .tolist()
    ]

# -------------------- wildcard constraints --------------------#


wildcard_constraints:
    id="|".join(set(["-rep".join(map(str, idx)) for idx in samples_sheet.index])),
    group="1|2",
    sample="|".join(
        set(
            [
                "-rep".join(map(str, idx))
                for idx in samples_sheet[
                    ~samples_sheet["antibody"].isna()
                ].index.unique()
            ]
        )
    ),
    unique_rep="|".join(
        (reps_dict_narrow | reps_dict_broad | reps_dict_verybroad).keys()
    ),


# -------------------- Sample sheet Sanity checks function ---------------#
def perform_checks(input_df):
    def check_index_files(folder_path, prefix):
        # Expected filenames
        expected_files = [
            "{}.1.bt2",
            "{}.2.bt2",
            "{}.3.bt2",
            "{}.4.bt2",
            "{}.rev.1.bt2",
            "{}.rev.2.bt2",
        ]
        # Check if the folder exists
        if not os.path.exists(folder_path):
            raise FileNotFoundError(
                "The genome index folder {} does not exist. \nPlease check that the folder is present and contains the indexing files".format(
                    folder_path
                )
            )
        # List all files in the directory to check for the presence of index files
        files_in_directory = os.listdir(folder_path)
        missing_files = []  # Check for each expected file
        for file_pattern in expected_files:
            expected_file = file_pattern.format(prefix)
            if expected_file not in files_in_directory:
                missing_files.append(expected_file)
        # Report missing files
        if missing_files:
            raise FileNotFoundError(
                """It appears that the genome index folder you provided is missing one/more indexing files.
                \nPlease check that the index prefix is correct and the index files are present in {}""".format(
                    folder_path
                )
            )

    # config file header
    header = [
        "sample",
        "replicate",
        "antibody",
        "control",
        "control_replicate",
        "peak_type",
        "fastq_1",
        "fastq_2",
    ]

    # 1. check if the header has not been changed by the user
    if list(input_df.columns) != header:
        raise ValueError("Please check samplesheet header")

    # 2. check extension of fastq files and whether the path exists
    # fastq_1
    if not all(input_df.fastq_1.map(os.path.exists)):
        raise FileNotFoundError(
            "Please check fastq_1 files paths, a file do not exist "
        )

    if not all(input_df.fastq_1.str.endswith(".fastq.gz")):
        raise ValueError("Please check fastq_1 files extension, it has to be .fastq.gz")

    # fastq_2
    if not all(input_df.fastq_2.isnull()):
        pairedEndSamp = input_df.loc[pd.notna(input_df.fastq_2), :]

        if not all(pairedEndSamp.fastq_2.map(os.path.exists)):
            raise FileNotFoundError(
                "Please check fastq_2 files paths, a file do not exist "
            )

        if not all(pairedEndSamp.fastq_2.str.endswith(".fastq.gz")):
            raise ValueError(
                "Please check fastq_2 files extension, it has to be .fastq.gz"
            )

    # 3. -check whether replicates from the same samples are all single-end or both paired-end
    #   -check if runs of the same sample   have same data type (single-end or paired -end)

    for sample in input_df.index.get_level_values("sample").unique():
        if all(input_df.loc[[sample]].fastq_2.notna()):
            pass
        elif any(input_df.loc[[sample]].fastq_2.notna()):
            raise Exception(
                "For sample {}, all replicates and runs should be either single or paired end".format(
                    sample
                )
            )
        # 4. check if all replicates have the same peak type
        if input_df.loc[[sample]].peak_type.nunique() > 1:
            raise Exception(
                "For sample {}, all replicates should have the same peak type".format(
                    sample
                )
            )
    # 5. Control identifier and replicate has to match a provided sample identifier and replicate
    input_df_controls = input_df[
        "antibody"
    ].isna()  # control samples (those with antibody to null)

    pairs_to_check = input_df[["control", "control_replicate"]]
    pairs_to_compare = input_df[["sample", "replicate"]].apply(tuple, axis=1)
    result_rows = ~pairs_to_check.apply(tuple, axis=1).isin(pairs_to_compare)

    noControl = input_df_controls ^ result_rows

    samplesNoControl = noControl[noControl == True].index.unique().tolist()
    if len(samplesNoControl) > 0:
        raise Exception(
            "The indicated control is missing in the samples column for these samples: {}".format(
                samplesNoControl
            )
        )

    # 6. in case an index is provided for the ref genome (different than ""), check whether it actually exists
    if config["resources"]["ref"]["index"] != "":
        check_index_files(
            os.path.dirname(config["resources"]["ref"]["index"]),
            os.path.basename(config["resources"]["ref"]["index"]),
        )

    # 7. check if the chromsome sizes file exists and if the blacklist file exists
    if not os.path.exists(config["params"]["peakCalling"]["chrom_sizes"]):
        raise FileNotFoundError(
            "The provided path to the chromosome sizes file does not exist. \nPlease check that the file is present"
        )
    if not os.path.exists(config["resources"]["ref"]["blacklist"]):
        raise FileNotFoundError(
            "The provided path to the blacklist file does not exist. \nPlease check that the file is present"
        )

    # 8. check if the contrast is correctly defined
    if config["diffPeakAnalysis"]["perform_diff_analysis"]:
        contrastsToCheck = config["diffPeakAnalysis"]["contrasts"]
        # diff bind analysis is performed only on the samples with antibody (and with same antibody value)
        for antibodyItem in contrastsToCheck:
            # we check that the antibody is present in the samplesheet
            if antibodyItem not in input_df[['antibody']].values:
                raise ValueError(
                    "Please indicate a valid antibody in the contrasts (config file) for the differential binding analysis\n"+
                    "The antibody has to be defined in the samplesheet for each sampel (not input)"
                )
            subdf = input_df[input_df['antibody'] == antibodyItem]
            sample_groups_antibody =  [sample.rsplit("_",1)[1] for sample in subdf.index.get_level_values("sample").unique().tolist()]
            #we check that the conditions in the contrast are present in the samplesheet and with right format
            for contrast in contrastsToCheck[antibodyItem]:
                if "_vs_" in contrast:
                    contrastElem = contrast.split("_vs_")
                    if len(contrastElem) != 2:
                        raise ValueError(
                            "The contrast should be defined as 'groupA_vs_groupB'. Please check the contrast definition in the config file\n"+
                            "Group has to be defined as 'sampleName_groupA' in the sample column of sample_sheet.csv"
                        )
                    if contrastElem[0] not in sample_groups_antibody or contrastElem[1] not in sample_groups_antibody:
                        raise ValueError(
                            "One of the group in the contrast is not present in the samplesheet. Please check the contrast definition in the config file\n"+
                            "Group has to be defined as 'sampleName_groupA' in the sample column of sample_sheet.csv"

                        )
                    if contrastElem[0] == contrastElem[1]:
                        raise ValueError(
                            "The groups in the contrast should be different. Please check the contrast definition in the config file\n"+
                            "Group has to be defined as 'sampleName_groupA' in the sample column of sample_sheet.csv"
                        )
                else:
                    raise ValueError(
                        "The contrast should be defined as 'groupA_vs_groupB'. Please check the contrast definition in the config file\n"+
                        "Group has to be defined as 'sampleName_groupA' in the sample column of sample_sheet.csv"
                    )
            # 9. if contrast is okay, we need to check the peak type are the same for all samples
            if len(subdf.peak_type.unique()) > 1:
                raise ValueError("The peak type is not the same for all samples with antibody {}. For differential peaks please specify same peak type".format(antibodyItem))
            if subdf.peak_type.unique() == "very-broad":
                raise ValueError(
                    "The differential binding analysis can not be performed on very-broad peaks.\n"+
                    "Consider to change peak type or disable the differential binding analysis in the config file!"
                )
                

# -------------------- Sample sheet Sanity checks ---------------#

perform_checks(samples_sheet)

# -------------------- Define input files for rule all ---------------#


def input_toget():
    # bigwigs
    bigWigs = expand(
        "{}results/bigWigs/{{id}}.bw".format(outdir), id=sample_to_input.keys()
    )

    # qc and peak calling
    QCfiles = [
        "{}results/QC/multiqc/SpikeFlow_multiqc_report.html".format(outdir),
        "{}results/QC/peaks_annotation_mqc.tsv".format(outdir),
    ]
    peak_files = []
    annot_files = []

    if narrowSamples:
        QCfiles.append("{}results/QC/macs2_peaks_mqc.tsv".format(outdir))
        peak_files += expand(
            "{}results/peakCalling/macs2/{{sample}}_peaks.narrowPeak".format(
                outdir
            ),
            sample=narrowSamples,
        )
        #for testing Spiker peak calling, we also add those output files
        peak_files += expand(
            "{}results/peakCallingNorm/{{sample}}.treat.pileup.SpikeIn_scaled.bdg".format(
                outdir
            ),
            sample=narrowSamples,
        )
        annot_files += expand(
            "{}results/peakCalling/peakAnnot/{{sample}}_annot.txt".format(outdir),
            sample=narrowSamples,
        )
    if broadSamples:
        QCfiles.append("{}results/QC/epic2_peaks_mqc.tsv".format(outdir))
        peak_files += expand(
            "{}results/peakCalling/epic2/{{sample}}_broadPeaks.bed".format(outdir),
            sample=broadSamples,
        )
        #for testing Spiker peak calling, we also add those output files
        peak_files += expand(
            "{}results/peakCallingNorm/{{sample}}.treat.pileup.SpikeIn_scaled.bdg".format(
                outdir
            ),
            sample=broadSamples,
        )
        annot_files += expand(
            "{}results/peakCalling/peakAnnot/{{sample}}_annot.txt".format(outdir),
            sample=broadSamples,
        )
    if veryBroadSamples:
        QCfiles.append("{}results/QC/edd_peaks_mqc.tsv".format(outdir))
        peak_files += expand(
            "{}results/peakCalling/edd/{{sample}}/{{sample}}_peaks.bed".format(outdir),
            sample=veryBroadSamples,
        )

    # if there are reps, we need to add to the peak files the merged samples
    if narrowSamples or broadSamples:
        merged_dicts = reps_dict_narrow | reps_dict_broad  # merge dicts
        if merged_dicts:
            peak_files += expand(
                "{}results/peakCalling/mergedPeaks/{{unique_rep}}_merged_optimal.bed".format(
                    outdir
                ),
                unique_rep=list(merged_dicts.keys()),
            )
            annot_files += expand(
                "{}results/peakCalling/peakAnnot/{{unique_rep}}_annot.txt".format(
                    outdir
                ),
                unique_rep=list(merged_dicts.keys()),
            )
    # for very broad peaks we shall use only intersect, whereas for narrow and broad we shall use chip-r
    if veryBroadSamples:
        if reps_dict_verybroad:
            peak_files += expand(
                "{}results/peakCalling/mergedPeaks/{{unique_rep}}_merged_intersected.bed".format(
                    outdir
                ),
                unique_rep=list(reps_dict_verybroad.keys()),
            )

    # we need the consensus peaks for the differential analysis
    if config["diffPeakAnalysis"]["perform_diff_analysis"]:
        diff_peak_files = []
        # we add to otputs the different combinations of antibody and contrast
        for antibody, contrasts in config["diffPeakAnalysis"]["contrasts"].items():
            for contrast in contrasts:
                path = "{outdir}results/differentialAnalysis/{antibody}/{antibody}_{contrast}_diffPeaks.tsv".format(
                    outdir=outdir, antibody=antibody, contrast=contrast
                )
                diff_peak_files.append(path)

        return bigWigs + peak_files + QCfiles + annot_files + diff_peak_files  
    else:
        return bigWigs + peak_files + QCfiles + annot_files


# -------------------- Other helpers functions ---------------#


def is_single_end(id):
    samp, rep = retrieve_index(id)
    check = pd.isnull(samples_sheet.loc[(samp, rep), "fastq_2"])
    # in case a sample has multiple lanes, we get a series instead of str
    if isinstance(check, pd.Series):
        return check.iloc[0]
    return check


# --------------------  Rules Input Functions ---------------#


def get_fastq(wildcards):
    """Function called by merged lanes. It is executed only when a sample has multiple lanes only"""
    samp, rep = retrieve_index(**wildcards)

    if is_single_end(**wildcards):
        return (
            samples_sheet.loc[(samp, rep), "fastq_1"]
            if wildcards.id in multiLanes_samp
            else []
        )
    else:
        u = samples_sheet.loc[(samp, rep), ["fastq_1", "fastq_2"]].dropna()
        return (
            {"fw": u.fastq_1.tolist(), "rv": u.fastq_2.tolist()}
            if wildcards.id in multiLanes_samp
            else {"fw": "", "rv": ""}
        )


def get_fastq_trimming(wildcards):
    """Function called by fastp_pe or se. Only called when trimming is activated"""

    samp, rep = retrieve_index(**wildcards)
    if is_single_end(**wildcards):
        # to run merge only on samples that have multiple lanes
        if wildcards.id in multiLanes_samp:
            return expand("{}results/fastq/{id}.fastq.gz".format(outdir, **wildcards))
        else:
            toret = samples_sheet.loc[
                (samp, rep), "fastq_1"
            ]  # we need this check beacuse if multindex has duplicated, loc returns a series not a str
            return toret.tolist() if isinstance(toret, pd.Series) else [toret]
    else:
        if wildcards.id in multiLanes_samp:
            return expand(
                "{}results/fastq/{{id}}_{{group}}.fastq.gz".format(outdir),
                group=[1, 2],
                **wildcards,
            )
        else:
            u = samples_sheet.loc[
                (samp, rep), ["fastq_1", "fastq_2"]
            ].dropna()  # we need this check because if multindex has duplicated, loc returns a df not a series
            return (
                [u.fastq_1, u.fastq_2]
                if isinstance(u, pd.Series)
                else [u.fastq_1.tolist()[0], u.fastq_2.tolist()[0]]
            )


def get_reads(wildcards):
    """Function called by aligners."""

    samp, rep = retrieve_index(**wildcards)
    # if trimming is performed, the trimmed fastqs are all in trimmed folder
    if config["trimming"]:
        if is_single_end(**wildcards):
            return expand("{}results/trimmed/{id}.fastq.gz".format(outdir, **wildcards))
        else:
            return expand(
                "{}results/trimmed/{{id}}_{{group}}.fastq.gz".format(outdir),
                group=[1, 2],
                **wildcards,
            )

    else:
        if is_single_end(**wildcards):
            # to run merge only on samples that have multiple lanes
            if wildcards.id in multiLanes_samp:
                return expand(
                    "{}results/fastq/{id}.fastq.gz".format(outdir, **wildcards)
                )
            else:
                toret = samples_sheet.loc[(samp, rep), "fastq_1"]
                return toret.tolist() if isinstance(toret, pd.Series) else [toret]
        else:
            if wildcards.id in multiLanes_samp:
                return expand(
                    "{}results/fastq/{{id}}_{{group}}.fastq.gz".format(outdir),
                    group=[1, 2],
                    **wildcards,
                )
            else:
                u = samples_sheet.loc[
                    (samp, rep), ["fastq_1", "fastq_2"]
                ].dropna()  # we need this check because if multindex has duplicated, loc returns a df not a series
                return (
                    [u.fastq_1, u.fastq_2]
                    if isinstance(u, pd.Series)
                    else [u.fastq_1.tolist()[0], u.fastq_2.tolist()[0]]
                )


def get_replicate_peaks(wildcards):
    """Function that returns the input files to merge replicates (if any) both narrow and broad peaks"""
    if wildcards.unique_rep in reps_dict_narrow:
        return expand(
            "{}results/peakCalling/macs2/{{sample}}_peaks.narrowPeak".format(
                outdir
            ),
            sample=reps_dict_narrow[wildcards.unique_rep],
        )
    elif wildcards.unique_rep in reps_dict_broad:
        return expand(
            "{}results/logs/peakCalling/epic2/{{sample}}.bed".format(outdir),
            sample=reps_dict_broad[wildcards.unique_rep],
        )


def get_replicate_peaks_edd(wildcards):
    """Function that returns the input files to merge replicates (if any) for very broad  peaks"""
    return expand(
        "{}results/peakCalling/edd/{{sample}}/{{sample}}_peaks.bed".format(outdir),
        sample=reps_dict_verybroad[wildcards.unique_rep],
    )


def get_singelRep_peaks(wildcards):
    """Function that returns the input files to annot single sample peak files both narrow and broad peaks"""
    if wildcards.sample in narrowSamples:
        return "{}results/peakCalling/macs2/{{sample}}_peaks.narrowPeak".format(
            outdir, sample=wildcards.sample
        )
    elif wildcards.sample in broadSamples:
        return "{}results/peakCalling/epic2/{{sample}}_broadPeaks.bed".format(
            outdir, sample=wildcards.sample
        )


def get_bams_by_antibody(wildcards):
    """Function that returns the bam files for the antibody for the consensus peaks"""
    return expand(
        "{}results/bam/{{sample}}_ref.sorted.bam".format(outdir),
        sample=antibody_dict[wildcards.antibody],
    )

def get_normFactor_by_antibody(wildcards):
    """Function that returns the normalization factors for the antibody for the diff peaks analysis"""
    return expand(
        "{}results/logs/spike/{{sample}}.normFactor".format(outdir),
        sample=antibody_dict[wildcards.antibody],
    )
    
def get_contrasts_by_antibody(wildcards):
    """Function that returns the contrasts for the antibody for the diff peaks analysis"""
    return expand(
            "{}results/differentialAnalysis/{{antibody}}/{{contrast}}_diffPeaks.tsv".format(outdir),
            antibody=wildcards.antibody,
            contrast=config["diffPeakAnalysis"]["contrasts"][wildcards.antibody],  
        )

# --------------------  Rules Functions ---------------#
def normalization_factor(wildcards):
    """
    Read the log message from the cleaning of bam files and compute normalization factor
    By using the log file we avoid to read the bam file with pysam just to get # of aligned reads
    """
    norm_type = config["normalization_type"]
    # open sample log file
    samp = wildcards.id
    with open(
        "{}results/logs/spike/{}.removeSpikeDups".format(outdir, samp), "r"
    ) as file:
        info_sample = file.read().strip().split("\n")

        Nsample = int(
            info_sample[1].split(":")[-1]
        )  # number of aligned reads in sample
        Nspike = int(info_sample[2].split(":")[-1])  # number of spike reads in sample

        # we need the information also from the input (if it is not the sample an input itself)
        inputSamp = sample_to_input[wildcards.id]
        if pd.isna(inputSamp) or norm_type != "RX-Input":
            if norm_type == "Orlando":
                alpha = (1 / Nspike) * 1000000  # From Orlando et. al 2014 (RRPM)
            else:
                alpha = (1 / Nsample) * 1000000  # RPM
        else:
            # open input log file
            with open(
                "{}results/logs/spike/{}.removeSpikeDups".format(outdir, inputSamp), "r"
            ) as input_file:
                info_input = input_file.read().strip().split("\n")

                gamma = int(info_input[2].split(":")[-1]) / int(
                    info_input[1].split(":")[-1]
                )  # ratio spike/samples in input
                alpha = gamma / Nspike * 1000000  # normalization factor

        # TO DO: add log file with the norm factors stored
    with open("{}results/logs/spike/{}.normFactor".format(outdir, samp), "w") as file:
        file.write("Normalization factor:{}\n".format(round(alpha, 4)))

    if is_single_end(wildcards.id):
        return "--scaleFactor {} --extendReads {}".format(
            str(round(alpha, 4)), str(config["params"]["deeptools"]["read_extension"])
        )
    else:
        return "--scaleFactor {} --extendReads ".format(str(round(alpha, 4)))

def spiker_normalization_factor(wildcards):
    """
    Function called by Spiker peak calling
    It returns the normalization factors for treatment and control samples once calculated by the function above
    """
    treatment_file = "{}results/logs/spike/{}.normFactor".format(outdir, wildcards.sample)
    #since the peak calling is done only on treatment samples, we need to get the input sample (and it has to have a control)
    control_file = "{}results/logs/spike/{}.normFactor".format(outdir, sample_to_input[wildcards.sample])
    with open(treatment_file) as tf, open(control_file) as cf:
        treatment_norm_factor = tf.read().strip().split(":")[-1].strip()
        control_norm_factor = cf.read().strip().split(":")[-1].strip()
    return "--csf {} --tsf {}".format(control_norm_factor, treatment_norm_factor)


def check_peak_type(wildcards):
    """Function to check the peak type of the sample"""
    if wildcards.sample in narrowSamples:
        return ""
    elif wildcards.sample in broadSamples:
        return "--broad"