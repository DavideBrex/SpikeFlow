
import os
from snakemake.utils import validate
import pandas as pd
import numpy as np

#we start by checking the input files (samples_sheet and config.yaml) to ensure that their format is correct

samples_sheet = pd.read_csv(config["samples_sheet"], dtype={
    "replicate": "Int64",
    "control_replicate": "Int64",
    'spike': "boolean"}, sep=",").set_index(["sample","replicate"], drop=False).sort_index()


#-------------------- validation of config and sample sheet --------------------#

validate(samples_sheet, schema="../schemas/sampleSheet.schema.yaml")

validate(config, schema="../schemas/config.schema.yaml")

#print(samples_sheet)

#let's get the samples that need to be merged due to presence of multiple lanes
duplicated_indices = samples_sheet.index[samples_sheet.index.duplicated(keep=False)].unique()
multiLanes_samp = [f"{a}-rep{b}" for a, b in duplicated_indices] 

#-------------------- wildcard constraints --------------------#

wildcard_constraints:
    id = "|".join(set(['-rep'.join(map(str, idx)) for idx in samples_sheet.index])),
    group = "1|2"

#-------------------- Sample sheet Sanity checks function ---------------#
def perform_checks(input_df):
    

    header=["sample","replicate","antibody","control","control_replicate","spike","peak_type","fastq_1","fastq_2"]

    #1. check if the header has not been changed by the user
    if list(input_df.columns) != header:
        raise ValueError("Please check samplesheet header")
    
    #2. check extension of fastq files and whether the path exists
    #fastq_1
    if not all(input_df.fastq_1.map(os.path.exists)):
        raise FileNotFoundError("Please check fastq_1 files paths, a file do not exist ")

    if not all(input_df.fastq_1.str.endswith(".fastq.gz")):
        raise ValueError("Please check fastq_1 files extension, it has to be .fastq.gz")

    #fastq_2
    if not all(input_df.fastq_2.isnull()):

        pairedEndSamp=input_df.loc[ pd.notna(input_df.fastq_2), :]

        if not all(pairedEndSamp.fastq_2.map(os.path.exists)):
            raise FileNotFoundError("Please check fastq_2 files paths, a file do not exist ")

        if not all(pairedEndSamp.fastq_2.str.endswith(".fastq.gz")):
            raise ValueError("Please check fastq_2 files extension, it has to be .fastq.gz")

    #3. -check whether replicates from the same samples are all single-end or both paired-end
    #   -check if runs of the same sample   have same data type (single-end or paired -end)
    #   -and also that spike column has the same values for all the reps of a sample
    
    for sample in input_df.index.get_level_values('sample').unique():
        if all(input_df.loc[[sample]].fastq_2.notna()):
            if all(input_df.loc[[sample]].spike == True) or all(input_df.loc[[sample]].spike == False):
                pass
            else:
                raise ValueError("For sample {}, all replicates should have the same value for spike column".format(sample))
            pass 
        elif any(input_df.loc[[sample]].fastq_2.notna()):
            raise Exception("For sample {}, all replicates and runs should be either single or paired end".format(sample))


    #4. Control identifier and replicate has to match a provided sample identifier and replicate
    input_df_controls = input_df['antibody'].isna() #control sames (those with antibody to null)

    pairs_to_check = input_df[['control', 'control_replicate']]
    pairs_to_compare = input_df[['sample', 'replicate']].apply(tuple, axis=1)
    result_rows = ~pairs_to_check.apply(tuple, axis=1).isin(pairs_to_compare)

    noControl=input_df_controls ^ result_rows

    samplesNoControl=noControl[noControl == True].index.unique().tolist()
    if len(samplesNoControl) > 0:
        raise Exception("ERROR: The indicated control is missing in the samples column for these samples: {}".format(samplesNoControl))
    
    #5. in case an index is provided for the ref genome (different than ""), check whether it actually exists
    if config["resources"]["ref"]["index"] != "":
        if not os.path.exists(os.path.dirname(config["resources"]["ref"]["index"])):
            raise FileNotFoundError("The provided path to the reference genome index does not exist. \nPlease check that the folder is present and contains the indexing files")
    #same for spike
    if config["resources"]["ref_spike"]["index_spike"] != "":
        if not os.path.exists(os.path.dirname(config["resources"]["ref_spike"]["index_spike"])):
            raise FileNotFoundError("The provided path to the spike genome index does not exist. \nPlease check that the folder is present and contains the indexing files")
        

#-------------------- Sample sheet Sanity checks ---------------#

perform_checks(samples_sheet)

#-------------------- Define input files for rule all ---------------#


#     return wanted_inputs
def input_toget():

    wanted_inputs=[]
    for (sample, replicate) in samples_sheet.index.unique():

        wanted_inputs += [f"{sample}-rep{replicate}"]


    return  expand("results/bam/{id}.bam", id=wanted_inputs)


#-------------------- Other helpers functions ---------------#

def retrieve_index(id):
    samp, rep = id.split("-rep")
    return (samp, int(rep))


def is_single_end(id):
    samp, rep = retrieve_index(id)
    check = pd.isnull(samples_sheet.loc[(samp, rep), "fastq_2"])
    #in case a sample has multiple lanes, we get a series instead of str
    if isinstance(check, pd.Series):
        return check[0]
    return check

def is_spike(id):
    samp, rep = retrieve_index(id)
    check = samples_sheet.loc[(samp, rep), "spike"]
    #in case a sample has multiple lanes, we get a series instead of str
    if isinstance(check, pd.Series):
        return check[0]
    return check
#--------------------  Rules Input Functions ---------------#


def get_fastq(wildcards):
    """  Function called by merged lanes. It is executed only when a sample has multiple lanes only """
    samp, rep = retrieve_index(**wildcards)

    if is_single_end(**wildcards):
        return samples_sheet.loc[(samp, rep), "fastq_1"] if wildcards.id in multiLanes_samp else  []
    else:
        u = samples_sheet.loc[ (samp, rep), ["fastq_1", "fastq_2"] ].dropna()
        return  {"fw": u.fastq_1.tolist(), "rv": u.fastq_2.tolist()}  if wildcards.id in multiLanes_samp else {"fw": "", "rv": ""} 
        

def get_fastq_trimming(wildcards):
    """  Function called by fastp_pe or se. Only called when trimming is activated """

    samp, rep = retrieve_index(**wildcards)

    if is_single_end(**wildcards):
        # to run merge only on samples that have multiple lanes
        if wildcards.id in multiLanes_samp: 
            return expand("results/fastq/{id}.fastq.gz".format(**wildcards))
        else:
            return samples_sheet.loc[(samp, rep), "fastq_1"]
    else:
        if wildcards.id in multiLanes_samp: 
            return expand("results/fastq/{id}_{group}.fastq.gz", group=[1, 2], **wildcards)
        else:
            u = samples_sheet.loc[ (samp, rep), ["fastq_1", "fastq_2"] ].dropna()
            return [ u.fastq_1.tolist()[0], u.fastq_2.tolist()[0] ]


def get_reads(wildcards):
    """  Function called by aligners. """

    samp, rep = retrieve_index(**wildcards)
    #if trimming is performed, the trimmed fastqs are all in 
    if config["trimming"]:
        if is_single_end(**wildcards):
            return expand("results/trimmed/{id}.fastq.gz".format(**wildcards))
        else:
            return expand("results/trimmed/{id}_{group}.fastq.gz", group=[1, 2], **wildcards)

    else:
        if is_single_end(**wildcards):
            # to run merge only on samples that have multiple lanes
            if wildcards.id in multiLanes_samp: 
                return expand("results/fastq/{id}.fastq.gz".format(**wildcards))
            else:
                return samples_sheet.loc[(samp, rep), "fastq_1"]
        else:
            if wildcards.id in multiLanes_samp: 
                return expand("results/fastq/{id}_{group}.fastq.gz", group=[1, 2], **wildcards)
            else:
                u = samples_sheet.loc[ (samp, rep), ["fastq_1", "fastq_2"] ].dropna()
                return [ u.fastq_1.tolist()[0], u.fastq_2.tolist()[0] ]


def get_reads_spike(wildcards):
    """  Function called by aligners. """

    samp, rep = retrieve_index(**wildcards)
    if (is_spike(**wildcards)):
        #if trimming is performed, the trimmed fastqs are all in 
        if config["trimming"]:
            if is_single_end(**wildcards):
                return expand("results/trimmed/{id}.fastq.gz".format(**wildcards))
            else:
                return expand("results/trimmed/{id}_{group}.fastq.gz", group=[1, 2], **wildcards)

        else:
            if is_single_end(**wildcards):
                # to run merge only on samples that have multiple lanes
                if wildcards.id in multiLanes_samp: 
                    return expand("results/fastq/{id}.fastq.gz".format(**wildcards))
                else:
                    return samples_sheet.loc[(samp, rep), "fastq_1"]
            else:
                if wildcards.id in multiLanes_samp: 
                    return expand("results/fastq/{id}_{group}.fastq.gz", group=[1, 2], **wildcards)
                else:
                    u = samples_sheet.loc[ (samp, rep), ["fastq_1", "fastq_2"] ].dropna()
                    return [ u.fastq_1.tolist()[0], u.fastq_2.tolist()[0] ]

