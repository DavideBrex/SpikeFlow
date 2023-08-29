
import os
from snakemake.utils import validate
import pandas as pd
import numpy as np
#TO DO: 
#  - add validation schema for config file 
#  - add check that all lanes have same settings
#
#

#we start by checking the input files (samples_sheet and config.yaml) to ensure that their format is correct

samples_sheet = pd.read_csv(config["samples_sheet"], dtype={
    "replicate": "Int64",
    "control_replicate": "Int64",
    'spike': "boolean"}, sep=",").set_index(["sample","replicate"], drop=False).sort_index()


validate(samples_sheet, schema="../schemas/sampleSheet.schema.yaml")


#let's get the samples that need to be merged due to presence of multiple lanes
duplicated_indices = samples_sheet.index[samples_sheet.index.duplicated(keep=False)].unique()
multiLanes_samp = [f"{a}-rep{b}" for a, b in duplicated_indices] 


##### wildcard constraints #####

wildcard_constraints:
    id = "|".join(set(['-rep'.join(map(str, idx)) for idx in samples_sheet.index])),
    group = "1|2"



#-------------------- Sample sheet Sanity checks    function ---------------#
def perform_checks(input_df):
    

    header=["sample","replicate","antibody","control","control_replicate","spike","peak_type","fastq_1","fastq_2"]

    #1. check if the header has not been changed by the user
    if list(input_df.columns) != header:
        print("ERROR: Please check samplesheet header")
        sys.exit(1)


    #2. check extension of fastq files and whether the path exists

    #fastq_1
    if not all(input_df.fastq_1.map(os.path.exists)):
        print("ERROR: Please check fastq_1 files paths, a file do not exist ")
        sys.exit(1)

    if not all(input_df.fastq_1.str.endswith(".fastq.gz")):
        print("ERROR: Please check fastq_1 files extension, it has to be .fastq.gz")
        sys.exit(1)

    #fastq_2
    if not all(input_df.fastq_2.isnull()):

        pairedEndSamp=input_df.loc[ pd.notna(input_df.fastq_2), :]
        #print(pairedEndSamp)

        if not all(pairedEndSamp.fastq_2.map(os.path.exists)):
            print("ERROR: Please check fastq_2 files paths, a file do not exist ")
            sys.exit(1)

        if not all(pairedEndSamp.fastq_2.str.endswith(".fastq.gz")):
            print("ERROR: Please check fastq_2 files extension, it has to be .fastq.gz")
            sys.exit(1)


    #3. check whether replicates from the same samples are all single-end or both paired-end
    #   also other runs of the same samples  must have same data type (single-end or paired -end)
    
    for sample in input_df.index.get_level_values('sample').unique():
        #print(input_df.loc[[sample]])
        if all(input_df.loc[[sample]].fastq_2.notna()):
            continue
        elif any(input_df.loc[[sample]].fastq_2.notna()):
            print("ERROR: for sample {}, all replicates and runs should be either single or paired end".format(sample))
            sys.exit(1)

    #4. Control identifier and replicate has to match a provided sample identifier and replicate
    input_df_controls = input_df['antibody'].isna() #control sames (those with antytbody to null)

    pairs_to_check = input_df[['control', 'control_replicate']]
    pairs_to_compare = input_df[['sample', 'replicate']].apply(tuple, axis=1)
    result_rows = ~pairs_to_check.apply(tuple, axis=1).isin(pairs_to_compare)

    noControl=input_df_controls ^ result_rows

    samplesNoControl=noControl[noControl == True].index.unique().tolist()
    if len(samplesNoControl) > 0:
            print("ERROR: The indicated control is missing in the samples column for these samples: {}".format(samplesNoControl))
            sys.exit(1)

#-------------------- Sample sheet Sanity checks ---------------#

perform_checks(samples_sheet)

#-------------------- Define input files for rule all ---------------#

# def input_toget():

#     wanted_inputs=[]
#     for (sample, replicate) in samples_sheet.index:

#         wanted_inputs.extend(
#             expand(
#                 [
#                     "results/trimmed/{sample}.{replicate}.fastq.gz"
#                 ],
#                 sample = sample,
#                 replicate = replicate
#             )
#         )


#     return wanted_inputs

def input_toget():

    wanted_inputs=[]
    for (sample, replicate) in samples_sheet.index.unique():

        wanted_inputs += [f"{sample}-rep{replicate}"]


    return  expand("results/bam/{id}.bam", id=wanted_inputs)


#-------------------- Other useful functions ---------------#

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

# get input files -- helper functions

def get_fastq(wildcards):
    """  Rule called by merged lanes. It is executed when a sample has multiple lanes only """
    samp, rep = retrieve_index(**wildcards)

    if is_single_end(**wildcards):
        return samples_sheet.loc[(samp, rep), "fastq_1"] if wildcards.id in multiLanes_samp else  []
    else:
        u = samples_sheet.loc[ (samp, rep), ["fastq_1", "fastq_2"] ].dropna()
        return  {"fw": u.fastq_1.tolist(), "rv": u.fastq_2.tolist()} if wildcards.id in multiLanes_samp else  []
        

def get_fastq_trimming(wildcards):
    """  Rule called by fastp_pe or se. Only called when trimming is activated """
    samp, rep = retrieve_index(**wildcards)

    if is_single_end(**wildcards):
        # to run merge only on samples that have multiple lanes
        if wildcards.id in multiLanes_samp: 
            return "results/fastq/{id}.fastq.gz".format(**wildcards)
        else:
            return samples_sheet.loc[(samp, rep), "fastq_1"]
    else:
        if wildcards.id in multiLanes_samp: 
            return expand("results/fastq/{id}_{group}.fastq.gz", group=[1, 2], **wildcards)
        else:
            u = samples_sheet.loc[ (samp, rep), ["fastq_1", "fastq_2"] ].dropna()
            return [ u.fastq_1.tolist()[0], u.fastq_2.tolist()[0] ]


def get_reads(wildcards):
    """  Rule called by aligners. """

    samp, rep = retrieve_index(**wildcards)
    print(wildcards.id)
    #if trimming is performed, the trimmed fastqs are all in 
    if config["trimming"]:
        if is_single_end(**wildcards):
            return "results/trimmed/{id}.fastq.gz".format(**wildcards)
        else:
            return expand("results/trimmed/{id}_{group}.fastq.gz", group=[1, 2], **wildcards)

    else:
        if is_single_end(**wildcards):
            # to run merge only on samples that have multiple lanes
            if wildcards.id in multiLanes_samp: 
                return "results/fastq/{id}.fastq.gz".format(**wildcards)
            else:
                return samples_sheet.loc[(samp, rep), "fastq_1"]
        else:
            if wildcards.id in multiLanes_samp: 
                return expand("results/fastq/{id}_{group}.fastq.gz", group=[1, 2], **wildcards)
            else:
                u = samples_sheet.loc[ (samp, rep), ["fastq_1", "fastq_2"] ].dropna()
                return [ u.fastq_1.tolist()[0], u.fastq_2.tolist()[0] ]


