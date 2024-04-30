import os
import pandas as pd
import sys

sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout
#read input files and config
logFiles=snakemake.input.logFiles

print("Input files: {}".format(logFiles))

normType=snakemake.params.typeOfNorm

dictToMatchInput=snakemake.params.sampleToInput
dictAntibody=snakemake.params.antibody_dict

print(dictAntibody)

outdir=snakemake.params.outdir
#For downsampling or median normalization we need to compute the scaling factors
#together for the samples from the same antibody
#------------------------------------------------#
#------------------------------------------------#
#For Orlando, RAW and RX-Input we calculate single-sample norm facotors
#------------------------------------------------#
#------------------------------------------------#
if normType in ['RAW', 'Orlando', 'RX-Input']:
    print("Computing single-sample normalization factors")
    for idSample in logFiles:
        idName = os.path.basename(idSample).rsplit(".",1)[0]
        with open(idSample, "r") as file:
            info_sample = file.read().strip().split("\n")

            Nsample = int(
                info_sample[0].split(":")[-1]
            )  # number of aligned reads in sample
            Nspike = int(info_sample[1].split(":")[-1])  # number of spike reads in sample
            print("Sample: {} - Nsample: {} - Nspike: {}".format(idSample, Nsample, Nspike))

            # we need the information also from the input (if it is not the sample an input itself)
            inputSamp = dictToMatchInput[idName]
            if pd.isna(inputSamp) or normType != "RX-Input":
                if normType == "Orlando":
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

        # store normalization factor
        with open(idSample.replace('.removeSpikeDups','.normFactor'), "w") as file:
            file.write("Normalization factor:{}\n".format(round(alpha, 4)))

#------------------------------------------------#
#------------------------------------------------#
#For Downsampling and median normalization we calculate the scaling factors booling together the samples from the same antibody        
elif normType in ['Downsampling', 'Median']:
    print("Computing pooled scaling factors for each antibody")


    #we check whether all samples in dictAntibody are present in logFiles
    #if not we raise an error
    logFilesNames = [os.path.basename(x).rsplit(".",1)[0] for x in logFiles]
    for antibody in dictAntibody:
        samples = dictAntibody[antibody]
        for sampleName in samples:
            if sampleName not in logFilesNames:
                raise ValueError("Sample {} not found in logFiles".format(sampleName))

    #we calculate the normalization factors for each antibody
    sampleInfo = {}
    for idSample in logFiles:
        idName = os.path.basename(idSample).rsplit(".",1)[0]
        with open(idSample, "r") as file:
            info_sample = file.read().strip().split("\n")
            Nsample = int(info_sample[0].split(":")[-1])
            Nspike = int(info_sample[1].split(":")[-1])
            sampleInfo[idName] = [Nsample, Nspike]

    print(sampleInfo)
    perAntibodyScale = {}
    #now we calculate the scaling factors for each antibody
    for antibody in dictAntibody:
        print("Computing scaling factors for antibody: {}\n".format(antibody))
        samples = dictAntibody[antibody]
        
        spikeInReads = [elem[1] for elem in [sampleInfo[sampleName] for sampleName in samples]]
        if normType == "Median":
            commonNorm = sum(spikeInReads) / len(spikeInReads)
            print("Median spike reads: {}".format(commonNorm))
        if normType == "Downsampling":
            commonNorm = min(spikeInReads)
            print("Minimum spike reads: {}".format(commonNorm))

        perAntibodyScale[antibody] = commonNorm

    #now we calculate per sample normalization factors
    for idSample in logFiles:
        
        idName = os.path.basename(idSample).rsplit(".",1)[0]
        
        inputSamp = dictToMatchInput[idName]
        if pd.isna(inputSamp):
            #if the sample is an input, we use RAW normalization
            Nsample = sampleInfo[idName][0]
            alpha = (1 / Nsample) * 1000000  # RPM
        else:
            antibody = [key for key, value in dictAntibody.items() if idName in value][0]
            commonNormFactor = perAntibodyScale[antibody]
            Nspike = sampleInfo[idName][1]            
            alpha = (commonNorm / Nspike) 
        
        with open(idSample.replace('.removeSpikeDups','.normFactor'), "w") as file:
            file.write("Normalization factor:{}\n".format(round(alpha, 4)))

#------------------------------------------------#
#------------------------------------------------#