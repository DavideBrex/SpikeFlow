import pandas as pd
import os

dict_allsamp = {}
#for each sample, we read the log file and extract the number of reads
for f in snakemake.input["logFile"]:
    idName = os.path.basename(f).split(".")[0]
    with open(f, "r") as file:
        info_sample = file.read().strip().split("\n")
        list_splitBam_values = []
        for i in info_sample:
            list_splitBam_values.append(int(i.split(":")[-1]))

        dict_allsamp[idName] = [
            (list_splitBam_values[1] + list_splitBam_values[2]),
            list_splitBam_values[1],
            list_splitBam_values[2],
            list_splitBam_values[0],
            list_splitBam_values[3],
            list_splitBam_values[4],
            list_splitBam_values[5],
            list_splitBam_values[6],
            ((list_splitBam_values[2] / (list_splitBam_values[1] + list_splitBam_values[2])) * 100),
        ]


df_info = pd.DataFrame.from_dict(dict_allsamp, orient="index")

if df_info.empty:
    df_info = pd.DataFrame(columns=[
        "Tot. mapped Reads",
        "Tot. Sample Reads",
        "Tot. spikeIn Reads",
        "Tot. Common Reads",
        "Tot. unmapped Reads",
        "Tot. QC fail Reads",
        "Tot. secondary Reads",
        "Tot. low mapQ Reads",
        "Percentage spikeIn",
    ])
else:
    df_info.columns = [
        "Tot. mapped Reads",
        "Tot. Sample Reads",
        "Tot. spikeIn Reads",
        "Tot. Common Reads",
        "Tot. unmapped Reads",
        "Tot. QC fail Reads",
        "Tot. secondary Reads",
        "Tot. low mapQ Reads",
        "Percentage spikeIn",
    ]
    #if some columns only contain zeros, we remove them
    df_info = df_info.loc[:, (df_info != 0).any(axis=0)]


df_info.to_csv(snakemake.output["tab"], float_format="%.3f")
df_info.drop(["Tot. mapped Reads", "Percentage spikeIn"], axis=1).to_csv(snakemake.output["tab_multiqc"], float_format="%.2f")