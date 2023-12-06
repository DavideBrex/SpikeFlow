import pandas as pd
import os

dict_allsamp = {}
#for each sample, we read the log file and extract the number of reads
for f in snakemake.input["logFile"]:
    idName = os.path.basename(f).split(".")[0]
    with open(f, "r") as file:
        info_sample = file.read().strip().split("\n")
        common_nr = int(info_sample[0].split(":")[-1])
        sample_nr = int(info_sample[1].split(":")[-1])
        spike_nr = int(info_sample[2].split(":")[-1])
        dict_allsamp[idName] = [
            (sample_nr + spike_nr),
            sample_nr,
            spike_nr,
            common_nr,
            ((spike_nr / (sample_nr + spike_nr)) * 100),
        ]

df_info = pd.DataFrame.from_dict(dict_allsamp, orient="index")
df_info.columns = [
    "Total mapped Reads",
    "Tot. Sample Reads",
    "Tot. spikeIn Reads",
    "Common Reads",
    "Percentage spikeIn",
]
df_info.to_csv(snakemake.output["tab"], float_format="%.3f")
df_info.iloc[:, 1:4].to_csv(snakemake.output["tab_multiqc"], float_format="%.2f")