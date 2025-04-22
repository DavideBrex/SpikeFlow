import pandas as pd
import os
import sys

sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout


dict_allsamp = {}
# for each sample, we read the log file and extract the annotations
for f in snakemake.input:
    #check if the peakCallingNorm is in the path of file, then we attach "PeakNorm" to the idName
    if "peakCallingNorm" in f:
        idName = os.path.basename(f).split("_annotInfo.")[0] + "_PeakNorm"
    else:
        idName = os.path.basename(f).split("_annotInfo.")[0]
    with open(f, "r") as file:
        headerRow = file.readline()
        annotRow = file.readline()
        dict_allsamp[idName] = annotRow.split()

df_info = pd.DataFrame.from_dict(dict_allsamp, orient="index")

if df_info.empty:
    df_info = pd.DataFrame(columns=[""])
else:
    df_info.columns = headerRow.split()

# save DataFrame to the file
df_info.to_csv(snakemake.output["tab"], sep="\t", index_label="Sample")