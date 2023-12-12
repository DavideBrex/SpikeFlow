import pandas as pd
import os

dict_allsamp = {}
# for each sample, we read the log file and extract the number of called peaks
for f in snakemake.input["logFile"]:
    idName = os.path.basename(f).split(".")[0]
    with open(f, "r") as file:
        calledPeaks = file.read().strip().split(":")[1].strip()
        dict_allsamp[idName] = [calledPeaks]

df_info = pd.DataFrame.from_dict(dict_allsamp, orient="index")
if df_info.empty:
    df_info = pd.DataFrame(columns=["Called Peaks"])
else:
    df_info.columns = ["Called Peaks"]

# Custom headers
headers = [
    "# id: \"Output from Epic2\"",
    "# parent_id: peakSection",
    "# parent_name: \"Peak calling\"",
    "# parent_description: \"Barplots of peak calling tools\"",
    "# section_name: \"Epic2 Peaks\"",
    "# description: \"Peak calling on samples set to broad. Peak caller: Epic2 (https://github.com/biocore-ntnu/epic2) \"",
    "# format: \"tsv\"",
    "# plot_type: \"bargraph\"",
    "# pconfig:",
    "#    id: \"custom_bargraph_w_header\"",
    "#    ylab: \"Number of Peaks\""
]

output_file = snakemake.output["tab"]

# Write headers to the file
with open(output_file, 'w') as file:
    for line in headers:
        file.write(line + '\n')

# Append DataFrame to the file
df_info.to_csv(output_file, sep="\t", index_label="Sample", mode='a')
