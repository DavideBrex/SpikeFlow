import pandas as pd
import os

norm_factor_files = snakemake.input["normFactorsFiles"]
normType = snakemake.params['normType']


dict_allsamp = {}
#for each sample, we read the log file and extract the number of reads
for f in snakemake.input["logFile"]:
    idName = os.path.basename(f).split(".")[0]
    with open(f, "r") as file:
        info_sample = file.read().strip().split("\n")

        #let's also get the norm factor
        with open(f.replace('.removeSpikeDups','.normFactor'), "r") as file:
            norm_factor = file.read().strip()
            norm_factor = round(float(norm_factor.split(":")[-1]),2)

        list_splitBam_values = []
        for i in info_sample:
            list_splitBam_values.append(int(i.split(":")[-1]))

        dict_allsamp[idName] = [
            (list_splitBam_values[0] + list_splitBam_values[1]),
            list_splitBam_values[0],
            list_splitBam_values[1],
            list_splitBam_values[2],
            ((list_splitBam_values[1] / (list_splitBam_values[0] + list_splitBam_values[1])) * 100),
            norm_factor,
        ]


df_info = pd.DataFrame.from_dict(dict_allsamp, orient="index")

if df_info.empty:
    df_info = pd.DataFrame(columns=[
        "Tot. mapped Reads",
        "Tot. Sample Reads",
        "Tot. spikeIn Reads",
        "Tot. low mapQ Reads",
        "Percentage spikeIn",
        "Normalization factor",
    ])
else:
    df_info.columns = [
        "Tot. mapped Reads",
        "Tot. Sample Reads",
        "Tot. spikeIn Reads",
        "Tot. low mapQ Reads",
        "Percentage spikeIn",
        "Normalization factor",
    ]
    #if some columns only contain zeros, we remove them
    df_info = df_info.loc[:, (df_info != 0).any(axis=0)]

# Custom headers
headers = [
    "# id: \"Mapped_reads_tab\"",
    "# parent_description: \"Table of mapped reads\"",
    "# section_name: \"Reads Table\"",
    "# description: \"Per sample reads and normalization factor calculated with {} normalisation.\"".format(normType),
    "# format: \"tsv\"",
    "# plot_type: \"table\"",
]

output_file_mqc = snakemake.output["tab_multiqc"]

# Write headers to the file
with open(output_file_mqc, 'w') as file:
    for line in headers:
        file.write(line + '\n')


df_info.to_csv(snakemake.output["tab"], float_format="%.3f")
df_info.to_csv(output_file_mqc,  sep="\t", index_label="Sample", mode='a')

#.drop(["Tot. mapped Reads", "Percentage spikeIn"], axis=1)