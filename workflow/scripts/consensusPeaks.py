import subprocess
import os
from collections import defaultdict
import shutil
import sys

sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout

# Assuming snakemake.input is a list of peak files and snakemake.output is the output file path
peak_files = snakemake.input
output_path = snakemake.output[0]
antibody = snakemake.params.antibody
sampleNamesToUse = snakemake.params.sampleNamesToUse

# Minimum number of replicates a peak must be present in to be included
min_num_reps = snakemake.params.min_num_reps

#check if bedtools is installed
assert shutil.which('bedtools') is not None, "ERROR! bedtools is not installed, can't create consensus peaks. Please check conda env."

#we want to merge only the peaks that are from the same antibody
peak_files_subset= [p for p in peak_files if os.path.basename(p).rsplit('_',1)[0] in sampleNamesToUse]
print("Merging peaks from the following files:{}".format(peak_files_subset))

# Create a temporary file for merged peaks
temp_merged_peaks_path = output_path + ".tmp_merged_peaks.bed"
with open(temp_merged_peaks_path, 'w') as temp_file:
    cmd = f"cat {' '.join(peak_files_subset)} | sort -k1,1 -k2,2n | bedtools merge -i stdin -d 150 -c 4,5,6,7,8,9,10 -o collapse,mean,collapse,mean,collapse,collapse,collapse"
    subprocess.run(cmd, shell=True, stdout=temp_file)

# Process the merged peaks
with open(temp_merged_peaks_path) as mergedPeaks, open(output_path, 'w') as outfile:
    for line in mergedPeaks:
        cols = line.strip().split('\t')
        # The peak names are in the 4th column (0-based indexing), split by commas
        peak_names = cols[3].split(',')
        # Using a set to count unique replicates, assuming replicate info is encoded in the peak names
        unique_replicates = set(peak_name.rsplit('_', 2)[0] for peak_name in peak_names)  
        # Filter based on min_num_reps
        if len(unique_replicates) >= min_num_reps:
            outfile.write(line)

# Optionally, remove the temporary merged peaks file
os.remove(temp_merged_peaks_path)

print("Consensus peaks generation completed.")