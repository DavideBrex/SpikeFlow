
import os
from snakemake.utils import validate
import pandas as pd


#we start by checking the input files (samples_sheet and config.yaml) to ensure that their format is correct

samples = pd.read_csv(config["samples_sheet"], sep=",", dtype = str).set_index(["sample","replicate"], drop=False)

validate(samples, schema="../schemas/sampleSheet.schema.yaml")