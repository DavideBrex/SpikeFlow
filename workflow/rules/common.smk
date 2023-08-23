
import os
from snakemake.utils import validate
import pandas as pd

#we start by checking the input files (samples_sheet and config.yaml) to ensure that their format is correct

samples_sheet = pd.read_csv(config["samples_sheet"], dtype={
    "replicate": "Int64",
    "control_replicate": "Int64",
    'spike': "boolean"}, sep=",").set_index(["sample","replicate"], drop=False)

print(samples_sheet)

validate(samples_sheet, schema="../schemas/sampleSheet.schema.yaml")

