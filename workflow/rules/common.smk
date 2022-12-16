import pandas as pd
import os
import time
from snakemake.utils import validate


# Pipeline setup --------------------------------------
version = open(os.path.join(workflow.basedir, "..", "VERSION"), "r").read()
pipe_log = os.path.join(os.getcwd(), "PIPELINE_STATUS")


# Validating config ----------------------------------
# validate(config, schema="../schema/config.schema.yaml")


# Loading and validating samples ---------------------
sample_path = config["samples"]
samples = pd.read_csv(sample_path, index_col="sample", sep="\t", engine="python")
validate(samples, schema="../schema/samples.schema.yaml")
samples.index = samples.index.astype("str", copy=False)

# Loading and validating reference genomes -----------
# genomes_path = config['genomes']
# genomes = pd.read_csv(genomes_path, index_col="organism", sep="\t", engine="python")
# validate(genomes, schema="../schema/genomes.schema.yaml")
# genomes.index = genomes.index.astype("str", copy=False)

# General puprose functions --------------------------
def get_local_time():
    return time.asctime(time.localtime(time.time()))


# Input functions ------------------------------------
def get_fastq(wildcards, read_pair="fq1"):
    return samples.loc[(wildcards.sample), [read_pair]].dropna()[0]

