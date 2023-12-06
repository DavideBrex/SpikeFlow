# Snakemake workflow: `ChIP-Rx-snakemake`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/DavideBrex/ChIP-Rx-snakemake/workflows/Tests/badge.svg?branch=main)](https://github.com/DavideBrex/ChIP-Rx-snakemake/actions?query=branch%3Amain+workflow%3ATests)

## **Under development**

A Snakemake workflow for the anlysis of ChIP-Rx data, i.e ChIP-Seq with reference exogenous genome spike-in normalization

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository. 

## About

**ChIP-Rx-snakemake** is a Snakemake-based workflow designed for the analysis of ChIP-seq data with spike-in normalization (i.e. ChIP-Rx). Spike-in controls are used to provide a reference for normalizing sample-to-sample variation. These controls are typically DNA from a different species added to each ChIP and input sample in consistent amounts. This workflow facilitates accurate and reproducible chromatin immunoprecipitation studies by integrating state-of-the-art computational methodologies.

Key Features:
- **Diverse Normalization Techniques**: This workflow implements three normalization methods,to accommodate varying experimental designs and enhance data comparability.

- **Quality Control**:  Spike-in quality control, to ensure a proper comparison between different experimental conditions

- **Peak Calling**: The workflow incorporates three  algorithms for peak identification, crucial for delineating protein-DNA interaction sites. The user can choose the type of peak to be called: narrow (macs2), broad(epic2), or very-broad(edd). Moreover, the pipeline will merge the called peaks from the replicates and perform peak annotation.

- **BigWig Generation for Visualization**: Normalised BigWig files are generaate for genome-wide visualization, compatible with standard genomic browsers, thus facilitating detailed chromatin feature analyses.

- **Scalability**: Leveraging Snakemake, the workflow ensures an automated, error-minimized pipeline adaptable to both small and large-scale genomic datasets.


## Usage

### Set up Snakemake
To run this pipeline, you'll first need to install **Snakemake** (version >= 6.3.0).
If you already have conda installed, you can simply create a new environment:

```
conda create -c bioconda -c conda-forge -n snakemake snakemake
```

In case you do not have conda, we reccomand to install Mamba (a drop-in replacement for conda). See [here](https://github.com/conda-forge/miniforge#mambaforge) the installation steps.

Once Mamba is installed, run


```
 mamba create -c conda-forge -c bioconda --name snakemake
```

Once the environment is created, activate it with:
```
conda activate snakemake
```
or
```
mamba activate snakemake
```

### Download the workflow

To obtain the Snakemake workflow, you can:
- Clone it on your local machine:
    1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
    2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, in the folder where you want to perform the data analysis.

- Download the source code as zip file from this page (code button)


The usage of this workflow is also described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=DavideBrex%2FChIP-Rx-snakemake).


### Configuration

#### Sample Sheet Input Requirements

Prior to executing the pipeline, you need to prepare a sample sheet containing detailed information about the samples to analyse. You can find an example of this file under ```config/samples_sheet.csv```.
The required format is a comma-separated values (CSV) file, consisting of eight columns and including a header row.
For each sample (row), you need to specify:
1. **sample**: Unique sample name 
2. **replicate**: Integer indicating the number of replicate (if no replicate simply add 1)
3. **antibody**: antibody used for the experiment
4. **control**: Unique sample name of the control (it has to be specified also in the sample column, but in another row)
5. **control_replicate**: Integer indicating the number of replicate for the control sample
6. **peak_type**: can only be equal to: narrow, broad, very-broad. It indicates the type of peak calling to perform
7. **fastq_1**: path to the fastq file of the sample (if paired-end, here goes the forward mate, i.e. R1)
8. **fastq_2**: ONLY for paired-end, otherwise leave empty. Path to the fastq file of the reverse mate (i.e. R2)

For the input samples, leave empty the values of the all the columns except for sample, replicate and fastq path(s).

##### Example 1 (single end)

|sample           |replicate|antibody|control        |control_replicate|peak_type|fastq_1                                               |fastq_2|
|-----------------|---------|--------|---------------|-----------------|---------|------------------------------------------------------|-------|
|H3K4me3_untreated|1        |H3K4me3 |Input_untreated|1                |narrow   |fastq/H3K4me3_untreated-1_L1.fastq.gz    |       |
|H3K4me3_untreated|1        |H3K4me3 |Input_untreated|1                |narrow   |fastq/H3K4me3_untreated-1__L2.fastq.gz    |       |
|Input_untreated  |1        |        |               |                 |         |fastq/Input-untreated-1_fastq.gz|       |


> **_⚠️ NOTE:_**  If your sample has multiple lanes, you can simple add a new row with the same values in all the columns except for fastq_1 (and fastq_2 if PE). In the table above, H3K4me3_untreated has two lanes

##### Example 2 (paired end)

|sample           |replicate|antibody|control        |control_replicate|peak_type|fastq_1                                               |fastq_2                              |
|-----------------|---------|--------|---------------|-----------------|---------|------------------------------------------------------|-------------------------------------|
|H3K9me2_untreated|1        |H3K9me2 |Input_untreated|1                |very-broad|fastq/H3K9me2_untreated-1_R1.fastq.gz                 |fastq/H3K9me2_untreated-1_R2.fastq.gz|
|H3K9me2_untreated|2        |H3K9me2 |Input_untreated|1                |very-broad|fastq/H3K9me2_untreated-2_R1.fastq.gz                 |fastq/H3K9me2_untreated-2_R2.fastq.gz|
|H3K9me2_EGF      |1        |H3K9me2 |Input_EGF      |1                |very-broad|fastq/H3K9me2_EGF-1_R1.fastq.gz                       |fastq/H3K9me2_EGF-1_R2.fastq.gz      |
|H3K9me2_EGF      |2        |H3K9me2 |Input_EGF      |1                |very-broad|fastq/H3K9me2_EGF-2_R1.fastq.gz                       |fastq/H3K9me2_EGF-2_R2.fastq.gz      |
|Input_untreated  |1        |        |               |                 |         |fastq/Input-untreated-1_R1.fastq.gz                   |fastq/Input-untreated-1_R2.fastq.gz  |
|Input_EGF        |1        |        |               |                 |         |fastq/Input-EGF-1_R1.fastq.gz                         |fastq/Input-EGF-1_R2.fastq.gz        |

> **_⚠️ NOTE:_**  In this case, we have two replicates per condition (untreated and EGF) and the samples are paired-end. However, mixed situations (some single and some paired-end samples) are also accepted by the pipeline.

#### Config file

The last step before running the analysis is to adjust of the config file (```config/config.yaml```), which stores the parameters used throught the workflow. The file is written in YAML (Yet Another Markup Language), which is a human-readable data serialization format. It contains key-value pairs that can be nested to multiple leves.

##### Reference and exogenous genomes
To run the pipline you need to set both endogenous and exogenous species (e.g. exogenous is drosophila and endogenous human).

For the alignment, you can set whether to use bowtie or chromap (default is bowtie).

In case you already have the genome indexes, you should add their paths under the field resources

```yaml
resources:
    ref:
        index: /path/to/hg38.bowtie.index
    ref_spike:
        index: /path/to/dm16.bowtie.index
```

In case you do not have the indexes at hand, do not worry, the pipeline will generate them for you.
To do so, you need to set the species, the ensembl release and build parameteres both for reference and spike-in. You can find these info on the ensembl [website](https://www.ensembl.org/index.html).

```yaml
    # Ensembl species name
    species: homo_sapiens
    # Ensembl release
    release: 101
    # Genome build
    build: GRCh38 
```

> **_⚠️ NOTE:_**  For the endogenous genome, please also add the path to blacklisted regions, that you can easily download from [here](endogehttps://github.com/Boyle-Lab/Blacklist/tree/master/listsnous)


##### Normalization type
In this field you can choose the type of normalization to perform on the samples. The available options are:
- **RAW**: This is a RPM normalization, i.e. it normalizes the read counts to the total number of reads in a sample, measured per million reads. This method is straightforward but does not account for spike-in controls. 
- **Orlando**: Standard Spike-in normalization as described in Orlando et al 2014. It does not incorporate input data in the normalization process, 
- **RX-Input**: RX-Input is a modified version of the Orlando normalization that accounts for the total number of reads mapped to the spike-in in both the ChIP and input samples. This approach allows for more accurate normalization by accounting for variations in both immunoprecipitation efficiency and background noise (as represented by the input).


```yaml
normalization_type: "RX-Input"
```


##### Required options

Based on your reference genome (e.g. mm10, hg38),there are twp required options to set:

- **effective genome length**: This is required by deeptools to generate the bigWig files. The value of this parameter is used by the program to adjust the mappable portion of the genome, ensuring that the read densities represented in the BigWig files accurately reflect the underlying biological reality. See [here](https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html) the possibile values. 

- **chrom sizes**: for an accurate peak calling, please change the chromosome sizes file accordingly. The supported genomes chrom sizes are avaialbe under ```resources/chrom_size```


##### Other (optional) parameters

- In case broad peak calling is set for some samples, please change the effective genome fraction according to this [page]( https://github.com/biocore-ntnu/epic2/blob/master/epic2/effective_sizes/hg38_50.txt). In the github page, it is indicated as effective genome size, and you need to know the read length  of your samples in order to set it

- It is possible to do not perform trimmign by setting the flag to false

- The p-values for the peak calling can be changed in the config file (narrow: macs2, broad:epic2, or very-broad:edd).

- While merging the peaks from the replicates, size from Chip-R can be adjusted (see [here](https://github.com/rhysnewell/ChIP-R))

- The threshold to annotate a peak is set to ±2500 bp around the promoter


### Run the workflow

To run the pipline type on the terminal (from within the folder with all the workflow subfolder):

```
snakemake --cores all --use-conda 
```

For running the workflow while using a combination of conda and singularity for software deployment, run Snakemake with

```
snakemake --cores all --use-conda --use-singularity
```

### Output files

to write

## Authors

- Davide Bressan (@DavideBrex)

