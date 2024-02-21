<img align="left"  width="55%" src="LogoSpikeFlow.png">

<br clear="left"/>

### A snakemake pipeline for the analysis of ChIP-Rx data 

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/DavideBrex/SpikeFlow/workflows/Tests/badge.svg?branch=main)](https://github.com/DavideBrex/SpikeFlow/actions?query=branch%3Amain+workflow%3ATests)
[![DOI](https://zenodo.org/badge/665587139.svg)](https://zenodo.org/doi/10.5281/zenodo.10686979)



If you use this workflow in a paper, don't forget to give credits to the authors. See the [citation](#citation) section.

## About

**SpikeFlow** is a Snakemake-based workflow designed for the analysis of ChIP-seq data with spike-in normalization (i.e. ChIP-Rx). Spike-in controls are used to provide a reference for normalizing sample-to-sample variation. These controls are typically DNA from a different species added to each ChIP and input sample in consistent amounts. This workflow facilitates accurate and reproducible chromatin immunoprecipitation studies by integrating state-of-the-art computational methodologies.
 
**Key Features**:

- **Diverse Normalization Techniques**: This workflow implements three normalization methods,to accommodate varying experimental designs and enhance data comparability.

- **Quality Control**:  Spike-in quality control, to ensure a proper comparison between different experimental conditions

- **Peak Calling**: The workflow incorporates three  algorithms for peak identification, crucial for delineating protein-DNA interaction sites. The user can choose the type of peak to be called: narrow (macs2), broad (epic2), or very-broad (edd). Moreover, the pipeline will merge the called peaks from the replicates and perform peak annotation.

- **BigWig Generation for Visualization**: Normalised BigWig files are generaate for genome-wide visualization, compatible with standard genomic browsers, thus facilitating detailed chromatin feature analyses.

- **Scalability**: Leveraging Snakemake, the workflow ensures an automated, error-minimized pipeline adaptable to both small and large-scale genomic datasets.


**Currently supported genomes**:

- Endogenous: mm9, mm10, hg19, hg38
- Exogenous (spike-in): dm3, dm6, mm10, mm9, hg19, hg38, hs1


## Table of Contents


1. [Installation](#install)
2. [Configuration](#config)
3. [Run the workflow](#run)
4. [Output files](#output)
5. [Troubleshooting](#troubleshooting)
6. [Citation](#citation)


<a name="install"></a>
## Installation 
### Step 1 - Install Snakemake
To run this pipeline, you'll first need to install **Snakemake** (version >= 6.3.0).
Please check the Snakemake documentation on [how to install](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
If you already have *conda* installed, you can simply create a new environment and install Snakemake with:

```
conda create -c bioconda -c conda-forge -n snakemake snakemake
```

For mamba, use the following code:

```
 mamba create -c conda-forge -c bioconda -n snakemake snakemake
```


Once the environment is created, activate it with:
```
conda activate snakemake
```
or
```
mamba activate snakemake
```

### Step 2 - Install Singularity (recommended)

For a fast installation of the workflow, it is recommended to use **Singularity** (compatible with version 3.9.5). This bypasses the need for *Conda* to set up required environments, as these are already present within the container that will be pulled from [dockerhub](https://hub.docker.com/repository/docker/davidebrex/spikeflow/general) with the use of the ```--use-singularity``` flag.

To install singularity check [its website](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).

### Step 3 - Download SpikeFlow

To obtain SpikeFlow, you can:

1.  Download the source code as zip file from the latest [version](https://github.com/DavideBrex/SpikeFlow/releases/latest).

2.  Clone the repository on your local machine. See [here](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) the instructions.

### Step 4 - Test the workflow

Once you obtained the latest version of SpikeFlow, the ```config.yaml```  and the ```samples_sheet.csv```  files are already set to run an installation test. 
You can open them to have an idea about their structure. 
All the files needed for the test are in the ```.test``` folder (on ubuntu, type *ctrl + h* to see hidden files and folders).
To test whether SpikeFlow is working properly, jump directly to the [Run the workflow](#run) section of the documentation.


The usage of this workflow is also described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=DavideBrex%2FSpikeFlow).

<a name="config"></a>
## Configuration

### 1. **Sample Sheet Input Requirements**

Prior to executing the pipeline, you need to prepare a sample sheet containing detailed information about the samples to analyse. You can find an example of this file under ```config/samples_sheet.csv```.
The required format is a comma-separated values (CSV) file, consisting of eight columns and including a header row.
For each sample (row), you need to specify:

| Column Name           | Description                                                                                             	|
|-------------------	|----------------------------------------------------------------------------------------------------------	|
| sample            	| Unique sample name                                                                                       	|
| replicate         	| Integer indicating the number of replicate (if no replicate simply add 1)                                	|
| antibody          	| Antibody used for the experiment (leave empty for Input samples)                                         	|
| control           	| Unique sample name of the control (it has to be specified also in the sample column, but in another row) 	|
| control_replicate 	| Integer indicating the number of replicate for the control sample (if no replicate simply add 1)         	|
| peak_type         	| Can only be equal to: narrow, broad, very-broad. It indicates the type of peak calling to perform        	|
| fastq_1           	| Path to the fastq file of the sample (if paired-end, here goes the forward mate, i.e. R1)                	|
| fastq_2           	| ONLY for paired-end, otherwise leave empty. Path to the fastq file of the reverse mate (i.e. R2)         	|

For the input samples, leave empty the values of the all the columns except for sample, replicate and fastq path(s).

#### Example 1 (single end)

|sample           |replicate|antibody|control        |control_replicate|peak_type|fastq_1                                               |fastq_2|
|-----------------|---------|--------|---------------|-----------------|---------|------------------------------------------------------|-------|
|H3K4me3_untreated|1        |H3K4me3 |Input_untreated|1                |narrow   |fastq/H3K4me3_untreated-1_L1.fastq.gz    |       |
|H3K4me3_untreated|1        |H3K4me3 |Input_untreated|1                |narrow   |fastq/H3K4me3_untreated-1_L2.fastq.gz    |       |
|Input_untreated  |1        |        |               |                 |         |fastq/Input-untreated-1_fastq.gz|       |


> **_⚠️ NOTE:_**  If your sample has **multiple lanes**, you can simple add a new row with the same values in all the columns except for fastq_1 (and fastq_2 if PE). In the table above, H3K4me3_untreated has two lanes

#### Example 2 (paired end)

|sample           |replicate|antibody|control        |control_replicate|peak_type|fastq_1                                               |fastq_2                              |
|-----------------|---------|--------|---------------|-----------------|---------|------------------------------------------------------|-------------------------------------|
|H3K9me2_untreated|1        |H3K9me2 |Input_untreated|1                |very-broad|fastq/H3K9me2_untreated-1_R1.fastq.gz                 |fastq/H3K9me2_untreated-1_R2.fastq.gz|
|H3K9me2_untreated|2        |H3K9me2 |Input_untreated|1                |very-broad|fastq/H3K9me2_untreated-2_R1.fastq.gz                 |fastq/H3K9me2_untreated-2_R2.fastq.gz|
|H3K9me2_EGF      |1        |H3K9me2 |Input_EGF      |1                |very-broad|fastq/H3K9me2_EGF-1_R1.fastq.gz                       |fastq/H3K9me2_EGF-1_R2.fastq.gz      |
|H3K9me2_EGF      |2        |H3K9me2 |Input_EGF      |1                |very-broad|fastq/H3K9me2_EGF-2_R1.fastq.gz                       |fastq/H3K9me2_EGF-2_R2.fastq.gz      |
|Input_untreated  |1        |        |               |                 |         |fastq/Input-untreated-1_R1.fastq.gz                   |fastq/Input-untreated-1_R2.fastq.gz  |
|Input_EGF        |1        |        |               |                 |         |fastq/Input-EGF-1_R1.fastq.gz                         |fastq/Input-EGF-1_R2.fastq.gz        |

> **_⚠️ NOTE:_**  In this case, we have two replicates per condition (untreated and EGF) and the samples are paired-end. However, **mixed situations (some single and some paired-end samples) are also accepted by the pipeline.**

### 2. **Config file**

The last step before running the workflow is to adjust the parameters in the config file (```config/config.yaml```). The file is written in YAML (Yet Another Markup Language), which is a human-readable data serialization format. It contains key-value pairs that can be nested to multiple leves.

#### *Reference and exogenous (spike-in) genomes*

To execute the pipeline, it's essential to specify both *endogenous* and *exogenous* species; for example, use Drosophila (dm16) as the exogenous and Human (hg38) as the endogenous species.
If bowtie v1 genome indexes are already available, you should input their paths (ending with the index files prefix) in the 'resources' section of the pipeline configuration. This setup ensures proper alignment and processing of your genomic data. **PLEASE NOTE**, the index must be created with bowtie  v1.3.0.


```yaml
resources:
    ref:
        index: /path/to/hg38.bowtie.index/indexFilesPrefix
        #blacklist regions 
        blacklist: ".test/data/hg38-blacklist.v2.bed"

    ref_spike:
        index: /path/to/dm16.bowtie.index/indexFilesPrefix
```

If you don't have the bowtie genome indexes readily available, the pipeline can generate them for you. To facilitate this, you'll need to specify the ucsc genome name (e.g. hg38, mm10, etc) for both the reference and spike-in species. You can find the the genome assembly on the [UCSC Genome Browser](https://genome-euro.ucsc.edu/cgi-bin/hgGateway).

```yaml
    ref:
        assembly: hg38
    ref_spike:
        spike_assembly: dm6
```

> **_⚠️ NOTE:_**  For the endogenous genome,  it's important to also include the path to blacklisted regions.  These regions, often associated with sequencing artifacts or other anomalies, can be downloaded from the Boyle Lab's Blacklist repository on GitHub. You can access these blacklisted region files [here](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)


#### *Normalization*

In this field you can choose the type of normalization to perform on the samples. The available options are:

- **RAW**: This is a RPM normalization, i.e. it normalizes the read counts to the total number of reads in a sample, measured per million reads. This method is straightforward but does not account for spike-in. 

- **Orlando**: Spike-in normalization as described in [Orlando et al 2014](https://pubmed.ncbi.nlm.nih.gov/25437568/). Also reffered as Reference-adjusted Reads Per Million (RRPM). It does not incorporate input data in the normalization process.

- **RX-Input** (default): RX-Input is a modified version of the Orlando normalization that accounts for the total number of reads mapped to the spike-in in both the ChIP and input samples. This approach allows for more accurate normalization by accounting for variations in both immunoprecipitation efficiency and background noise (as represented by the input). See [Fursova et al 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6561741/#bib42) for further details.


```yaml
normalization_type: "RX-Input"
```


#### *Required options*

When configuring your pipeline based on the chosen reference/endogenous genome (like mm10 or hg38), two essential options need to be set:

- **effective genome length**: This is required by deeptools to generate the bigWig files. The value of this parameter is used by the program to adjust the mappable portion of the genome, ensuring that the read densities represented in the BigWig files accurately reflect the underlying biological reality. You can find the possible values for this parameter in the deeptools [documentation](https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html).

- **chrom sizes**: To achieve accurate peak calling, it's important to use the correct chromosome sizes file. The supported genomes' chromosome sizes are available under ```resources/chrom_size```.  **Make sure to select the file that corresponds to your chosen genome**.


#### *Other (optional) parameters*

-  To direct Snakemake to save all outputs in a specific directory, add the desired path in the config file: ```output_path: "path/to/directory"```.

- **Broad Peak Calling:** For samples requiring broad peak calling, adjust the effective genome fraction as per the guidelines on this [page]( https://github.com/biocore-ntnu/epic2/blob/master/epic2/effective_sizes/hg38_50.txt). The *'effective genome size'* mentioned on the GitHub page depends on the read length of your samples.

- **Very Broad Peak Callling:** If you have samples that will undergo very-broad peak calling, please check the log files produced by EDD. This because the tool might fail if it can not accurately estimate the parameters for the peak calling. In this case, you can tweak the parameters in the EDD config file, which is in the config directory (```config/edd_parameters.conf```). For more information about EDD parameters tuning see the [documentation](https://github.com/CollasLab/edd).

- **Trimming Option:** Trimming can be skipped by setting the respective flag to false.

- **P-Value Adjustment for Peak Calling:** Modify the p-values for peak calling in the config file. This applies to different peak calling methods: narrow (macs2), broad (epic2), or very-broad (edd).

- **Peak Merging from Replicates:** While merging peaks from replicates, the size can be adjusted as described [here](https://github.com/rhysnewell/ChIP-R))

- **Peak Annotation Threshold:** The default setting annotates a peak within ±2500 bp around the promoter region.

<a name="run"></a>
## Run the workflow

To execute the pipeline, make sure to be in the main  **Snakemake working directory**, which includes subfolders like 'workflow', 'resources', and 'config'. 

The workflow can be operated in two ways: using Conda alone, or a combination of Conda and Singularity (**recommended**). 
After obtaining a copy of the workflow on your machine, you can verify its proper functioning by executing one of the two commands below. 
The ```config``` and ```sample_sheet``` files come pre-configured for a test run.

#### Conda and Singularity (recommended)

```bash
snakemake --cores 10 --use-conda --use-singularity
```

First, the singularity container will be pulled from DockerHub and then the workflow will be executed. To install sigularity, see the [installation](#install) section.

#### Conda only

```bash
snakemake --cores 10 --use-conda 
```
This will install all the required conda envs (it might take a while, just for the first execution).

#### Snakemake flags

- ```--cores```: adjust this number (here set to 10) based on your machine configuration
- ```-n```: add this flag  to the command line for a "dry run," which allows Snakemake to display the rules that it would execute, without actually running them. 
- ```--singularity-args "-B /shares,/home -e"```: only with ```--use-singularity``` and if you are running SpikeFlow from a server that mounts /shares and /home disks (where you have your working dir and files).

To execute the pipeline on a HPC cluster, please follow [these guidelines](https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#cluster-execution).

If you are using **Snakemake version $\ge$ 8**, the command line arguments have [different names](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#containerization-of-conda-based-workflows). In this case, run the workflow with:

```bash
snakemake --cores --software-deployment-method conda apptainer
# or the shorthand version
snakemake --cores 10 --sdm conda apptainer
```
<a name="output"></a>
## Output files


All the outputs of the workflow are stored in the ```results``` folder. Additionally, in case of any errors during the workflow execution, the log files are stored within the ```results/logs``` directory.


The main outputs of the workflow are:

- **MultiQC Report** 

    Consolidates various quality control (QC) metrics:
    - Basic QC with FastQC: Evaluates basic quality metrics of the sequencing data.
    - Phantom Peak Qual Tools: Provides NSC and RSC values, indicating the quality and reproducibility of ChIP samples. NSC measures signal-to-noise ratio, while RSC assesses enrichment strength.
    - Fingerprint Plots: Visual representation of the sample quality, showing how reads are distributed across the genome.
    - Spike-in QC: Reveals the number of reads aligned to both the endogenous and exogenous genomes, crucial for evaluating spike-in normalization effectiveness.
    - Peak Calling Data: Displays the number of peaks called per sample for each method (MACS2, EPIC2, EDD).

     ```results/QC/multiqc/multiqc_report.html```

- **Normalized BigWig Files**: 
    - Essential for visualizing read distribution and creating detailed heatmaps.

    ```results/bigWigs/```

-  **Peak Files and Annotation**:
    - Provides called peaks for each peak calling method, along with merged peaks for replicate samples.
    - Peak annotation using ChIPseeker, resulting in two files for promoter and distal peaks for each sample

    ```/results/peakCalling/```


<a name="troubleshooting"></a>
## Troubleshooting

1. When you run SpikeFlow with Singularity ( ```--use-singularity```), you might get an error if you set the  ```-n``` flag. This happens ONLY at the first execution of the workflow. Remove the flag and it should work.

2. The ``` --use-singularity```  option temporary requires about 7 GB of disk space to build the image from Docker Hub. If your ```/tmp``` directory is full, you'll encounter a ```No space left on device``` error. To avoid this, change the Singularity temp directory to a different disk by setting the ```SINGULARITY_TMPDIR``` environment variable. More details are available in the [Singularity guide on temporary folders](https://docs.sylabs.io/guides/latest/user-guide/build_env.html#temporary-folders).

<a name="citation"></a>
## Citation

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository:
https://github.com/DavideBrex/SpikeFlow

**Author**:
- Davide Bressan ([@DavideBrex](https://twitter.com/BrexDavide))

