## Configuration

### 1. **Sample Sheet Input Requirements**

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

### *Reference and exogenous (spike-in) genomes*

To execute the pipeline, it's essential to specify both *endogenous* and *exogenous* species; for example, use Drosophila (dm16) as the exogenous and Human (hg38) as the endogenous species.
Regarding alignment, you have the option to select between Bowtie and Chromap for the alignment process, with Bowtie set as the default. If genome indexes are already available, you should input their paths in the 'resources' section of the pipeline configuration. This setup ensures proper alignment and processing of your genomic data.

```yaml
resources:
    ref:
        index: /path/to/hg38.bowtie.index
    ref_spike:
        index: /path/to/dm16.bowtie.index
```

If you don't have the genome indexes readily available, the pipeline can generate them for you. To facilitate this, you'll need to specify the parameters shown below for both the reference and spike-in species. These parameters include the species name, the Ensembl release, and the genome build version. You can find all of this information on the Ensembl [website](https://www.ensembl.org/index.html).

```yaml
    # Ensembl species name
    species: homo_sapiens
    # Ensembl release
    release: 101
    # Genome build
    build: GRCh38 
```

> **_⚠️ NOTE:_**  For the endogenous genome,  it's important to also include the path to blacklisted regions.  These regions, often associated with sequencing artifacts or other anomalies, can be downloaded from the Boyle Lab's Blacklist repository on GitHub. You can access these blacklisted region files [here](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)


### *Normalization*

In this field you can choose the type of normalization to perform on the samples. The available options are:

- **RAW**: This is a RPM normalization, i.e. it normalizes the read counts to the total number of reads in a sample, measured per million reads. This method is straightforward but does not account for spike-in. 

- **Orlando**: Standard Spike-in normalization as described in [Orlando et al 2014](https://pubmed.ncbi.nlm.nih.gov/25437568/). It does not incorporate input data in the normalization process.

- **RX-Input** (default): RX-Input is a modified version of the Orlando normalization that accounts for the total number of reads mapped to the spike-in in both the ChIP and input samples. This approach allows for more accurate normalization by accounting for variations in both immunoprecipitation efficiency and background noise (as represented by the input).


```yaml
normalization_type: "RX-Input"
```


### *Required options*

When configuring your pipeline based on the chosen reference/endogenous genome (like mm10 or hg38), two essential options need to be set:

- **effective genome length**: This is required by deeptools to generate the bigWig files. The value of this parameter is used by the program to adjust the mappable portion of the genome, ensuring that the read densities represented in the BigWig files accurately reflect the underlying biological reality. You can find the possible values for this parameter in the deeptools [documentation](https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html).

- **chrom sizes**: To achieve accurate peak calling, it's important to use the correct chromosome sizes file. The supported genomes' chromosome sizes are available under ```resources/chrom_size```.  **Make sure to select the file that corresponds to your chosen genome**.


### *Other (optional) parameters*

-  To direct Snakemake to save all outputs in a specific directory, add the desired path in the config file: ```output_path: "path/to/directory"```.

- **Broad Peak Calling:** For samples requiring broad peak calling, adjust the effective genome fraction as per the guidelines on this [page]( https://github.com/biocore-ntnu/epic2/blob/master/epic2/effective_sizes/hg38_50.txt). The *'effective genome size'* mentioned on the GitHub page depends on the read length of your samples.

- **Trimming Option:** Trimming can be skipped by setting the respective flag to false.

- **P-Value Adjustment for Peak Calling:** Modify the p-values for peak calling in the config file. This applies to different peak calling methods: narrow (macs2), broad (epic2), or very-broad (edd).

- **Peak Merging from Replicates:** While merging peaks from replicates, the size can be adjusted as described [here](https://github.com/rhysnewell/ChIP-R))

- **Peak Annotation Threshold:** The default setting annotates a peak within ±2500 bp around the promoter region.