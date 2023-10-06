def normalization_factor(wildcards):
    """
    Read the log message from the cleaning of bam files and compute normalization factor
    By using the log file we avoid to read the bam file with pysam just to get # of aligned reads
    """

    with open(f"results/logs/spike/{wildcards.id}.removeSpikeDups", 'r') as file:
        info_sample=file.read().strip().split("\n")

        #we need the information also from the input (if it is not the sample an input itself)
        inputSamp=sample_to_input[wildcards.id]
        if  pd.isna(inputSamp): 
            Nsample=int(info_sample[1].split(":")[-1]) #number of aligned reads in sample
            alpha=1/(Nsample)*1000000 #normalization factsor
        else:
            #open input log file
            with open(f"results/logs/spike/{inputSamp}.removeSpikeDups", 'r') as file:
                info_input=file.read().strip().split("\n")

                gamma=int(info_input[2].split(":")[-1])/int(info_input[1].split(":")[-1]) #ratio spike/samples in input 
                Nspike=int(info_sample[2].split(":")[-1]) #number of spike reads in sample
                alpha=gamma/Nspike*1000000 #normalization factor
    


    if(is_single_end(wildcards.id)):
        return "--scaleFactor "+str(round(alpha, 4))+" --extendReads "+str(config["params"]["deeptools"]["read_extension"])
    else:
        return "--scaleFactor "+str(round(alpha, 4))+" --extendReads"
    #TO DO: add log file with the norm factors stored



rule bam2bigwig_general:
    input: 
        bam="results/bam/{id}.clean.bam",
        logFile="results/logs/spike/{id}.removeSpikeDups",
        logFileInput=lambda wildcards: "results/logs/spike/{}.removeSpikeDups".format(sample_to_input[wildcards.id])
        if not pd.isna(sample_to_input[wildcards.id])
        else "results/logs/spike/{id}.removeSpikeDups",
        blacklist=config["resources"]["ref"]["blacklist"],
    output:
        out='results/bigWigs/{id}.bw',
    params: 
        extra=lambda wildcards: normalization_factor(wildcards),
        effective_genome_size=config["params"]["deeptools"]["effective_genome_length"],
    log: 'results/logs/bam2bigwig/{id}.log'
    threads: 5
    wrapper: "v2.6.0/bio/deeptools/bamcoverage"


#ruleorder: clean_spike > bam2bigwig_general
