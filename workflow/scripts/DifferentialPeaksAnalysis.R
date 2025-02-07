# Load the DESeq2 package
log <- file(snakemake@log[[1]], open= "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages(require(DESeq2))
suppressPackageStartupMessages(require(tidyverse))


# Read the params
antibody_contrastToApply <- gsub("_diffPeaks.tsv", "", basename(snakemake@output[['diffTab']]))
antibody <-sub("_.*", "", antibody_contrastToApply)
contrastToApply <-  sub("^[^_]*_", "", antibody_contrastToApply)
cat("Antibody: ")
print(antibody)
cat("Contrast: ")
print(contrastToApply)
padjCutoff <- snakemake@params[["padjCutoff"]]
log2FCcutoff <- snakemake@params[["log2FCcutoff"]]
outdir <- snakemake@params[["outdir"]]
normMethod <- snakemake@params[['normMethod']]

normFactorsFiles <- snakemake@input[['logFile']]

#we check whether we are dealing with peaks obtained from raw or normalised peak calling
if (tail(strsplit(outdir, '/')[[1]], n=2)[1] == "NormalisedPeaks"){
  peak_type="Norm"
} else{
  peak_type="Raw"
}
  
# Read the file with raw counts
cat("Reading raw counts\n")
countsTab <- read.table(snakemake@input[['rawReadsOnPeaks']], check.names = F, header = TRUE, sep = "\t")
#remove the chr, star, end columns
countsTab <- countsTab[,-c(2,3,4)]
#set rownames to region
countsTab <- countsTab %>% column_to_rownames("region")

colNames <- colnames(countsTab)

#if the number of samples is 0 or 1 we stop..no sense to perform the analysis
if (length(colNames) == 0){
  write.table(data.frame(), file = snakemake@output[['diffTab']], sep = "\t", quote = F, row.names = F, col.names = T)
  cat("No counts found")
  quit(save =  "no",  status=0)
} else if (length(colNames) == 1){
  write.table(data.frame(), file = snakemake@output[['diffTab']], sep = "\t", quote = F, row.names = F, col.names = T)
  cat("Only one sample found, skipping differential analysis")
  quit(save =  "no",  status=0)
}
#if the count table is empty, we stop
if (nrow(countsTab) == 0){
  write.table(data.frame(), file = snakemake@output[['diffTab']], sep = "\t", quote = F, row.names = F, col.names = T)
  cat("Count table is empty, skipping differential analysis")
  quit(save =  "no",  status=0)
}

# Extract the parameters
leftContrast <- strsplit(contrastToApply, "_vs_")[[1]][1]
rightContrast <- strsplit(contrastToApply, "_vs_")[[1]][2]


# Extract the group info
#remove rep
colNames <- gsub("-rep\\d+$", "", colNames)
#extract group
group <- gsub(".*_", "", colNames)

# now read normFactors from spike-in
cat("Reading normFactors\n")


# Read the first line from each file
normFactorsList <- sapply(normFactorsFiles, function(path) {
  # Open the file and read the first line
  line <- readLines(path, n = 1)
  # Extract the normFactor
  normFactor <- as.numeric(trimws(strsplit(line, ":")[[1]][2]))
  normFactor
})

# Print the first lines
names(normFactorsList) <- gsub(".normFactor", "", basename(normFactorsFiles))
print(normFactorsList)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#DESeq2
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
cat("Running DESeq2\n")
# Create the colData with the group information
colData <- data.frame(sample = colnames(countsTab), condition = group)
rownames(colData) <- colData$sample

stopifnot(all(rownames(colData) %in% colnames(countsTab)))

print(colData)
#in case there is just one sample per group we cannot run deseq2, we simply calculate the logFC

if (sum(leftContrast == colData$condition) == 1 && sum(rightContrast == colData$condition) == 1) {
  cat("Only one sample per group, skipping DESeq2 and only calculate log2FC!\n")
  dds <- DESeqDataSetFromMatrix(countData = countsTab, colData = colData, design = ~ 1)
  stopifnot(colnames(assay(dds)) == colnames(normFactorsList))
  #we invert the normFactorsList, since this is the similar to adding a row to the countTab with the 
  #number of spike-ins reads and then run:  dds <- estimateSizeFactors(dds, controlGenes=1)
  #see here: https://support.bioconductor.org/p/9147716/
  sizeFactors(dds) <- 1/normFactorsList

  # Calculate vst and logFC
  tableVST <- counts(dds,normalized = T) %>% as.data.frame() 
  sample1 <- colData[colData$condition == leftContrast,]$sample
  sample2 <- colData[colData$condition == rightContrast,]$sample
  tableVST$log2FC <- log2(tableVST[,sample1]+0.01) - log2(tableVST[,sample2]+0.01)
  cat("Saving results\n")
  tableVST <- tableVST %>% rownames_to_column("region")
  write.table(tableVST[,c('region',sample1, sample2, "log2FC")], file = snakemake@output[['diffTab']], sep = "\t", quote = F, row.names = F, col.names = T)

  #------------------------------------------------------------------------------#
  #------------------------------------------------------------------------------#
  #scatter plot of logFC
  #------------------------------------------------------------------------------#
  #------------------------------------------------------------------------------#
  minLim <- log2(min(c(tableVST[,sample1], tableVST[,sample2]))+0.01)
  maxLim <- log2(max(c(tableVST[,sample1], tableVST[,sample2]))+0.01)
  
  plotScatter <- ggplot(tableVST, aes(x = log2(tableVST[,sample2]+0.01), y = log2(tableVST[,sample1]+0.01))) +
    geom_point(aes(color = log2FC), alpha = 0.6) +
    xlab(paste0("log2(", sample2, "+0.01)")) +
    ylab(paste0("log2(", sample1, "+0.01)")) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
    xlim(c(minLim, maxLim)) + ylim(c(minLim, maxLim)) +
    theme(legend.position = "none") +
    theme_bw(base_size=15)

  if (peak_type == "Raw"){
    print(paste0(outdir, contrastToApply, "_log2_scatterPlot.pdf"))
    pdf(paste0(outdir, antibody, '_',contrastToApply, "_log2_scatterPlot.pdf"), width = 10, height = 10)
    print(plotScatter)
    dev.off()
    #png for multiqc
    png(paste0(outdir, antibody, '_', contrastToApply, "_log2_scatterPlot_mqc.png"), width=800, height=800)
    print(plotScatter)
    dev.off()
  } else if (peak_type == "Norm"){
    print(paste0(outdir, contrastToApply,'_',normMethod, "_log2_scatterPlot.pdf"))
    pdf(paste0(outdir, antibody, '_',contrastToApply, '_',normMethod, "_log2_scatterPlot_NormPeaks.pdf"), width = 10, height = 10)
    print(plotScatter)
    dev.off()
    #png for multiqc
    png(paste0(outdir, antibody, '_', contrastToApply,'_',normMethod,  "_log2_scatterPlot_NormPeaks_mqc.png"), width=800, height=800)
    print(plotScatter)
    dev.off()
  }
  #------------------------------------------------------------------------------#
  #------------------------------------------------------------------------------#
  
} else {
  # Create the DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = countsTab, colData = colData, design = ~ condition)

  # Perform differential analysis
  stopifnot(colnames(assay(dds)) == colnames(normFactorsList))
  #if (mean(normFactorsList) > 5000){
  #  normFactorsList <- normFactorsList/10000
  #}
  sizeFactors(dds) <- 1/normFactorsList
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  # Get the results
  results <- results(dds, contrast = c("condition", leftContrast, rightContrast))


  #------------------------------------------------------------------------------#
  #------------------------------------------------------------------------------#
  #PCA PLOT
  #------------------------------------------------------------------------------#
  #------------------------------------------------------------------------------#
  cat("Plotting PCA\n")

  pcaPLot <- DESeq2::plotPCA(rlog(dds), intgroup = "condition", returnData = F) +
    theme_bw(base_size = 17) +
    labs(title = "PCA plot", color='Group') +
    ggrepel::geom_text_repel(aes(label=name),show.legend = FALSE )+
    scale_color_manual(values = c("firebrick", "steelblue3"))

  if (peak_type == "Raw"){
    #check if the pca plot already exists
    if (file.exists(paste0(outdir,antibody, '_pcaPlot.pdf'))){
      cat("PCA plot for this antibody already exists, skipping\n")
    } else {     
      pdf(paste0(outdir,antibody, '_pcaPlot.pdf'), width = 10, height = 10)
      print(pcaPLot)
      dev.off()
      #png for multiqc
      png(paste0(outdir,antibody, '_pcaPlot_mqc.png'), width=800, height=800)
      print(pcaPLot)
      dev.off()
    }

  } else if (peak_type == "Norm"){
    if (file.exists(paste0(outdir,antibody, '_pcaPlot_', normMethod,'_NormPeaks.pdf'))){
      cat("PCA plot for this antibody already exists, skipping\n")
    } else {
      pdf(paste0(outdir,antibody, '_pcaPlot_', normMethod,'_NormPeaks.pdf'), width = 10, height = 10)
      print(pcaPLot)
      dev.off()
      #png for multiqc
      png(paste0(outdir,antibody, '_pcaPlot_', normMethod,'_NormPeaks_mqc.png'), width=800, height=800)
      print(pcaPLot)
      dev.off()
    }
  }
  #------------------------------------------------------------------------------#
  #------------------------------------------------------------------------------#
  #VOLCANO PLOT
  #------------------------------------------------------------------------------#
  #------------------------------------------------------------------------------#
  cat("Plotting Volcano\n")
  resultsToPlot <- results %>%  
    as.data.frame() %>%
    mutate(
      Levels = case_when(log2FoldChange >= log2FCcutoff & padj <= padjCutoff ~ paste0("Increased-binding in ", leftContrast),
                        log2FoldChange <= -log2FCcutoff & padj <= padjCutoff ~ paste0("Decreased-binding in ", leftContrast) ,
                        TRUE ~ "Unchanged")
    )

  colors_value <- setNames(
    c("firebrick3", "dodgerblue3", "gray50"),
    c(paste0("Increased-binding in ", leftContrast),
      paste0("Decreased-binding in ", leftContrast),
      "Unchanged")
  )
  p2 <- ggplot(resultsToPlot, aes(log2FoldChange, -log(padj,10))) +
    geom_point(aes(color = Levels)) +
    xlab(expression("log2FC")) + 
    ylab(expression("-log10(p-adjusted)")) +
    scale_color_manual(values = colors_value)+
    theme_light(base_size = 17)+
    labs(title = "Volcano plot")

  if (peak_type == "Raw"){
    pdf(paste0(outdir,antibody, '_', contrastToApply,'_volcanoPlot.pdf'), width = 11, height = 10)
    print(p2)
    dev.off()
    #png for multiqc
    png(paste0(outdir,antibody, '_', contrastToApply,'_volcanoPlot_mqc.png'), width=1000, height=800)
    print(p2)
    dev.off()
  } else if (peak_type == "Norm"){
    pdf(paste0(outdir,antibody, '_', contrastToApply,'_',normMethod, '_volcanoPlot_NormPeaks.pdf'), width = 11, height = 10)
    print(p2)
    dev.off()
    #png for multiqc
    png(paste0(outdir,antibody, '_', contrastToApply,'_', normMethod, '_volcanoPlot_NormPeaks_mqc.png'), width=1000, height=800)
    print(p2)
    dev.off()
  }
  #------------------------------------------------------------------------------#
  #------------------------------------------------------------------------------#
  #------------------------------------------------------------------------------#
  # save table
  cat("Saving results\n")

  results <- results %>%
              data.frame %>% 
              mutate(
                baseMean = round(baseMean, 2),
                log2FoldChange = round(log2FoldChange, 3),
                lfcSE = round(lfcSE, 3),
                stat = round(stat, 3),
              ) %>% 
              rownames_to_column(var = "region")

  write.table(results, file = snakemake@output[['diffTab']], sep = "\t", quote = F, row.names = F, col.names = T)
}

cat("Done\n")
