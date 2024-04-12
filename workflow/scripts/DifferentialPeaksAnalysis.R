# Load the DESeq2 package
log <- file(snakemake@log[[1]], open= "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages(require(DESeq2))
suppressPackageStartupMessages(require(tidyverse))

# Read the params
contrastToApply <- snakemake@params[["contrast"]]
padjCutoff <- snakemake@params[["padjCutoff"]]
log2FCcutoff <- snakemake@params[["log2FCcutoff"]]
outdir <- snakemake@params[["outdir"]]

# Read the file with raw counts
cat("Reading raw counts\n")
countsTab <- read.table(snakemake@input[['rawReadsOnPeaks']], check.names = F, header = TRUE, sep = "\t")
#remove the chr, star, end columns
countsTab <- countsTab[,-c(2,3,4)]
#set rownames to region
countsTab <- countsTab %>%  column_to_rownames("region")

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

normFactorsFiles <- snakemake@input[['logFile']]

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
# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countsTab, colData = colData, design = ~ condition)

# Perform differential analysis
stopifnot(colnames(assay(dds)) == colnames(normFactorsList))
if (mean(normFactorsList) > 5000){
  normFactorsList <- normFactorsList/10000
}

sizeFactors(dds) <- normFactorsList
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
  theme_bw() 

  
pdf(paste0(outdir, 'pcaPlot.pdf'), width = 10, height = 10)
pcaPLot
dev.off()
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
                       log2FoldChange <= -log2FCcutoff & padj <= padjCutoff ~ paste0("Increased-binding in", rightContrast) ,
                       TRUE ~ "Unchanged")
  )

p2 <- ggplot(resultsToPlot, aes(log2FoldChange, -log(padj,10))) +
  geom_point(aes(color = Levels)) +
  xlab(expression("log2FC")) + 
  ylab(expression("-log10(p-adjusted)")) +
  scale_color_manual(values = c( "firebrick3","dodgerblue3",  "gray50"))+
  theme_minimal()

pdf(paste0(outdir, 'volcanoPlot.pdf'), width = 10, height = 10)
p2
dev.off()
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# save table
cat("Saving results\n")

results <- results %>%
            data.frame %>% 
            round(3) %>% 
            rownames_to_column(var = "region")

write.table(results, file = snakemake@output[['diffTab']], sep = "\t", quote = F, row.names = F, col.names = T)

cat("Done\n")