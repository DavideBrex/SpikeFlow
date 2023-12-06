log <- file(snakemake@log[[1]], open= "wt")
sink(log)
sink(log, type = "message")

require(dplyr)

####################                                                                                                              
# DEFINE FUNCTIONS #                                                                                              
####################

# Function to annotate the peaks using ChIPseeker
Peak_Annot <- function(infile, tssRegion = c(-2500, 2500), TxDb, annoDb) {
  
  # Load packages and set parameters
  require(ChIPseeker)
  require(TxDb, character.only = TRUE)  
  
  # read file and annotate peaks
  df <- readPeakFile(infile)
  #just for the testing of the pipeline
  if(as.character(seqnames(df))[1] == "chr10:14000000-25000000"){
    newname <- c("chr10")
    names(newname) <- "chr10:14000000-25000000"
    df <- renameSeqlevels(df, newname)
  }

  print(df)
  annot <- annotatePeak(df, TxDb = get(TxDb), annoDb = annoDb, tssRegion = tssRegion)
  
  # The program classifies promotres in 1kb, 2kb... this line removes that annotation and leaves just "Promoters"
  annot@anno$annotation <- sub(" \\(.*\\)", "", annot@anno$annotation)
  
  # The program changes the type of object, so go back to GRanges (useful for downstream analysis)
  final <- as.GRanges(annot)
  return(final)
}


#################
## Read inputs ##
#################
input <- as.character(snakemake@input[["peaks"]])
before <- as.numeric(snakemake@params[["before"]])
after <- as.numeric(snakemake@params[["after"]])
out1 <- as.character(snakemake@output[["annotations"]])
out2 <- as.character(snakemake@output[["promoBed"]])
out3 <- as.character(snakemake@output[["distalBed"]])
genome <- as.character(snakemake@params[["genome"]])

#------ If the input is empty the script will stop here ------#
if (file.size(input) == 0){

  file.create(out1)
  file.create(out2)
  file.create(out3)

} else {
  
#------------ Define the genome that is going to be used for the annotation of peaks (to add more genomes they need to be first installed in R) ------------
  if (genome == "mm9") { 
    txdb <- "TxDb.Mmusculus.UCSC.mm9.knownGene"
    annodb <- "org.Mm.eg.db"
  } else if (genome == "mm10") {
    txdb <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
    annodb <- "org.Mm.eg.db"
  } else if  (genome == "hg19") {
    txdb <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
    annodb <- "org.Hs.eg.db"
  } else if  (genome == "hg38") {
    txdb <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    annodb <- "org.Hs.eg.db"
  }

####################
## Annot bed file ##
####################
  annot <- Peak_Annot(input, tssRegion = c(-before, after), TxDb = txdb, annoDb = annodb) %>% as.data.frame()

  # Filter annotared bed file to obtain bed files from peaks overlaping with promoters/distal regions
  distal.peaks <- annot %>% subset(annotation != "Promoter") %>% dplyr::select(c("seqnames", "start", "end", "V4", "V5")) 
  promo.peaks <- annot %>% subset(annotation == "Promoter") %>% dplyr::select(c("seqnames", "start", "end", "V4", "V5")) 
  
  annot <- annot[-c(6,7,8,9,10,11)]

##################
## Write output ##
##################
  write.table(annot, file = out1, sep = "\t", quote = F, row.names = F)
  write.table(promo.peaks, file = out2, sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(distal.peaks, file = out3, sep = "\t", quote = F, row.names = F, col.names = F)
}