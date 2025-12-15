# Run interactively using salloc --mem=30gb --time=10:00:00
library(limma)
library(rtracklayer)

# Find all directories containing expression files.
expressionDirLarge <- "/n/quackenbush_lab/Lab/teicher/RNASeq_Large/SSC/objLinks/seqLib/"

# Loop through all directories and files on the system, read in the gene counts
# as txt, and compile them.
dirsInLarge <- list.dirs(expressionDirLarge, recursive = FALSE)
geneExpressionLarge <- do.call(cbind, lapply(dirsInLarge, function(directory){
  cat(".")
  fileIn <- read.table(paste0(directory, 
                              "/star.p2-wasp_md.featureCounts.txt"),
                       header = TRUE)
  expression <- data.frame(v1 = fileIn[,7])
  dirSplit <- strsplit(directory, "/")[[1]]
  colnames(expression) <- dirSplit[length(dirSplit)]
  rownames(expression) <- fileIn[,1]
  return(expression)
}))
write.csv(geneExpressionLarge, paste0(expressionDirLarge, "geneExpressionRaw.csv"))

# Use VOOM to transform the data.
geneExpressionLogCPMLarge <- as.matrix(voom(geneExpressionLarge))
write.csv(geneExpressionLogCPMLarge, paste0(expressionDirLarge, "geneExpressionLogCPM.csv"))