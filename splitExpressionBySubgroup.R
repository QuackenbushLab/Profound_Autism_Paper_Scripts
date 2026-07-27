# Read SSC data.
phenoGrp <- NULL
profoundAutismModerateIDOnly <- read.csv(paste0(phenoGrp, "/profoundAutismModerateIDOnly_above8.csv"), row.names = 1)
profoundAutismNonverbalOnly <- read.csv(paste0(phenoGrp, "/profoundAutismNonverbalOnly_above8.csv"), row.names = 1)
profoundAutismBoth <- read.csv(paste0(phenoGrp, "/profoundAutismBoth_above8.csv"), row.names = 1)
verbalMildID <- read.csv(paste0(phenoGrp, "/verbalMildID_above8.csv"), row.names = 1)
verbalNoID <- read.csv(paste0(phenoGrp, "/verbalNoID_above8.csv"), row.names = 1)
verbalGifted <- read.csv(paste0(phenoGrp, "/verbalGifted_above8.csv"), row.names = 1)

# Split genomics data.
sourceDirGenomics <- NULL
genomics <- read.csv(paste0(sourceDirGenomics, "gene_expression_processed/geneExpressionLogCPMLarge.csv"),
                     row.names = 1)
# Split data into families.
familyInfo <- do.call(rbind, lapply(colnames(genomics), function(samp){
  
  # Get the family label. Remove the "X" added by R.
  sampSplitDot <- strsplit(samp, ".", fixed = TRUE)[[1]]
  family <- substr(sampSplitDot[1], 2, nchar(sampSplitDot[1]))
  
  # Get the position within the family.
  positionStr <- substr(sampSplitDot[2], 1, 3)
  
  # Construct a data frame.
  retval <- data.frame(family = family, position = positionStr)
  return(retval)
}))
familyInfoPasted <- paste(familyInfo$family, familyInfo$position, sep = ".")

splitGenomicsData <- function(genomics, familyInfoPasted, group){
  genomicsSubset <- genomics[,which(familyInfoPasted %in% rownames(group))]
  str(genomicsSubset)
  colnames(genomicsSubset) <- familyInfoPasted[which(familyInfoPasted %in% rownames(group))]
  shared <- intersect(colnames(genomicsSubset), rownames(group))
  genomicsSubset <- genomicsSubset[,shared]
  return(genomicsSubset)
}
# Split genomics data.
splitGenomicsProfoundModerateIDOnly <- splitGenomicsData(genomics, familyInfoPasted, profoundAutismModerateIDOnly)
splitGenomicsProfoundNonverbalOnly <- splitGenomicsData(genomics, familyInfoPasted, profoundAutismNonverbalOnly)
splitGenomicsProfoundBoth <- splitGenomicsData(genomics, familyInfoPasted, profoundAutismBoth)
splitGenomicsProfoundEither <- cbind(splitGenomicsProfoundModerateIDOnly, splitGenomicsProfoundNonverbalOnly)
splitGenomicsMildIDVerbal <- splitGenomicsData(genomics, familyInfoPasted, verbalMildID)
splitGenomicsNoIDVerbal <- splitGenomicsData(genomics, familyInfoPasted, verbalNoID)
splitGenomicsGiftedVerbal <- splitGenomicsData(genomics, familyInfoPasted, verbalGifted)
write.csv(splitGenomicsProfoundModerateIDOnly, paste0(sourceDirGenomics, "/splitGenomicsProfoundModerateIDOnly.csv"))
write.csv(splitGenomicsProfoundNonverbalOnly, paste0(sourceDirGenomics, "/splitGenomicsProfoundNonverbalOnly.csv"))
write.csv(splitGenomicsProfoundBoth, paste0(sourceDirGenomics, "/splitGenomicsProfoundBoth.csv"))
write.csv(splitGenomicsProfoundEither, paste0(sourceDirGenomics, "/splitGenomicsProfoundEither.csv"))
write.csv(splitGenomicsMildIDVerbal, paste0(sourceDirGenomics, "/splitGenomicsMildIDVerbal.csv"))
write.csv(splitGenomicsNoIDVerbal, paste0(sourceDirGenomics, "/splitGenomicsNoIDVerbal.csv"))
write.csv(splitGenomicsGiftedVerbal, paste0(sourceDirGenomics, "/splitGenomicsGiftedVerbal.csv"))