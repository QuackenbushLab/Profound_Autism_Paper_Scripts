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
genomics <- read.csv(paste0(sourceDirGenomics, "diffExpression.csv"),
                     row.names = 1)
newColNames <- unlist(lapply(colnames(genomics), function(family){
  return(paste(strsplit(family, "X")[[1]][2], "p1", sep = "."))
}))
colnames(genomics) <- newColNames

splitGenomicsData <- function(genomics, group){
  shared <- intersect(colnames(genomics), rownames(group))
  print(length(shared))
  genomics <- genomics[,shared]
  return(genomics)
}
# Split genomics data.
splitGenomicsProfoundModerateIDOnly <- splitGenomicsData(genomics, profoundAutismModerateIDOnly)
splitGenomicsProfoundNonverbalOnly <- splitGenomicsData(genomics, profoundAutismNonverbalOnly)
splitGenomicsProfoundBoth <- splitGenomicsData(genomics, profoundAutismBoth)
splitGenomicsProfoundEither <- cbind(splitGenomicsProfoundModerateIDOnly, splitGenomicsProfoundNonverbalOnly)
splitGenomicsMildIDVerbal <- splitGenomicsData(genomics, verbalMildID)
splitGenomicsNoIDVerbal <- splitGenomicsData(genomics, verbalNoID)
splitGenomicsGiftedVerbal <- splitGenomicsData(genomics, verbalGifted)
write.csv(splitGenomicsProfoundModerateIDOnly, paste0(sourceDirGenomics, "/splitGenomicsProfoundModerateIDOnlyDiff.csv"))
write.csv(splitGenomicsProfoundNonverbalOnly, paste0(sourceDirGenomics, "/splitGenomicsProfoundNonverbalOnlyDiff.csv"))
write.csv(splitGenomicsProfoundBoth, paste0(sourceDirGenomics, "/splitGenomicsProfoundBothDiff.csv"))
write.csv(splitGenomicsProfoundEither, paste0(sourceDirGenomics, "/splitGenomicsProfoundEitherDiff.csv"))
write.csv(splitGenomicsMildIDVerbal, paste0(sourceDirGenomics, "/splitGenomicsMildIDVerbalDiff.csv"))
write.csv(splitGenomicsNoIDVerbal, paste0(sourceDirGenomics, "/splitGenomicsNoIDVerbalDiff.csv"))
write.csv(splitGenomicsGiftedVerbal, paste0(sourceDirGenomics, "/splitGenomicsGiftedVerbalDiff.csv"))