# Read all files containing sex and age info.
clinDir <- "../profoundAutism/"
clinFiles <- list("profoundAutismBoth_above8.csv",
                  "profoundAutismModerateIDOnly_above8.csv",
                  "verbalNoID_above8.csv",
                  "verbalMildID_above8.csv",
                  "verbalGifted_above8.csv")
clinData <- do.call(rbind, lapply(clinFiles, function(f){
  return(read.csv(paste0(clinDir, f), header = TRUE, row.names =1))
}))

# Read current covariate file.
covar <- read.table("../PLINK/Omni_covar_2.txt", sep = " ")
colnames(covar) <- c("FID", "IID", paste0("PC", 1:6))

# Map the samples.
mergedFile <- "../eQTL/merged.txt"
map <- read.table(mergedFile, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, sep = "\t")
clinDataSamps <- unlist(lapply(rownames(clinData), function(iid){
  val <- NA
  if(length(which(map$individual_id == iid)) > 0){
    val <- map[which(map$individual_id == iid), "array_id"]
    if(length(val) > 1)
      if(length(which(val %in% covar$IID) > 0)){
        val <- val[which(val %in% covar$IID)[1]]
      }else{
        val <- val[1]
      }
  }
  if(length(val) == 0){
    print(which(map$individual_id == iid))
  }
  return(val)
}))

# Add columns for the sex and age and label the columns.
clinDataCovar <- clinData[which(clinDataSamps %in% covar$IID),]
rownames(clinDataCovar) <- clinDataSamps[which(clinDataSamps %in% covar$IID)]
clinDataCovar$age_at_ados <- clinDataCovar$age_at_ados / 12
rownames(covar) <- covar$IID
covarNew <- cbind(covar[rownames(clinDataCovar),], clinDataCovar[,c("sex", "age_at_ados")])

# Convert sex to numeric (female = 0, male = 1)
covarNew[which(covarNew$sex == "female"), "sex"] <- "0"
covarNew[which(covarNew$sex == "male"), "sex"] <- "1"
covarNew$sex <- as.numeric(covarNew$sex)

formatCovar <- function(expressionBed, covarFile, filtBedFile){
  # Filter to include only the samples also in the BED file.
  bed <- read.table(expressionBed, sep = "\t",
                    comment.char = "", header = TRUE, check.names = FALSE)
  str(bed)
  shared <- intersect(colnames(bed), covarNew$IID)
  rownames(covarNew) <- covarNew$IID
  covarNew <- covarNew[shared,]
  bed <- bed[,c("#chr", "start", "end", "phenotype_id", shared)]
  
  # Modify chromosome names.
  bed[,"#chr"] <- unlist(lapply(bed[,"#chr"], function(line){
    return(strsplit(line, split = "chr")[[1]][2])
  }))
  str(covarNew)
  str(bed)
  str(shared)
  # Remove the family ID.
  covarNew <- covarNew[,2:ncol(covarNew)]
  colnames(covarNew)[1] <- "ID"
  
  # Transpose and write.
  covarNew <- t(covarNew)
  write.table(covarNew, covarFile, sep = "\t", quote = FALSE,
              col.names = FALSE)
  write.table(bed, filtBedFile, sep = "\t", quote = FALSE, row.names = FALSE)
}
formatCovar(expressionBed = "../eQTL/expression.pheno.profoundBoth.bed",
            covarFile = "../eQTL/covarNewProfoundBoth.txt",
            filtBedFile = "../eQTL/expression.pheno.profoundBoth.filt.bed")
formatCovar(expressionBed = "../eQTL/expression.pheno.profoundModerateIDOnly.bed",
            covarFile = "../eQTL/covarNewProfoundModerateID.txt",
            filtBedFile = "../eQTL/expression.pheno.profoundModerateID.filt.bed")
formatCovar(expressionBed = "../eQTL/expression.pheno.verbalMildID.bed",
            covarFile = "../eQTL/covarNewVerbalMildID.txt",
            filtBedFile = "../eQTL/expression.pheno.verbalMildID.filt.bed")
formatCovar(expressionBed = "../eQTL/expression.pheno.verbalNoID.bed",
            covarFile = "../eQTL/covarNewVerbalNoID.txt",
            filtBedFile = "../eQTL/expression.pheno.verbalNoID.filt.bed")
formatCovar(expressionBed = "../eQTL/expression.pheno.verbalGifted.bed",
            covarFile = "../eQTL/covarNewVerbalGifted.txt",
            filtBedFile = "../eQTL/expression.pheno.verbalGifted.filt.bed")
