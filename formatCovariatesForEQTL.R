# Read all files containing sex and age info.
clinDir <- "../profoundAutism/"
clinFiles <- list.files(clinDir)
clinData <- do.call(rbind, lapply(clinFiles, function(f){
  return(read.csv(paste0(clinDir, f), header = TRUE, row.names =1))
}))

# Read current covariate file.
covar <- read.table("../PLINK/Omni_covar_2.txt", sep = " ")
colnames(covar) <- c("FID", "IID", paste0("PC", 1:6))

# Map the samples.
mergedFile <- "../eQTL/merged.txt"
expressionBed <- "../eQTL/expression.pheno.profoundBoth.bed"
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

# Filter to include only the samples also in the BED file.
bed <- read.table(expressionBed, sep = "\t",
                  comment.char = "", header = TRUE, check.names = FALSE)
shared <- intersect(colnames(bed), covarNew$IID)
rownames(covarNew) <- covarNew$IID
covarNew <- covarNew[shared,]
bed <- bed[,c("#chr", "start", "end", "phenotype_id", shared)]

# Remove the family ID.
covarNew <- covarNew[,2:ncol(covarNew)]
colnames(covarNew)[1] <- "ID"

# Transpose and write.
covarNew <- t(covarNew)
write.table(covarNew, "../eQTL/covarNewProfoundBoth.txt", sep = "\t", quote = FALSE,
            col.names = FALSE)
write.table(bed, "../eQTL/expression.pheno.profoundBoth.filt.bed", sep = "\t", quote = FALSE, row.names = FALSE)