# Read all files containing sex, sibling sex, and age info.
clinFiles <- list.files("../profoundAutism/")
clinData <- do.call(rbind, lapply(clinFiles, function(f){
  return(read.csv(paste0("../profoundAutism/", f), header = TRUE, row.names =1))
}))
siblingData <- read.csv("../profoundAutism/sibling.csv",
                        row.names = 1)
rownames(siblingData) <- unlist(lapply(rownames(siblingData), function(row){
  return(paste0(strsplit(row, ".s1")[[1]][1], ".p1"))
}))

# Read current covariate file.
covar <- read.table("../PLINK/Omni_covar_2.txt", sep = " ")
colnames(covar) <- c("FID", "IID", paste0("PC", 1:6))

# Map the samples.
map <- read.table("../eQTL/merged.txt", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, sep = "\t")
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

# Subset to those that have sibling data.
whichHasSib <- which(rownames(clinData) %in% rownames(siblingData))
clinDataSib <- clinData[whichHasSib,]
clinDataSib$siblingSex <- siblingData[rownames(clinDataSib),"sex"]
clinDataSampsSib <- clinDataSamps[whichHasSib]

# Add columns for the sex and age and label the columns.
clinDataCovar <- clinDataSib[which(clinDataSampsSib %in% covar$IID),]
rownames(clinDataCovar) <- clinDataSampsSib[which(clinDataSampsSib %in% covar$IID)]
clinDataCovar$age_at_ados <- clinDataCovar$age_at_ados / 12
rownames(covar) <- covar$IID
covarNew <- cbind(covar[rownames(clinDataCovar),], clinDataCovar[,c("sex", "age_at_ados", "siblingSex")])
covarNew <- covarNew[,c(1,2,9,10,11)]

# Convert sex to numeric (female = 0, male = 1)
# Need to do this for all combinations.
covarNew[which(covarNew$sex == "female"), "sex"] <- "0"
covarNew[which(covarNew$sex == "male"), "sex"] <- "1"
covarNew[which(covarNew$siblingSex == "female"), "siblingSex"] <- "0"
covarNew[which(covarNew$siblingSex == "male"), "siblingSex"] <- "1"
covarNew$sex <- as.numeric(covarNew$sex)
covarNew$siblingSex <- as.numeric(covarNew$siblingSex)

# Filter to include only the samples also in the BED file.
bed <- read.table("../eQTL/diffExpression.pheno.profoundBoth.bed", sep = "\t",
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
write.table(covarNew, "../eQTL/covarNewProfoundBothDiff.txt", sep = "\t", quote = FALSE,
            col.names = FALSE)
write.table(bed, "../eQTL/diffExpression.pheno.profoundBoth.filt.bed", sep = "\t", quote = FALSE, row.names = FALSE)