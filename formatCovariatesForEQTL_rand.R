# Read all files containing sex and age info.
clinDir <- "../profoundAutism/"
clinData <- read.csv(paste0(clinDir,"profoundAutismBoth_above8.csv"), header = TRUE, row.names =1)

# Read current covariate file.
covar <- read.table("../PLINK/Omni_covar_2.txt", sep = " ")
colnames(covar) <- c("FID", "IID", paste0("PC", 1:6))

# Map the samples.
mergedFile <- "../eQTL/merged.txt"
map <- read.table(mergedFile, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, sep = "\t")

# Clinical samples.
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
clinDataCovar <- clinData[which(clinDataSamps %in% covar$IID),]
rownames(clinDataCovar) <- clinDataSamps[which(clinDataSamps %in% covar$IID)]
clinDataCovar$age_at_ados <- clinDataCovar$age_at_ados / 12
rownames(covar) <- covar$IID
covarNew <- cbind(covar[rownames(clinDataCovar),], clinDataCovar[,c("sex", "age_at_ados")])

# Convert sex to numeric (female = 0, male = 1)
covarNew[which(covarNew$sex == "female"), "sex"] <- "0"
covarNew[which(covarNew$sex == "male"), "sex"] <- "1"
covarNew$sex <- as.numeric(covarNew$sex)
str(covarNew)

for(i in 1:10){
  for(j in 1:2){
    expressionBed <- paste0("../eQTL/expression.pheno.profoundBoth_", i, "_", j, ".bed")
    snps <- read.table(paste0("../eQTL/profoundAutismBoth_rand_", i, "_", j, ".psam"),
                       sep = "\t", header = TRUE)
    
    # Filter to include only the samples also in the BED file.
    bed <- read.table(expressionBed, sep = "\t",
                      comment.char = "", header = TRUE, check.names = FALSE)
    shared <- Reduce(intersect, list(colnames(bed), covarNew$IID, snps$IID))
    covarNewLoc <- covarNew
    rownames(covarNewLoc) <- covarNewLoc$IID
    covarNewLoc <- covarNewLoc[shared,]
    bed <- bed[,c("#chr", "start", "end", "phenotype_id", shared)]
    
    # Remove the family ID.
    covarNewLoc <- covarNewLoc[,2:ncol(covarNewLoc)]
    colnames(covarNewLoc)[1] <- "ID"
    
    # Fix the chromosome formatting.
    bed[,"#chr"] <- unlist(lapply(bed[,"#chr"], function(chr){
      return(strsplit(chr, "chr")[[1]][2])
    }))
    print(unique(bed[,"#chr"]))
    str(bed)
    
    # Transpose and write.
    covarNewLoc <- t(covarNewLoc)
    write.table(covarNewLoc, paste0("../eQTL/covarNewProfoundBoth_", i, "_", j, ".txt"), sep = "\t",
                quote = FALSE, col.names = FALSE)
    write.table(bed, paste0("../eQTL/expression.pheno.profoundBoth_", i, "_", j, "_filt.bed"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}