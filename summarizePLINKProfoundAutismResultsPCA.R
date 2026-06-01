# Read in PLINK results with covariates.
sourceDir <- NULL

# Find significant SNPs.
findSignificant <- function(fileName){
  res <- read.table(fileName, header = TRUE, comment.char = "")
  resSNP <- res[which(res$V11 == "ADD"),]
  resSNP$padj <- p.adjust(resSNP$V16, method = "fdr")
  print(min(resSNP$padj))
  sig <- resSNP[which(resSNP$padj < 0.05),]
  str(sig)
  write.csv(sig, paste0(fileName, ".sig"))
  return(sig)
}

# Check the "no sex" models.
profoundAutismAllGWAS <- findSignificant(paste0(sourceDir, "profoundAll_noSex_pcs.PHENO1.glm.logistic.hybrid"))
profoundAutismBothGWAS <- findSignificant(paste0(sourceDir, "profoundBoth_noSex_pcs.PHENO1.glm.logistic.hybrid"))
profoundAutismNonverbalGWAS <- findSignificant(paste0(sourceDir, "profoundNonverbalOnly_noSex_pcs.assoc.logistic")) #Variance inflation factor on COVAR1 is too high.
profoundAutismIDGWAS <- findSignificant(paste0(sourceDir, "profoundModerateIDOnly_noSex_pcs.assoc.logistic")) #Variance inflation factor on COVAR1 is too high.
profoundAutismBothIDGWAS <- findSignificant(paste0(sourceDir, "profoundBothID_noSex_pcs.PHENO1.glm.logistic.hybrid"))
profoundAutismBothNonverbalGWAS <- findSignificant(paste0(sourceDir, "profoundBothNonverbal_noSex_pcs.PHENO1.glm.logistic.hybrid"))
verbalAutismAllGWAS <- findSignificant(paste0(sourceDir, "notProfoundAll_noSex_pcs.PHENO1.glm.logistic.hybrid"))
verbalAutismNotGiftedGWAS <- findSignificant(paste0(sourceDir, "notProfoundNoGifted_noSex_pcs.PHENO1.glm.logistic.hybrid"))
verbalAutismNoIDGWAS <- findSignificant(paste0(sourceDir, "notProfoundNoID_noSex_pcs.PHENO1.glm.logistic.hybrid"))
verbalMildIDGWAS <- findSignificant(paste0(sourceDir, "verbalMildID_noSex_pcs.PHENO1.glm.logistic.hybrid"))
verbalNoIDGWAS <- findSignificant(paste0(sourceDir, "verbalNoID_noSex_pcs.PHENO1.glm.logistic.hybrid"))
verbalGiftedGWAS <- findSignificant(paste0(sourceDir, "verbalGifted_noSex_pcs.assoc.logistic")) #Variance inflation factor on COVAR1 is too high.


# Check the male-only models.
profoundAutismAllGWAS_M <- findSignificant(paste0(sourceDir, "profoundAll_male_pcs.PHENO1.glm.logistic.hybrid"))
verbalAutismAllGWAS_M <- findSignificant(paste0(sourceDir, "notProfoundAll_male_pcs.PHENO1.glm.logistic.hybrid"))

# Check the female-only models.
profoundAutismAllGWAS_F <- findSignificant(paste0(sourceDir, "profoundAll_female_pcs.PHENO1.glm.logistic.hybrid")) 
verbalAutismAllGWAS_F <- findSignificant(paste0(sourceDir, "notProfoundAll_female_pcs.assoc.logistic")) #Variance inflation factor on COVAR2 is too high.

# Find where these genes rank in each of the lists. Also, what about PA genes?
snps <- verbalAutismNotGiftedGWAS$ID
getRanking <- function(res, resFile, sig){
  
  # Read and adjust.
  str(res)
  resSNP <- res[which(res$TEST == "ADD"),]
  str(resSNP)
  resSNP <- resSNP[order(resSNP$P),]
  
  # Get the rankings for everything not shared.
  ranks <- unlist(lapply(sig, function(snp){
    whichSNP <- which(resSNP$ID == snp)
    ranking <- whichSNP / nrow(resSNP)
    return(ranking)
  }))
  names(ranks) <- sig
  
  # Write the result.
  write.csv(ranks, paste0(resFile, ".ranking"))
  print(range(ranks))
  return(ranks)
}
sigASD <- read.csv(paste0(sourceDir, "notProfoundAll_noSex_pcs.PHENO1.glm.logistic.hybrid.sig"))$ID
sigPA <- read.csv(paste0(sourceDir, "profoundAll_noSex_pcs.PHENO1.glm.logistic.hybrid.sig"))$ID
res <- read.table(paste0(sourceDir, "profoundAll_noSex_pcs.PHENO1.glm.logistic.hybrid"), header = TRUE, comment.char = "", sep = "\t")
profoundAutismAllRanks <- getRanking(res, paste0(sourceDir, "profoundAll_noSex_pcs.PHENO1.glm.logistic.hybrid"), setdiff(sigASD, sigPA))