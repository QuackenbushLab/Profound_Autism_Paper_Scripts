# Read data.
sourceDirGenomics <- NULL
phenoGrp <- NULL
outDirFinal <- NULL
profoundAutismModerateIDOnly <- read.csv(paste0(phenoGrp, "/profoundAutismModerateIDOnly_above8.csv"), row.names = 1)
profoundAutismNonverbalOnly <- read.csv(paste0(phenoGrp, "/profoundAutismNonverbalOnly_above8.csv"), row.names = 1)
profoundAutismBoth <- read.csv(paste0(phenoGrp, "/profoundAutismBoth_above8.csv"), row.names = 1)
verbalMildID <- read.csv(paste0(phenoGrp, "/verbalMildID_above8.csv"), row.names = 1)
verbalNoID <- read.csv(paste0(phenoGrp, "/verbalNoID_above8.csv"), row.names = 1)
verbalGifted <- read.csv(paste0(phenoGrp, "/verbalGifted_above8.csv"), row.names = 1)

# Read in genomics data.
splitGenomicsProfoundModerateIDOnly <- read.csv(paste0(sourceDirGenomics, "/splitGenomicsProfoundModerateIDOnly.csv"), row.names = 1)
splitGenomicsProfoundNonverbalOnly <- read.csv(paste0(sourceDirGenomics, "/splitGenomicsProfoundNonverbalOnly.csv"), row.names = 1)
splitGenomicsProfoundBoth <- read.csv(paste0(sourceDirGenomics, "/splitGenomicsProfoundBoth.csv"), row.names = 1)
splitGenomicsMildIDVerbal <- read.csv(paste0(sourceDirGenomics, "/splitGenomicsMildIDVerbal.csv"), row.names = 1)
splitGenomicsNoIDVerbal <- read.csv(paste0(sourceDirGenomics, "/splitGenomicsNoIDVerbal.csv"), row.names = 1)
splitGenomicsGiftedVerbal <- read.csv(paste0(sourceDirGenomics, "/splitGenomicsGiftedVerbal.csv"), row.names = 1)

# Formula
formulaAll <- "gene ~ subtype + sex + race + ethnicity + age8_10 + age10_12 + age12_14 + age14_16 + age16_18"

# Run linear models.
dir.create(outDirFinal)
# Do the comparisons for IQ-related, speech-related, and intersectional relationships.
runLinearModels <- function(sscGroup1, sscGroup2, genomicsGroup1, genomicsGroup2, 
                            subtype1, subtype2, fileName){
  formula <- formulaAll
  ssc <- rbind(sscGroup1, sscGroup2)
  gen <- rbind(t(genomicsGroup1), t(genomicsGroup2))
  ssc$subtype <- c(rep(subtype1, nrow(sscGroup1)), 
                   rep(subtype2, nrow(sscGroup2)))
  shared <- intersect(rownames(ssc), rownames(gen))
  shared <- Reduce(intersect, list(shared, rownames(ssc)[which(ssc$ethnicity != "")],
                                   rownames(ssc)[which(ssc$race != "not-specified")]))
  ssc <- ssc[shared,]
  gen <- gen[shared,]
  pvaluesList <- lapply(colnames(gen), function(gene){
    fullDataSet <- ssc
    fullDataSet$gene <- gen[,gene]
    toreturn <- NULL
    tryCatch({
      model <- lm(formula = formula, data = fullDataSet)
      toreturn <- as.data.frame(t(data.frame(model[["coefficients"]])))
      toreturn$gene <- gene
      toreturn$pval <- summary(model)$coefficients[2,4]
      toreturn$stdError <- summary(model)$coefficients[2,2]
      toreturn$rsq <- summary(model)$r.squared
      co <- summary(model)$coefficients
      term <- grep("subtypeB", rownames(co), value = TRUE)[1]
      tval <- unname(co[term, "t value"])
      toreturn$tval <- tval
    }, error = function(cond){print(cond)})
    return(toreturn)
  })
  pvalues <- do.call(rbind, pvaluesList)
  pvalues$padj <- stats::p.adjust(pvalues$pval, method = "fdr")
  str(pvalues)
  write.csv(pvalues, fileName)
}
runLinearModels(sscGroup1 = profoundAutismNonverbalOnlySubsetSSC, sscGroup2 = profoundBothSubsetSSC, 
                genomicsGroup1 = splitGenomicsProfoundNonverbalOnly, genomicsGroup2 = splitGenomicsProfoundBoth, 
                subtype1 = "A_profoundNonverbal", subtype2 = "B_profoundBoth",
                fileName = paste0(outDirFinal, "IQ1_nonverbalToBoth.csv"))
runLinearModels(sscGroup1 = verbalMildIDSubsetSSC, sscGroup2 = profoundAutismModerateIDOnlySubsetSSC, 
                genomicsGroup1 = splitGenomicsMildIDVerbal, genomicsGroup2 = splitGenomicsProfoundModerateIDOnly, 
                subtype1 = "A_mildID", subtype2 = "B_profoundModerateID",
                fileName = paste0(outDirFinal, "IQ2_mildToModerate.csv"))
runLinearModels(sscGroup1 = verbalNoIDSubsetSSC, sscGroup2 = verbalMildIDSubsetSSC, 
                genomicsGroup1 = splitGenomicsNoIDVerbal, genomicsGroup2 = splitGenomicsMildIDVerbal, 
                subtype1 = "A_noID", subtype2 = "B_mildID",
                fileName = paste0(outDirFinal, "IQ3_noToMild.csv"))
runLinearModels(sscGroup1 = verbalGiftedSubsetSSC, sscGroup2 = verbalNoIDSubsetSSC, 
                genomicsGroup1 = splitGenomicsGiftedVerbal, genomicsGroup2 = splitGenomicsNoIDVerbal, 
                subtype1 = "A_gifted", subtype2 = "B_noID",
                fileName = paste0(outDirFinal, "IQ5_giftedToNo.csv"))
runLinearModels(sscGroup1 = profoundAutismModerateIDOnlySubsetSSC, sscGroup2 = profoundBothSubsetSSC, 
                genomicsGroup1 = splitGenomicsProfoundModerateIDOnly, genomicsGroup2 = splitGenomicsProfoundBoth, 
                subtype1 = "A_profoundModerateID", subtype2 = "B_profoundBoth",
                fileName = paste0(outDirFinal, "Speech1_moderateToBoth.csv"))
runLinearModels(sscGroup1 = verbalMildIDSubsetSSC, sscGroup2 = profoundAutismNonverbalOnlySubsetSSC, 
                genomicsGroup1 = splitGenomicsMildIDVerbal, genomicsGroup2 = splitGenomicsProfoundNonverbalOnly, 
                subtype1 = "A_mildID", subtype2 = "B_profoundNonverbal",
                fileName = paste0(outDirFinal, "Speech2_mildToNonverbal.csv"))
runLinearModels(sscGroup1 = verbalNoIDSubsetSSC, sscGroup2 = profoundAutismNonverbalOnlySubsetSSC, 
                genomicsGroup1 = splitGenomicsNoIDVerbal, genomicsGroup2 = splitGenomicsProfoundNonverbalOnly, 
                subtype1 = "A_noID", subtype2 = "B_profoundNonverbal",
                fileName = paste0(outDirFinal, "Speech3_noToNonverbal.csv"))
runLinearModels(sscGroup1 = verbalGiftedSubsetSSC, sscGroup2 = profoundAutismNonverbalOnlySubsetSSC, 
                genomicsGroup1 = splitGenomicsGiftedVerbal, genomicsGroup2 = splitGenomicsProfoundNonverbalOnly, 
                subtype1 = "A_gifted", subtype2 = "B_profoundNonverbal",
                fileName = paste0(outDirFinal, "Speech4_giftedToNonverbal.csv"))
runLinearModels(sscGroup1 = verbalMildIDSubsetSSC, sscGroup2 = profoundBothSubsetSSC, 
                genomicsGroup1 = splitGenomicsMildIDVerbal, genomicsGroup2 = splitGenomicsProfoundBoth, 
                subtype1 = "A_mildID", subtype2 = "B_profoundBoth",
                fileName = paste0(outDirFinal, "Intersectional1_mildToBoth.csv"))
runLinearModels(sscGroup1 = verbalNoIDSubsetSSC, sscGroup2 = profoundBothSubsetSSC, 
                genomicsGroup1 = splitGenomicsNoIDVerbal, genomicsGroup2 = splitGenomicsProfoundBoth, 
                subtype1 = "A_noID", subtype2 = "B_profoundBoth",
                fileName = paste0(outDirFinal, "Intersectional2_noToBoth.csv"))
runLinearModels(sscGroup1 = verbalGiftedSubsetSSC, sscGroup2 = profoundBothSubsetSSC, 
                genomicsGroup1 = splitGenomicsGiftedVerbal, genomicsGroup2 = splitGenomicsProfoundBoth, 
                subtype1 = "A_gifted", subtype2 = "B_profoundBoth",
                fileName = paste0(outDirFinal, "Intersectional3_giftedToBoth.csv"))