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

# Function for renaming genomics.
rename <- function(file){
  newColNames <- unlist(lapply(colnames(file), function(family){
    return(strsplit(family, "X")[[1]][2])
  }))
  colnames(file) <- newColNames
  return(file)
}

# Read in genomics data.
splitGenomicsProfoundModerateIDOnly <- rename(read.csv(paste0(sourceDirGenomics, "/splitGenomicsProfoundModerateIDOnlyDiff.csv"), row.names = 1))
splitGenomicsProfoundNonverbalOnly <- rename(read.csv(paste0(sourceDirGenomics, "/splitGenomicsProfoundNonverbalOnlyDiff.csv"), row.names = 1))
splitGenomicsProfoundBoth <- rename(read.csv(paste0(sourceDirGenomics, "/splitGenomicsProfoundBothDiff.csv"), row.names = 1))
splitGenomicsMildIDVerbal <- rename(read.csv(paste0(sourceDirGenomics, "/splitGenomicsMildIDVerbalDiff.csv"), row.names = 1))
splitGenomicsNoIDVerbal <- rename(read.csv(paste0(sourceDirGenomics, "/splitGenomicsNoIDVerbalDiff.csv"), row.names = 1))
splitGenomicsGiftedVerbal <- rename(read.csv(paste0(sourceDirGenomics, "/splitGenomicsGiftedVerbalDiff.csv"), row.names = 1))

# Bin age.
binAge <- function(data){
  age <- data$age_at_ados / 12
  age8_10 <- age10_12 <- age12_14 <- age14_16 <- age16_18 <- rep(0, length(age))
  age8_10[intersect(which(age >= 8), which(age < 10))] <- 1
  age10_12[intersect(which(age >= 10), which(age < 12))] <- 1
  age12_14[intersect(which(age >= 12), which(age < 14))] <- 1
  age14_16[intersect(which(age >= 14), which(age < 16))] <- 1
  age16_18[intersect(which(age >= 16), which(age <= 18))] <- 1
  addition <- data.frame(age8_10 = age8_10, age10_12 = age10_12, age12_14 = age12_14,
                         age14_16 = age14_16, age16_18 = age16_18)
  for(ageBracket in colnames(addition)){
    addition[,ageBracket] <- as.factor(addition[,ageBracket])
  }
  data <- cbind(data, addition)
  return(data)
}
profoundAutismModerateIDOnly <- binAge(profoundAutismModerateIDOnly)
profoundAutismNonverbalOnly <- binAge(profoundAutismNonverbalOnly)
profoundAutismBoth <- binAge(profoundAutismBoth)
profoundAutismEither <- rbind(profoundAutismModerateIDOnly, profoundAutismNonverbalOnly)
verbalMildID <- binAge(verbalMildID)
verbalNoID <- binAge(verbalNoID)
verbalGifted <- binAge(verbalGifted)

# Subset SSC data.
# We do not adjust for race or ethnicity because we are comparing against siblings.
# We do adjust for sex of sibling and proband.
siblingDir <- NULL
siblingData <- read.csv(paste0(siblingDir, "/ssc_core_descriptive.csv"),
                        row.names = 1)
rownames(siblingData) <- unlist(lapply(rownames(siblingData), function(row){
  return(paste0(strsplit(row, ".s1")[[1]][1], ".p1"))
}))
covariates <- c("sexCombination")
subsetData <- function(dataSSC, siblingData, g, subtypeName){
  subsetSSC <- dataSSC
  subsetSSC$siblingSex <- siblingData[rownames(subsetSSC), "sex"]
  subsetSSC$sexCombination <- paste(subsetSSC$sex, subsetSSC$siblingSex, sep = "_")
  subsetSSC <- subsetSSC[colnames(g),]
  str(subsetSSC)
  return(subsetSSC)
}
profoundAutismModerateIDOnlySubsetSSC <- subsetData(profoundAutismModerateIDOnly, siblingData, splitGenomicsProfoundModerateIDOnly, "profoundModerateIDOnly")
profoundAutismNonverbalOnlySubsetSSC <- subsetData(profoundAutismNonverbalOnly, siblingData, splitGenomicsProfoundNonverbalOnly, "profoundNonverbalOnly")
profoundBothSubsetSSC <- subsetData(profoundAutismBoth, siblingData, splitGenomicsProfoundBoth, "profoundBoth")
verbalMildIDSubsetSSC <- subsetData(verbalMildID, siblingData, splitGenomicsMildIDVerbal, "mildIDVerbal")
verbalNoIDSubsetSSC <- subsetData(verbalNoID, siblingData, splitGenomicsNoIDVerbal, "noIDVerbal")
verbalGiftedSubsetSSC <- subsetData(verbalGifted, siblingData, splitGenomicsGiftedVerbal, "giftedVerbal")

# Formula
formulaAll <- "gene ~ subtype + sexCombination + age8_10 + age10_12 + age12_14 + age14_16 + age16_18"

# Run linear models.
dir.create(outDirFinal)
runLinearModels <- function(sscGroup1, sscGroup2, genomicsGroup1, genomicsGroup2, 
                            subtype1, subtype2, fileName){
    formula <- formulaAll
    print(paste(subtype1, subtype2))
    print(length(intersect(rownames(sscGroup1), colnames(genomicsGroup1))))
    print(length(intersect(rownames(sscGroup2), colnames(genomicsGroup2))))
    
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
      model <- lm(formula = formula, data = fullDataSet)
      toreturn <- as.data.frame(t(data.frame(model[["coefficients"]])))
      toreturn$gene <- gene
      toreturn$pval <- summary(model)$coefficients[2,4]
      toreturn$stdError <- summary(model)$coefficients[2,2]
      toreturn$rsq <- summary(model)$r.squared
      # Beta value and confidence intervals
      parmName <- names(coef(model))[which(grepl("subtype", names(coef(model))) == TRUE)]
      ci <- confint(model, parmName, level = 0.95)
      toreturn$ciLow <- ci[1]
      toreturn$ciHigh <- ci[2]
      return(toreturn)
    })
    pvalues <- do.call(rbind, pvaluesList)
    pvalues$padj <- stats::p.adjust(pvalues$pval, method = "fdr")
    str(pvalues)
    str(pvalues[which(pvalues$padj < 0.05), "gene"])
    write.csv(pvalues, fileName)
}

# Compare the profound autism groups to the other groups.
runLinearModels(sscGroup1 = profoundAutismModerateIDOnlySubsetSSC, sscGroup2 = verbalMildIDSubsetSSC, 
                genomicsGroup1 = splitGenomicsProfoundModerateIDOnly, genomicsGroup2 = splitGenomicsMildIDVerbal, 
                subtype1 = "profoundModerateIDOnly", subtype2 = "mildIDVerbal",
                fileName = paste0(outDirFinal, "profoundModerateIDOnly_MildIDVerbal.csv"))
runLinearModels(profoundAutismModerateIDOnlySubsetSSC, verbalNoIDSubsetSSC, splitGenomicsProfoundModerateIDOnly,
                splitGenomicsNoIDVerbal, "profoundModerateIDOnly", "noIDVerbal",
                paste0(outDirFinal, "profoundModerateIDOnly_NoIDVerbal.csv"))
runLinearModels(profoundAutismModerateIDOnlySubsetSSC, verbalGiftedSubsetSSC, splitGenomicsProfoundModerateIDOnly,
                splitGenomicsGiftedVerbal, "profoundModerateIDOnly", "giftedVerbal",
                paste0(outDirFinal, "profoundModerateIDOnly_GiftedVerbal.csv"))
runLinearModels(sscGroup1 = profoundAutismNonverbalOnlySubsetSSC, sscGroup2 = verbalMildIDSubsetSSC, 
                genomicsGroup1 = splitGenomicsProfoundNonverbalOnly, genomicsGroup2 = splitGenomicsMildIDVerbal, 
                subtype1 = "profoundNonverbalOnly", subtype2 = "mildIDVerbal",
                fileName = paste0(outDirFinal, "profoundNonverbalOnly_MildIDVerbal.csv"))
runLinearModels(profoundAutismNonverbalOnlySubsetSSC, verbalNoIDSubsetSSC, splitGenomicsProfoundNonverbalOnly,
                splitGenomicsNoIDVerbal, "profoundNonverbalOnly", "noIDVerbal",
                paste0(outDirFinal, "profoundNonverbalOnly_NoIDVerbal.csv"))
runLinearModels(profoundAutismNonverbalOnlySubsetSSC, verbalGiftedSubsetSSC, splitGenomicsProfoundNonverbalOnly,
                splitGenomicsGiftedVerbal, "profoundNonverbalOnly", "giftedVerbal",
                paste0(outDirFinal, "profoundNonverbalOnly_GiftedVerbal.csv"))
runLinearModels(sscGroup1 = profoundBothSubsetSSC, sscGroup2 = verbalMildIDSubsetSSC, 
                genomicsGroup1 = splitGenomicsProfoundBoth, genomicsGroup2 = splitGenomicsMildIDVerbal, 
                subtype1 = "profoundBoth", subtype2 = "mildIDVerbal",
                fileName = paste0(outDirFinal, "profoundBoth_MildIDVerbal.csv"))
runLinearModels(profoundBothSubsetSSC, verbalNoIDSubsetSSC, splitGenomicsProfoundBoth,
                splitGenomicsNoIDVerbal, "profoundBoth", "noIDVerbal",
                paste0(outDirFinal, "profoundBoth_NoIDVerbal.csv"))
runLinearModels(profoundBothSubsetSSC, verbalGiftedSubsetSSC, splitGenomicsProfoundBoth,
                splitGenomicsGiftedVerbal, "profoundBoth", "giftedVerbal",
                paste0(outDirFinal, "profoundBoth_GiftedVerbal.csv"))

# Compare the profound autism groups to each other.
runLinearModels(sscGroup1 = profoundAutismModerateIDOnlySubsetSSC, sscGroup2 = profoundAutismNonverbalOnlySubsetSSC, 
                genomicsGroup1 = splitGenomicsProfoundModerateIDOnly, genomicsGroup2 = splitGenomicsProfoundNonverbalOnly, 
                subtype1 = "profoundModerateIDOnly", subtype2 = "profoundNonverbalOnly",
                fileName = paste0(outDirFinal, "profoundModerateIDOnly_NonverbalOnly.csv"))
runLinearModels(sscGroup1 = profoundAutismModerateIDOnlySubsetSSC, sscGroup2 = profoundBothSubsetSSC, 
                genomicsGroup1 = splitGenomicsProfoundModerateIDOnly, genomicsGroup2 = splitGenomicsProfoundBoth, 
                subtype1 = "profoundModerateIDOnly", subtype2 = "profoundBoth",
                fileName = paste0(outDirFinal, "profoundModerateIDOnly_ProfoundBoth.csv"))
runLinearModels(sscGroup1 = profoundAutismNonverbalOnlySubsetSSC, sscGroup2 = profoundBothSubsetSSC, 
                genomicsGroup1 = splitGenomicsProfoundNonverbalOnly, genomicsGroup2 = splitGenomicsProfoundBoth, 
                subtype1 = "profoundNonverbalOnly", subtype2 = "profoundBoth",
                fileName = paste0(outDirFinal, "profoundNonverbalOnly_ProfoundBoth.csv"))