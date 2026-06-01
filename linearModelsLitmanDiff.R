# Read data.
sourceDirGenomics <- NULL
phenoGrp <- NULL
outDir <- NULL
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
siblingData <- read.csv("/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/SFARI/SSC\ Version\ 15.3\ Phenotype\ Dataset/Designated\ Unaffected\ Sibling\ Data/ssc_core_descriptive.csv",
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
# Function for renaming genomics.
rename <- function(file){
  newColNames <- unlist(lapply(colnames(file), function(family){
    return(strsplit(family, "X")[[1]][2])
  }))
  colnames(file) <- newColNames
  return(file)
}


# Formula
formulaAll <- "gene ~ subtype + sexCombination + age8_10 + age10_12 + age12_14 + age14_16 + age16_18"

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