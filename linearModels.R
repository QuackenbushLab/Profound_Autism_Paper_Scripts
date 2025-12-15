# Read data.
sourceDirGenomics <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/SFARI/large_results/"
phenoGrp <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/SFARI/profoundAutism/phenotypeGroups/"
profoundAutismModerateIDOnly <- read.csv(paste0(phenoGrp, "/profoundAutismModerateIDOnly_above8.csv"), row.names = 1)
profoundAutismNonverbalOnly <- read.csv(paste0(phenoGrp, "/profoundAutismNonverbalOnly_above8.csv"), row.names = 1)
profoundAutismBoth <- read.csv(paste0(phenoGrp, "/profoundAutismBoth_above8.csv"), row.names = 1)
verbalMildID <- read.csv(paste0(phenoGrp, "/verbalMildID_above8.csv"), row.names = 1)
verbalNoID <- read.csv(paste0(phenoGrp, "/verbalNoID_above8.csv"), row.names = 1)
verbalGifted <- read.csv(paste0(phenoGrp, "/verbalGifted_above8.csv"), row.names = 1)
outDir <- NULL

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
outDirFinal <- paste0(outDir, "/geneExpression_ageBinnedModels/")
dir.create(outDirFinal)
newColNamesFunc <- function(genomics){
  newColNames <- unlist(lapply(colnames(genomics), function(family){
    return(paste(strsplit(family, "X")[[1]][2]))
  }))
  colnames(genomics) <- newColNames
  return(genomics)
}
runLinearModels <- function(sscGroup1, sscGroup2, genomicsGroup1, genomicsGroup2, 
                            subtype1, subtype2, fileName){
    formula <- formulaAll
    genomicsGroup1 <- t(newColNamesFunc(genomicsGroup1))
    genomicsGroup2 <- t(newColNamesFunc(genomicsGroup2))
    print(paste(subtype1, subtype2))
    sharedGroup1Init <- intersect(rownames(sscGroup1), rownames(genomicsGroup1))
    sharedGroup2Init <- intersect(rownames(sscGroup2), rownames(genomicsGroup2))
    sharedGroup1 <- Reduce(intersect, list(sharedGroup1Init, rownames(sscGroup1)[which(sscGroup1$ethnicity != "")],
                                     rownames(sscGroup1)[which(sscGroup1$race != "not-specified")]))
    sharedGroup2 <- Reduce(intersect, list(sharedGroup2Init, rownames(sscGroup2)[which(sscGroup2$ethnicity != "")],
                                           rownames(sscGroup2)[which(sscGroup2$race != "not-specified")]))
    print(paste(length(sharedGroup1), length(sharedGroup2)))
    # ssc <- rbind(sscGroup1, sscGroup2)
    # gen <- rbind(t(genomicsGroup1), t(genomicsGroup2))
    # str(gen)
    # ssc$subtype <- c(rep(subtype1, nrow(sscGroup1)), 
    #                  rep(subtype2, nrow(sscGroup2)))
    # shared <- intersect(rownames(ssc), rownames(gen))
    # shared <- Reduce(intersect, list(shared, rownames(ssc)[which(ssc$ethnicity != "")],
    #                                  rownames(ssc)[which(ssc$race != "not-specified")]))
    # ssc <- ssc[shared,]
    # gen <- gen[shared,]
    # print(dim(ssc))
    # pvaluesList <- lapply(colnames(gen), function(gene){
    #   fullDataSet <- ssc
    #   fullDataSet$gene <- gen[,gene]
    #   model <- lm(formula = formula, data = fullDataSet)
    #   toreturn <- as.data.frame(t(data.frame(model[["coefficients"]])))
    #   toreturn$gene <- gene
    #   toreturn$pval <- summary(model)$coefficients[2,4]
    #   toreturn$stdError <- summary(model)$coefficients[2,2]
    #   toreturn$rsq <- summary(model)$r.squared
    #   return(toreturn)
    # })
    # pvalues <- do.call(rbind, pvaluesList)
    # pvalues$padj <- stats::p.adjust(pvalues$pval, method = "fdr")
    # str(pvalues[which(pvalues$padj < 0.05), "gene"])
    # write.csv(pvalues, fileName)
}

# Compare the profound autism groups to the other groups.
runLinearModels(sscGroup1 = profoundAutismModerateIDOnly, sscGroup2 = verbalMildID, 
                genomicsGroup1 = splitGenomicsProfoundModerateIDOnly, genomicsGroup2 = splitGenomicsMildIDVerbal, 
                subtype1 = "profoundModerateIDOnly", subtype2 = "mildIDVerbal",
                fileName = paste0(outDirFinal, "profoundModerateIDOnly_MildIDVerbal.csv"))
runLinearModels(profoundAutismModerateIDOnly, verbalNoID, splitGenomicsProfoundModerateIDOnly,
                splitGenomicsNoIDVerbal, "profoundModerateIDOnly", "noIDVerbal",
                paste0(outDirFinal, "profoundModerateIDOnly_NoIDVerbal.csv"))
runLinearModels(profoundAutismModerateIDOnly, verbalGifted, splitGenomicsProfoundModerateIDOnly,
                splitGenomicsGiftedVerbal, "profoundModerateIDOnly", "giftedVerbal",
                paste0(outDirFinal, "profoundModerateIDOnly_GiftedVerbal.csv"))
runLinearModels(sscGroup1 = profoundAutismNonverbalOnly, sscGroup2 = verbalMildID, 
                genomicsGroup1 = splitGenomicsProfoundNonverbalOnly, genomicsGroup2 = splitGenomicsMildIDVerbal, 
                subtype1 = "profoundAutismNonverbalOnly", subtype2 = "mildIDVerbal",
                fileName = paste0(outDirFinal, "profoundNonverbalOnly_MildIDVerbal.csv"))
runLinearModels(profoundAutismNonverbalOnly, verbalNoID, splitGenomicsProfoundNonverbalOnly,
                splitGenomicsNoIDVerbal, "profoundAutismNonverbalOnly", "noIDVerbal",
                paste0(outDirFinal, "profoundNonverbalOnly_NoIDVerbal.csv"))
runLinearModels(profoundAutismNonverbalOnly, verbalGifted, splitGenomicsProfoundNonverbalOnly,
                splitGenomicsGiftedVerbal, "profoundAutismNonverbalOnly", "giftedVerbal",
                paste0(outDirFinal, "profoundNonverbalOnly_GiftedVerbal.csv"))
runLinearModels(sscGroup1 = profoundAutismBoth, sscGroup2 = verbalMildID, 
                genomicsGroup1 = splitGenomicsProfoundBoth, genomicsGroup2 = splitGenomicsMildIDVerbal, 
                subtype1 = "profoundAutismBoth", subtype2 = "mildIDVerbal",
                fileName = paste0(outDirFinal, "profoundAutismBoth_MildIDVerbal.csv"))
runLinearModels(profoundAutismBoth, verbalNoID, splitGenomicsProfoundBoth,
                splitGenomicsNoIDVerbal, "profoundAutismBoth", "noIDVerbal",
                paste0(outDirFinal, "profoundAutismBoth_NoIDVerbal.csv"))
runLinearModels(profoundAutismBoth, verbalGifted, splitGenomicsProfoundBoth,
                splitGenomicsGiftedVerbal, "profoundAutismBoth", "giftedVerbal",
                paste0(outDirFinal, "profoundAutismBoth_GiftedVerbal.csv"))
runLinearModels(sscGroup1 = profoundAutismEither, sscGroup2 = verbalMildID, 
                genomicsGroup1 = splitGenomicsProfoundEither, genomicsGroup2 = splitGenomicsMildIDVerbal, 
                subtype1 = "profoundAutismEither", subtype2 = "mildIDVerbal",
                fileName = paste0(outDirFinal, "profoundAutismEither_MildIDVerbal.csv"))
runLinearModels(profoundAutismEither, verbalNoID, splitGenomicsProfoundEither,
                splitGenomicsNoIDVerbal, "profoundAutismEither", "noIDVerbal",
                paste0(outDirFinal, "profoundAutismEither_NoIDVerbal.csv"))
runLinearModels(profoundAutismEither, verbalGifted, splitGenomicsProfoundEither,
                splitGenomicsGiftedVerbal, "profoundAutismEither", "giftedVerbal",
                paste0(outDirFinal, "profoundAutismEither_GiftedVerbal.csv"))

# Compare the profound autism groups to each other.
runLinearModels(sscGroup1 = profoundAutismModerateIDOnly, sscGroup2 = profoundAutismNonverbalOnly, 
                genomicsGroup1 = splitGenomicsProfoundModerateIDOnly, genomicsGroup2 = splitGenomicsProfoundNonverbalOnly, 
                subtype1 = "profoundModerateIDOnly", subtype2 = "profoundNonverbalOnly",
                fileName = paste0(outDirFinal, "profoundModerateIDOnly_NonverbalOnly.csv"))
runLinearModels(sscGroup1 = profoundAutismModerateIDOnly, sscGroup2 = profoundAutismBoth, 
                genomicsGroup1 = splitGenomicsProfoundModerateIDOnly, genomicsGroup2 = splitGenomicsProfoundBoth, 
                subtype1 = "profoundModerateIDOnly", subtype2 = "profoundBoth",
                fileName = paste0(outDirFinal, "profoundModerateIDOnly_ProfoundBoth.csv"))
runLinearModels(sscGroup1 = profoundAutismNonverbalOnly, sscGroup2 = profoundAutismBoth, 
                genomicsGroup1 = splitGenomicsProfoundNonverbalOnly, genomicsGroup2 = splitGenomicsProfoundBoth, 
                subtype1 = "profoundNonverbalOnly", subtype2 = "profoundBoth",
                fileName = paste0(outDirFinal, "profoundNonverbalOnly_ProfoundBoth.csv"))
runLinearModels(sscGroup1 = profoundAutismNonverbalOnly, sscGroup2 = profoundAutismEither, 
                genomicsGroup1 = splitGenomicsProfoundNonverbalOnly, genomicsGroup2 = splitGenomicsProfoundEither, 
                subtype1 = "profoundNonverbalOnly", subtype2 = "profoundEither",
                fileName = paste0(outDirFinal, "profoundNonverbalOnly_ProfoundEither.csv"))
runLinearModels(sscGroup1 = profoundAutismModerateIDOnly, sscGroup2 = profoundAutismEither, 
                genomicsGroup1 = splitGenomicsProfoundModerateIDOnly, genomicsGroup2 = splitGenomicsProfoundEither, 
                subtype1 = "profoundModerateIDOnly", subtype2 = "profoundEither",
                fileName = paste0(outDirFinal, "profoundNonverbalOnly_ProfoundEither.csv"))
runLinearModels(sscGroup1 = profoundAutismBoth, sscGroup2 = profoundAutismEither, 
                genomicsGroup1 = splitGenomicsProfoundBoth, genomicsGroup2 = splitGenomicsProfoundEither, 
                subtype1 = "profoundBoth", subtype2 = "profoundEither",
                fileName = paste0(outDirFinal, "profoundBoth_ProfoundEither.csv"))