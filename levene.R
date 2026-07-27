# Needed to run Levene's test.
if(!require("car")){
  install.packages("car")
}
library("car")

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
formulaAll <- "gene ~ subtype"

# Run linear models.
dir.create(outDirFinal)
levenesTest <- function(sscGroup1, sscGroup2, genomicsGroup1, genomicsGroup2, 
                            subtype1, subtype2, fileName){
    formula <- formulaAll
    ssc <- rbind(sscGroup1, sscGroup2)
    rownames(ssc) <- paste0("X", rownames(ssc))
    gen <- rbind(t(genomicsGroup1), t(genomicsGroup2))
    ssc$subtype <- as.factor(c(rep(subtype1, nrow(sscGroup1)), 
                     rep(subtype2, nrow(sscGroup2))))
    shared <- intersect(rownames(ssc), rownames(gen))
    shared <- Reduce(intersect, list(shared, rownames(ssc)[which(ssc$ethnicity != "")],
                                     rownames(ssc)[which(ssc$race != "not-specified")]))
    ssc <- ssc[shared,]
    gen <- gen[shared,]
    pvaluesList <- lapply(colnames(gen), function(gene){
      fullDataSet <- ssc
      fullDataSet$gene <- gen[,gene]
      model <- car::leveneTest(gene ~ subtype, fullDataSet)
      return(data.frame(pval = model[["Pr(>F)"]][1],
                        fStat = model[["F value"]][1],
                        gene = gene))
    })
    pvalues <- do.call(rbind, pvaluesList)
    pvalues$padj <- stats::p.adjust(pvalues$pval, method = "fdr")
    str(pvalues[which(pvalues$padj < 0.05), "gene"])
    write.csv(pvalues, fileName)
}


# Compare profound autism to other groups, combined.
levenesTest(sscGroup1 = profoundAutismModerateIDOnly, sscGroup2 = do.call(rbind, list(verbalMildID, verbalNoID, verbalGifted)),
            genomicsGroup1 = splitGenomicsProfoundModerateIDOnly, genomicsGroup2 = do.call(cbind, list(splitGenomicsMildIDVerbal, 
                                                                                                       splitGenomicsNoIDVerbal,
                                                                                                       splitGenomicsGiftedVerbal)),
            subtype1 = "profoundModerateIDOnly", subtype2 = "notProfound",
            fileName = paste0(outDirFinal, "profoundModerateIDOnly_NotProfound.csv"))
levenesTest(sscGroup1 = profoundAutismNonverbalOnly, sscGroup2 = do.call(rbind, list(verbalMildID, verbalNoID, verbalGifted)),
            genomicsGroup1 = splitGenomicsProfoundNonverbalOnly, genomicsGroup2 = do.call(cbind, list(splitGenomicsMildIDVerbal, 
                                                                                                      splitGenomicsNoIDVerbal,
                                                                                                      splitGenomicsGiftedVerbal)),
            subtype1 = "profoundNonverbalOnly", subtype2 = "notProfound",
            fileName = paste0(outDirFinal, "profoundNonverbalOnly_NotProfound.csv"))
levenesTest(sscGroup1 = profoundAutismBoth, sscGroup2 = do.call(rbind, list(verbalMildID, verbalNoID, verbalGifted)),
            genomicsGroup1 = splitGenomicsProfoundBoth, genomicsGroup2 = do.call(cbind, list(splitGenomicsMildIDVerbal, 
                                                                                             splitGenomicsNoIDVerbal,
                                                                                             splitGenomicsGiftedVerbal)),
            subtype1 = "profoundAutismBoth", subtype2 = "notProfound",
            fileName = paste0(outDirFinal, "profoundAutismBoth_NotProfound.csv"))

# Analyze the "not profound" group.
nonverbalOnly_notProfound <- read.csv(paste0(outDirFinal, "profoundNonverbalOnly_NotProfound.csv"), row.names = 1)
nonverbalOnly_notProfoundSig <- nonverbalOnly_notProfound[which(nonverbalOnly_notProfound$padj < 0.05),]
rownames(nonverbalOnly_notProfoundSig) <- nonverbalOnly_notProfoundSig$gene
profoundNonverbalGenomics <- splitGenomicsProfoundNonverbalOnly
notProfoundGenomics <- do.call(cbind, list(splitGenomicsMildIDVerbal, splitGenomicsNoIDVerbal, splitGenomicsGiftedVerbal))
genesVarianceLower <- lapply(rownames(nonverbalOnly_notProfoundSig), function(row){
  varProfound <- var(unlist(profoundNonverbalGenomics[row,]))
  varNotProfound <- var(unlist(notProfoundGenomics[row,]))
  retval <- NULL
  if(varNotProfound > varProfound){
    retval <- as.data.frame(nonverbalOnly_notProfoundSig[row,])
    retval$varProfound <- varProfound
    retval$varNotProfound <- varNotProfound
  }
  return(retval)
})
genesVarianceHigher <- lapply(rownames(nonverbalOnly_notProfoundSig), function(row){
  varProfound <- var(unlist(profoundNonverbalGenomics[row,]))
  varNotProfound <- var(unlist(notProfoundGenomics[row,]))
  retval <- NULL
  if(varNotProfound < varProfound){
    retval <- as.data.frame(nonverbalOnly_notProfoundSig[row,])
    retval$varProfound <- varProfound
    retval$varNotProfound <- varNotProfound
  }
  return(retval)
})
genesVarianceLowerDF <- do.call(rbind, genesVarianceLower)
genesVarianceHigherDF <- do.call(rbind, genesVarianceHigher)
write.csv(genesVarianceLowerDF, paste0(outDirFinal, "profoundNonverbalOnly_lowerVariance.csv"))
write.csv(genesVarianceHigherDF, paste0(outDirFinal, "profoundNonverbalOnly_higherVariance.csv"))