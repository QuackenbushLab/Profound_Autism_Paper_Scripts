# Needed to run Levene's test.
if(!require("car")){
  install.packages("car")
}
library("car")

# Read data.
sourceDirGenomics <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/SFARI/large_results/"
phenoGrp <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/SFARI/profoundAutism/phenotypeGroups/"
outDirFinal <- paste0(outDir, "/diffGeneExpressionVariance/")
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

# Formula
formulaAll <- "gene ~ subtype"

# Run linear models.
dir.create(outDirFinal)
levenesTest <- function(sscGroup1, sscGroup2, genomicsGroup1, genomicsGroup2, 
                            subtype1, subtype2, fileName){
    formula <- formulaAll
    ssc <- rbind(sscGroup1, sscGroup2)
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
                        gene = gene))
    })
    pvalues <- do.call(rbind, pvaluesList)
    pvalues$padj <- stats::p.adjust(pvalues$pval, method = "fdr")
    str(pvalues[which(pvalues$padj < 0.05), "gene"])
    write.csv(pvalues, fileName)
}

# Compare the profound autism groups to everything else.
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
 levenesTest(sscGroup1 = profoundAutismEither, sscGroup2 = do.call(rbind, list(verbalMildID, verbalNoID, verbalGifted)),
            genomicsGroup1 = splitGenomicsProfoundEither, genomicsGroup2 = do.call(cbind, list(splitGenomicsMildIDVerbal, 
                                                                                                      splitGenomicsNoIDVerbal,
                                                                                                      splitGenomicsGiftedVerbal)),
            subtype1 = "profoundAutismBoth", subtype2 = "notProfound",
            fileName = paste0(outDirFinal, "profoundAutismBoth_NotProfound.csv"))
