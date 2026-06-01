# Read SSC data.
outDir <- NULL

# Set colors.
bothCol <- rgb(red = 175 / 255, blue = 0 / 255, green = 93 / 255, alpha = 0.5)
nonverbalCol <- rgb(red = 187 / 255, blue = 255 / 255, green = 1 / 255, alpha = 0.5)
modIDCol <- rgb(red = 0 / 255, blue = 0 / 255, green = 160 / 255, alpha = 0.5)
mildIDCol <- rgb(red = 2 / 255, blue = 83 / 255, green = 125 / 255, alpha = 0.5)
noIDCol <- rgb(red = 0 / 255, blue = 170 / 255, green = 129 / 255, alpha = 0.5)
giftedCol <- rgb(red = 0 / 255, blue = 255 / 255, green = 0 / 255, alpha = 0.5)

# Subset to above 8.
profoundAutismModerateIDOnly <- read.csv(paste0(outDir, "/profoundAutismModerateIDOnly_above8.csv"), row.names = 1)
profoundAutismNonverbalOnly <- read.csv(paste0(outDir, "/profoundAutismNonverbalOnly_above8.csv"), row.names = 1)
profoundAutismBoth <- read.csv(paste0(outDir, "/profoundAutismBoth_above8.csv"), row.names = 1)
verbalMildID <- read.csv(paste0(outDir, "/verbalMildID_above8.csv"), row.names = 1)
verbalNoID <- read.csv(paste0(outDir, "/verbalNoID_above8.csv"), row.names = 1)
verbalGifted <- read.csv(paste0(outDir, "/verbalGifted_above8.csv"), row.names = 1)

# Read split genomics data.
genomicsDir <- NULL
splitGenomicsProfoundBoth <- read.csv(paste0(genomicsDir, "geneExpressionBoth.csv"), row.names = 1)
splitGenomicsProfoundModerateIDOnly <- read.csv(paste0(genomicsDir, "geneExpressionModerateID.csv"), row.names = 1)
splitGenomicsProfoundNonverbalOnly <- read.csv(paste0(genomicsDir, "geneExpressionNonverbal.csv"), row.names = 1)
splitGenomicsMildIDVerbal <- read.csv(paste0(genomicsDir, "geneExpressionMildIDVerbal.csv"), row.names = 1)
splitGenomicsNoIDVerbal <- read.csv(paste0(genomicsDir, "geneExpressionNoIDVerbal.csv"), row.names = 1)
splitGenomicsGiftedVerbal <- read.csv(paste0(genomicsDir, "geneExpressionGiftedVerbal.csv"), row.names = 1)
fullDataSet <- do.call(cbind, list(splitGenomicsProfoundModerateIDOnly, splitGenomicsProfoundNonverbalOnly,
                                   splitGenomicsProfoundBoth, splitGenomicsMildIDVerbal,
                                   splitGenomicsNoIDVerbal, splitGenomicsGiftedVerbal))
otherData <- do.call(cbind, list(splitGenomicsMildIDVerbal,
                                 splitGenomicsNoIDVerbal, splitGenomicsGiftedVerbal))


# Subset SSC data..
covariates <- c("sex")
subsetData <- function(dataSSC, g){
  subsetSSC <- dataSSC
  subsetSSC <- subsetSSC[unlist(lapply(colnames(g), function(samp){return(strsplit(samp, split = "X")[[1]][2])})),]
  return(subsetSSC)
}
profoundAutismModerateIDOnlySubsetSSC <- subsetData(profoundAutismModerateIDOnly, splitGenomicsProfoundModerateIDOnly)
profoundAutismNonverbalOnlySubsetSSC <- subsetData(profoundAutismNonverbalOnly, splitGenomicsProfoundNonverbalOnly)
profoundBothSubsetSSC <- subsetData(profoundAutismBoth, splitGenomicsProfoundBoth)
verbalMildIDSubsetSSC <- subsetData(verbalMildID, splitGenomicsMildIDVerbal)
verbalNoIDSubsetSSC <- subsetData(verbalNoID, splitGenomicsNoIDVerbal)
verbalGiftedSubsetSSC <- subsetData(verbalGifted, splitGenomicsGiftedVerbal)

plotFirstTwoPCs <- function(pcaSubset){
  # Get PCs.
  pc1 <- pcaSubset$x[,1]
  pc2 <- pcaSubset$x[,2]
  pca <- data.frame(pc1 = pc1, pc2 = pc2)
  rownames(pca) <- rownames(pcaSubset$x)  
  col = rep(rgb(red = 218 / 255, blue = 218 / 255, green = 218 / 255, alpha = 0.5), nrow(pca))
  col[which(make.names(rownames(pca)) %in% make.names(rownames(profoundBothSubsetSSC)))] <- bothCol
  col[which(make.names(rownames(pca)) %in% make.names(rownames(profoundAutismNonverbalOnlySubsetSSC)))] <- nonverbalCol
  col[which(make.names(rownames(pca)) %in% make.names(rownames(profoundAutismModerateIDOnlySubsetSSC)))] <- modIDCol
  col[which(make.names(rownames(pca)) %in% make.names(rownames(verbalMildIDSubsetSSC)))] <- mildIDCol
  col[which(make.names(rownames(pca)) %in% make.names(rownames(verbalNoIDSubsetSSC)))] <- noIDCol
  col[which(make.names(rownames(pca)) %in% make.names(rownames(verbalGiftedSubsetSSC)))] <- giftedCol
  
  
  # Calculate variances.
  eigs <- pcaSubset$sdev^2
  variance1 <- (eigs[1] / sum(eigs)) * 100
  variance2 <- (eigs[2] / sum(eigs)) * 100
  plot(pca[,c(1:2)], col = col, pch = 19, cex = 2,
       xlab = paste("PC 1 - % Variance:", format(round(variance1, 2), nsmall = 2)), 
       ylab = paste("PC 2 - % Variance:", format(round(variance2, 2), nsmall = 2)))
}
outDirPCA <- paste0(outDir, "../PCA_nodiff/")
pdf(paste0(outDirPCA, "fullResultPlot.pdf"))
par(mfrow = c(2, 2), mar = c(5,5,0,0))

# Do PCA.
dir.create(outDirPCA)
for(sex in unique(fullDataSSC$sex)){
  pcGenomics <- prcomp(t(fullDataSet[,paste0("X", rownames(fullDataSSC)[which(fullDataSSC$sex == sex)])]))
  str(pcGenomics)
  saveRDS(pcGenomics, paste0(outDir, "/expressionPCA_", sexCombo, ".RDS"))
  pcGenomics <- readRDS(paste0(outDir, "/expressionPCA_", sexCombo, ".RDS"))

  # For each pair of PC's, get the ratio of the Euclidean distance between the
  # profound autism samples and from the profound autism samples to the other samples.
  ratiosList <- lapply(1:ncol(pcGenomics$x), function(pc){
    
    # Get PCs.
    pcProfoundBoth <- pcGenomics$x[paste0("X", rownames(profoundBothSubsetSSC)[which(profoundBothSubsetSSC$sex == sex)]),pc]
    pcProfoundNonverbal <- pcGenomics$x[paste0("X", rownames(profoundAutismNonverbalOnlySubsetSSC)[which(profoundAutismNonverbalOnlySubsetSSC$sex == sex)]),pc]
    pcProfoundModerateID <- pcGenomics$x[paste0("X", rownames(profoundAutismModerateIDOnlySubsetSSC)[which(profoundAutismModerateIDOnlySubsetSSC$sex == sex)]),pc]
    pcOther <- pcGenomics$x[paste0("X", rownames(otherSSC)[which(otherSSC$sex == sex)]),pc]
    
    # Do a Wilcoxon test on the PC values.
    wilcoxBoth <- wilcox.test(x = pcProfoundBoth, y = pcOther)$p.value
    wilcoxNonverbal <- NA
    tryCatch({
      wilcoxNonverbal <- wilcox.test(x = pcProfoundNonverbal, y = pcOther)$p.value
    }, error = function(cond){})
    wilcoxModerateID <- wilcox.test(x = pcProfoundModerateID, y = pcOther)$p.value

    # Return.
    results <- data.frame(both = wilcoxBoth,
                          nonverbalOnly = wilcoxNonverbal,
                          moderateIDOnly = wilcoxModerateID)
    return(results)
  })
  pvals <- do.call(rbind, ratiosList)
  pvals$padjBoth <- p.adjust(pvals$both, method = "fdr")
  pvals$padjNonverbal <- p.adjust(pvals$nonverbalOnly, method = "fdr")
  pvals$padjModerateID <- p.adjust(pvals$moderateIDOnly, method = "fdr")

  # Save values.
  write.csv(pvals, paste0(outDirPCA, "/profoundPCADistribution_", sex, ".csv"))
  
  # Plot PCs.
  plotFirstTwoPCs(pcaSubset = pcGenomics)
  
  # Plot p-values.
  hist(as.numeric(pvals$padjBoth), breaks = seq(0, 1, by = 0.05), 
       xlab = "FDR-Adjusted P-Value for PC Separability",
       ylab = "Number of PCs", xlim = c(0, 1), ylim = c(0, length(pcGenomics$sdev)),
       col = bothCol, main = "")
  tryCatch({
    hist(as.numeric(pvals$padjNonverbal), breaks = seq(0, 1, by = 0.05), col = nonverbalCol, add = TRUE)
  }, error = function(cond){})
  hist(as.numeric(pvals$padjModerateID), breaks = seq(0, 1, by = 0.05), col = modIDCol, add = TRUE) 
}
dev.off()