library("fgsea")
library("org.Hs.eg.db")

# Paths
msigdbPath <- NULL
diffExpressionPath <- NULL
# Read in the MSigDB pathways.
msigdb <- gmtPathways(msigdbPath)

# Read in the set of genes with higher variance in the nonverbal group, for erential gene expression only.
profoundNonverbalData <- as.matrix(read.csv(paste0(diffExpressionPath, "diffGeneExpressionNonverbal.csv"),
                                  row.names = 1))
mildIdData <- read.csv(paste0(diffExpressionPath, "diffGeneExpressionMildIDVerbal.csv"),
                                  row.names = 1)
noIdData <- read.csv(paste0(diffExpressionPath, "diffGeneExpressionNoIDVerbal.csv"),
                       row.names = 1)
giftedData <- read.csv(paste0(diffExpressionPath, "diffGeneExpressionGiftedVerbal.csv"),
                     row.names = 1)
notProfoundData <- as.matrix(do.call(cbind, list(mildIdData, noIdData, giftedData)))
profoundNonverbal <- read.csv(paste0(diffExpressionPath, "diffGeneExpressionVariance/profoundNonverbalOnly_NotProfound.csv"),
                                        row.names = 1)
profoundNonverbal$logpval <- -1 * log10(profoundNonverbal$pval)
profoundNonverbal$profoundvar <- unlist(lapply(row.names(profoundNonverbalData), function(gene){
  return(var(profoundNonverbalData[gene,]))
}))
profoundNonverbal$notprofoundvar <- unlist(lapply(row.names(notProfoundData), function(gene){
  return(var(notProfoundData[gene,]))
}))
profoundNonverbal$score <- profoundNonverbal$logpval
profoundNonverbal$score[which(profoundNonverbal$profoundvar < profoundNonverbal$notprofoundvar)] <- -1 *
  profoundNonverbal$score[which(profoundNonverbal$profoundvar < profoundNonverbal$notprofoundvar)]

# Map to gene symbols.
SigSymbols <- mapIds(org.Hs.eg.db, keys = unlist(lapply(profoundNonverbal$gene, function(gene){return(strsplit(gene, split = ".", fixed = TRUE)[[1]][1])})), 
                                                column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
profoundNonverbalWithSymbols <- profoundNonverbal[which(!is.na(SigSymbols)),]
pathwayInput <- profoundNonverbalWithSymbols$score
names(pathwayInput) <- SigSymbols[which(!is.na(SigSymbols))]
toRemove <- unlist(lapply(unique(names(pathwayInput)), function(gene){
  whichMatch <- which(names(pathwayInput) == gene)
  retval <- c()
  if(length(whichMatch) > 1){
    whichMax <- which.max(pathwayInput[whichMatch])
    whichNotMax <- setdiff(1:length(whichMatch), whichMax)
    whichNotMaxMapped <- whichMatch[whichNotMax]
    retval <- whichNotMaxMapped
    #print(retval)
  }
  return(retval)
}))

# Run pathway analysis.
higherVarianceNonverbalPathways <- fgsea::fgsea(pathways = msigdb, stats = pathwayInput[setdiff(1:length(pathwayInput), toRemove)],
                                                   minSize = 5, maxSize = 100)
write.csv(higherVarianceNonverbalPathways[,c(1:5)],
          paste0(diffExpressionPath, "diffGeneExpressionVariance/profoundNonverbalOnly_higherVarianceGSEA_pvalWithDirection.csv"))

# Plot pathways.
logpvalOrdered <- -1 * log10(higherVarianceNonverbalPathways$pval)
names(logpvalOrdered) <- higherVarianceNonverbalPathways$pathway
logpvalOrdered <- logpvalOrdered[order(-logpvalOrdered)]
par(mar=c(5.1, 45, 4.1, 2.1))
p1 <- barplot(logpvalOrdered[1:10], horiz = TRUE, col = rgb(red = 174 / 255, blue = 79 / 255, green = 116 / 255), 
              las = 2, 
              xlab = "-log10(p-value)")



# Do the same, but for the non-proband-specific data.
profoundNonverbalDataNodiff <- as.matrix(read.csv(paste0(diffExpressionPath, "profoundNonverbalGenomics.csv"),
                                            row.names = 1))
mildIdDataNodiff <- read.csv(paste0(diffExpressionPath, "mildIDVerbalGenomics.csv"),
                       row.names = 1)
noIdDataNodiff <- read.csv(paste0(diffExpressionPath, "noIDVerbalGenomics.csv"),
                     row.names = 1)
giftedDataNodiff <- read.csv(paste0(diffExpressionPath, "giftedVerbalGenomics.csv"),
                       row.names = 1)
notProfoundDataNodiff <- as.matrix(do.call(cbind, list(mildIdDataNodiff, noIdDataNodiff, giftedDataNodiff)))
profoundNonverbalNodiff <- read.csv(paste0(diffExpressionPath, "geneExpressionVariance/profoundNonverbalOnly_NotProfound.csv"),
                                  row.names = 1)
profoundNonverbalNodiff$logpval <- -1 * log10(profoundNonverbalNodiff$pval)
profoundNonverbalNodiff$profoundvar <- unlist(lapply(row.names(profoundNonverbalData), function(gene){
  return(var(profoundNonverbalData[gene,]))
}))
profoundNonverbalNodiff$notprofoundvar <- unlist(lapply(row.names(notProfoundDataNodiff), function(gene){
  return(var(notProfoundDataNodiff[gene,]))
}))
profoundNonverbalNodiff$score <- profoundNonverbalNodiff$logpval
profoundNonverbalNodiff$score[which(profoundNonverbalNodiff$profoundvar < profoundNonverbalNodiff$notprofoundvar)] <- -1 *
  profoundNonverbalNodiff$score[which(profoundNonverbalNodiff$profoundvar < profoundNonverbalNodiff$notprofoundvar)]

# Map to gene symbols.
SigSymbolsNodiff <- mapIds(org.Hs.eg.db, keys = unlist(lapply(profoundNonverbalNodiff$gene, function(gene){return(strsplit(gene, split = ".", fixed = TRUE)[[1]][1])})), 
                         column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
profoundNonverbalWithSymbolsNodiff <- profoundNonverbalNodiff[which(!is.na(SigSymbolsNodiff)),]
pathwayInputNodiff <- profoundNonverbalWithSymbolsNodiff$score
names(pathwayInputNodiff) <- SigSymbolsNodiff[which(!is.na(SigSymbolsNodiff))]
toRemoveNodiff <- unlist(lapply(unique(names(pathwayInputNodiff)), function(gene){
  whichMatch <- which(names(pathwayInputNodiff) == gene)
  retval <- c()
  if(length(whichMatch) > 1){
    whichMax <- which.max(pathwayInputNodiff[whichMatch])
    whichNotMax <- setdiff(1:length(whichMatch), whichMax)
    whichNotMaxMapped <- whichMatch[whichNotMax]
    retval <- whichNotMaxMapped
    #print(retval)
  }
  return(retval)
}))

# Run pathway analysis.
higherVarianceNonverbalPathwaysNodiff <- fgsea::fgsea(pathways = msigdb, stats = pathwayInputNodiff[setdiff(1:length(pathwayInputNodiff), toRemoveNodiff)],
                                                    minSize = 5, maxSize = 100)
#write.csv(higherVarianceNonverbalPathwaysNodiff[,c(1:5)],
#          paste0(dir, "geneExpressionVariance/profoundNonverbalOnly_higherVarianceGSEA.csv"))
write.csv(higherVarianceNonverbalPathwaysNodiff[,c(1:5)],
          paste0(diffExpressionPath, "geneExpressionVariance/profoundNonverbalOnly_higherVarianceScaledGSEA.csv")))

# Plot pathways.
logpvalOrderedNodiff <- -1 * log10(higherVarianceNonverbalPathwaysNodiff$pval)
names(logpvalOrderedNodiff) <- higherVarianceNonverbalPathwaysNodiff$pathway
logpvalOrderedNodiff <- logpvalOrderedNodiff[order(-logpvalOrderedNodiff)]
par(mar=c(5.1, 30, 4.1, 2.1))
p1 <- barplot(logpvalOrderedNodiff[1:4], horiz = TRUE, col = rgb(red = 32 / 255, blue = 154 / 255, green = 95 / 255), las = 2, 
              xlab = "-log10(p-value)")


# Do the same, but for the gene expression data.
profoundNonverbalNodiff <- read.csv(paste0(diffExpressionPath, "geneExpression_ageBinnedModels/profoundNonverbalOnly_NotProfound.csv"),
                              row.names = 1)
profoundNonverbalNodiff$logpval <- -1 * log10(profoundNonverbalNodiff$pval)
profoundNonverbalNodiff$profoundvar <- unlist(lapply(row.names(profoundNonverbalNodiffData), function(gene){
  return(var(profoundNonverbalNodiffData[gene,]))
}))
profoundNonverbalNodiff$notprofoundvar <- unlist(lapply(row.names(notProfoundDataNodiff), function(gene){
  return(var(notProfoundDataNodiff[gene,]))
}))
profoundNonverbalNodiff$score <- profoundNonverbalNodiff$logpval
profoundNonverbalNodiff$score[which(profoundNonverbalNodiff$profoundvar < profoundNonverbalNodiff$notprofoundvar)] <- -1 *
  profoundNonverbalNodiff$score[which(profoundNonverbalNodiff$profoundvar < profoundNonverbalNodiff$notprofoundvar)]

# Map to gene symbols.
SigSymbolsNodiff <- mapIds(org.Hs.eg.db, keys = unlist(lapply(profoundNonverbalNodiff$gene, function(gene){return(strsplit(gene, split = ".", fixed = TRUE)[[1]][1])})), 
                     column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
profoundNonverbalNodiffWithSymbols <- profoundNonverbalNodiff[which(!is.na(SigSymbolsNodiff)),]
pathwayInputNodiff <- profoundNonverbalNodiffWithSymbols$score
names(pathwayInputNodiff) <- SigSymbolsNodiff[which(!is.na(SigSymbolsNodiff))]
toRemoveNodiff <- unlist(lapply(unique(names(pathwayInputNodiff)), function(gene){
  whichMatch <- which(names(pathwayInputNodiff) == gene)
  retval <- c()
  if(length(whichMatch) > 1){
    whichMax <- which.max(pathwayInputNodiff[whichMatch])
    whichNotMax <- setdiff(1:length(whichMatch), whichMax)
    whichNotMaxMapped <- whichMatch[whichNotMax]
    retval <- whichNotMaxMapped
    #print(retval)
  }
  return(retval)
}))

# Run pathway analysis.
higherVarianceNonverbalPathwaysNodiff <- fgsea::fgsea(pathways = msigdb, stats = pathwayInputNodiff[setdiff(1:length(pathwayInputNodiff), toRemoveNodiff)],
                                                minSize = 5, maxSize = 100)
write.csv(higherVarianceNonverbalPathwaysNodiff[,c(1:5)],
          paste0(diffExpressionPath, "geneExpression_ageBinnedModels/geneExpression.csv"))

# Plot pathways.
logpvalOrderedNodiff <- -1 * log10(higherVarianceNonverbalPathwaysNodiff$pval)
names(logpvalOrderedNodiff) <- higherVarianceNonverbalPathwaysNodiff$pathway
logpvalOrderedNodiff <- logpvalOrderedNodiff[order(-logpvalOrderedNodiff)]
par(mar=c(5.1, 30, 4.1, 2.1))
p1 <- barplot(logpvalOrderedNodiff[1:4], horiz = TRUE, col = rgb(red = 32 / 255, blue = 154 / 255, green = 95 / 255), las = 2, 
              xlab = "-log10(p-value)")
