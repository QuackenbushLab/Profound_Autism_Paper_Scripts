# Read in all of the Litman files.
sourceDir <- NULL
asdRiskGenes <- read.table(paste0(sourceDir, "/ASD_risk_genes_TADA_FDR0.3.bed"), sep = "\t",
                           header = TRUE)$name
constrainedPLIScore <- read.table(paste0(sourceDir, "/Constrained_PLIScoreOver0.9.bed"), sep = "\t",
                           header = TRUE)$name
developmentalDelay <- read.table(paste0(sourceDir, "/Developmental_delay_DDD.bed"), sep = "\t",
                                  header = TRUE)$name
FMRP_targets <- read.table(paste0(sourceDir, "/FMRP_targets_Darnell2011.bed"), sep = "\t",
                                 header = TRUE)$name
PSD <- read.table(paste0(sourceDir, "/PSD_Genes2Cognition.bed"), sep = "\t",
                           header = TRUE)$name

# Map to ENSEMBL IDs.
#Add SYMBOLS.
library(org.Hs.eg.db)
mapToIDs <- function(genes){
  ens2symbol <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                    key=genes, 
                                    column="ENSEMBL",
                                    keytype="SYMBOL")
  str(ens2symbol)
  return(unname(ens2symbol$ENSEMBL))
}
mapToIDs <- function(genes){
  ens2symbol <- unname(AnnotationDbi::mapIds(org.Hs.eg.db,
                                    key=genes, 
                                    column="ENSEMBL",
                                    keytype="SYMBOL"))
  return(ens2symbol[which(!is.na(ens2symbol))])
}
asdRiskGenesEnsembl <- mapToIDs(asdRiskGenes)
constrainedPLIScoreEnsembl <- mapToIDs(constrainedPLIScore)
developmentalDelayEnsembl <- mapToIDs(developmentalDelay)
FMRP_targetsEnsembl <- mapToIDs(FMRP_targets)
PSD_Ensembl <- mapToIDs(PSD)
toEnsembl <- list(asdRiskGenes = asdRiskGenesEnsembl, constrainedPLIScore = constrainedPLIScoreEnsembl,
                   developmentalDelay = developmentalDelayEnsembl, FMRP_targets = FMRP_targetsEnsembl,
                   PSD = PSD_Ensembl)

# Read ENSEMBL gene files.
satterstrom <- read.csv(paste0(sourceDir, "/satterstrom_2020_102_ASD_genes.csv"))
satterstromASD <- satterstrom[which(satterstrom$ASD_vs_DDID == "ASD"),"ensembl_gene_id"]
satterstromDDID <- satterstrom[which(satterstrom$ASD_vs_DDID == "DDID"),"ensembl_gene_id"]
sfari <- read.csv(paste0(sourceDir, "/SFARI_genes.csv"))[,"ensembl.id"]
ensembl <- append(toEnsembl, list(satterstromASD = satterstromASD, satterstromDDID = satterstromDDID,
                sfari = sfari))

# Arrange the cell marker data.
getCellMarkersData <- function(cellMarkers, groupName){
  cellMarkersUp <- cellMarkers[which(cellMarkers$trend_class == "up"), "Genes"]
  cellMarkersDown <- cellMarkers[which(cellMarkers$trend_class == "down"), "Genes"]
  cellMarkersTransUp <- cellMarkers[which(cellMarkers$trend_class == "trans_up"), "Genes"]
  cellMarkersTransDown <- cellMarkers[which(cellMarkers$trend_class == "trans_down"), "Genes"]
  results <- list(cellMarkersUp, cellMarkersDown, cellMarkersTransUp, cellMarkersTransDown)
  names <- paste(groupName, c("up", "down", "trans_up", "trans_down"), sep = "_")
  names(results) <- names
  return(results)
}
cellMarkersA <- append(ensembl, getCellMarkersData(read.csv(paste0(sourceDir, "/cell_markers_developmental_stages_A.csv")), "A"))
cellMarkersB <- append(cellMarkersA, getCellMarkersData(read.csv(paste0(sourceDir, "/cell_markers_developmental_stages_B.csv")), "B"))
cellMarkersC <- append(cellMarkersB, getCellMarkersData(read.csv(paste0(sourceDir, "/cell_markers_developmental_stages_C.csv")), "C"))
cellMarkersD <- append(cellMarkersC, getCellMarkersData(read.csv(paste0(sourceDir, "/cell_markers_developmental_stages_D.csv")), "D"))
cellMarkersE <- append(cellMarkersD, getCellMarkersData(read.csv(paste0(sourceDir, "/cell_markers_developmental_stages_E.csv")), "E"))
cellMarkersF <- append(cellMarkersE, getCellMarkersData(read.csv(paste0(sourceDir, "/cell_markers_developmental_stages_F.csv")), "F"))
cellMarkersG <- append(cellMarkersF, getCellMarkersData(read.csv(paste0(sourceDir, "/cell_markers_developmental_stages_G.csv")), "G"))
cellMarkersH <- append(cellMarkersG, getCellMarkersData(read.csv(paste0(sourceDir, "/cell_markers_developmental_stages_H.csv")), "H"))
cellMarkersI <- append(cellMarkersH, getCellMarkersData(read.csv(paste0(sourceDir, "/cell_markers_developmental_stages_I.csv")), "I"))
cellMarkersJ <- append(cellMarkersI, getCellMarkersData(read.csv(paste0(sourceDir, "/cell_markers_developmental_stages_J.csv")), "J"))
cellMarkersK <- append(cellMarkersJ, getCellMarkersData(read.csv(paste0(sourceDir, "/cell_markers_developmental_stages_K.csv")), "K"))
cellMarkersL <- append(cellMarkersK, getCellMarkersData(read.csv(paste0(sourceDir, "/cell_markers_developmental_stages_L.csv")), "L"))
cellMarkersM <- append(cellMarkersL, getCellMarkersData(read.csv(paste0(sourceDir, "/cell_markers_developmental_stages_M.csv")), "M"))
cellMarkersN <- append(cellMarkersM, getCellMarkersData(read.csv(paste0(sourceDir, "/cell_markers_developmental_stages_N.csv")), "N"))
saveRDS(cellMarkersN, paste0(sourceDir, "/litmanCellMarkers.RDS"))

# Aggregate.
cellMarkersAgg <- cellMarkersN
cellMarkersAgg[["Glia_up"]] <- c(cellMarkersN[["A_up"]], cellMarkersN[["H_up"]],
                                 cellMarkersN[["I_up"]], cellMarkersN[["J_up"]])
cellMarkersAgg[["Glia_down"]] <- c(cellMarkersN[["A_down"]], cellMarkersN[["H_down"]],
                                 cellMarkersN[["I_down"]], cellMarkersN[["J_down"]])
cellMarkersAgg[["Glia_trans_down"]] <- c(cellMarkersN[["A_trans_down"]], cellMarkersN[["H_trans_down"]],
                                   cellMarkersN[["I_trans_down"]], cellMarkersN[["J_trans_down"]])
cellMarkersAgg[["Glia_trans_up"]] <- c(cellMarkersN[["A_trans_up"]], cellMarkersN[["H_trans_up"]],
                                         cellMarkersN[["I_trans_up"]], cellMarkersN[["J_trans_up"]])
cellMarkersAgg[["Inhibitory_interneuron_up"]] <- c(cellMarkersN[["B_up"]], cellMarkersN[["G_up"]],
                                 cellMarkersN[["K_up"]], cellMarkersN[["L_up"]],
                                 cellMarkersN[["M_up"]], cellMarkersN[["N_up"]])
cellMarkersAgg[["Inhibitory_interneuron_down"]] <- c(cellMarkersN[["B_down"]], cellMarkersN[["G_down"]],
                                                   cellMarkersN[["K_down"]], cellMarkersN[["L_down"]],
                                                   cellMarkersN[["M_down"]], cellMarkersN[["N_down"]])
cellMarkersAgg[["Inhibitory_interneuron_trans_down"]] <- c(cellMarkersN[["B_trans_down"]], cellMarkersN[["G_trans_down"]],
                                                     cellMarkersN[["K_trans_down"]], cellMarkersN[["L_trans_down"]],
                                                     cellMarkersN[["M_trans_down"]], cellMarkersN[["N_trans_down"]])
cellMarkersAgg[["Inhibitory_interneuron_trans_up"]] <- c(cellMarkersN[["B_trans_up"]], cellMarkersN[["G_trans_up"]],
                                                           cellMarkersN[["K_trans_up"]], cellMarkersN[["L_trans_up"]],
                                                           cellMarkersN[["M_trans_up"]], cellMarkersN[["N_trans_up"]])
cellMarkersAgg[["Principal_excitatory_neuron_up"]] <- c(cellMarkersN[["C_up"]], cellMarkersN[["D_up"]],
                                                   cellMarkersN[["E_up"]], cellMarkersN[["F_up"]])
cellMarkersAgg[["Principal_excitatory_neuron_down"]] <- c(cellMarkersN[["C_down"]], cellMarkersN[["D_down"]],
                                                        cellMarkersN[["E_down"]], cellMarkersN[["F_down"]])
cellMarkersAgg[["Principal_excitatory_neuron_trans_down"]] <- c(cellMarkersN[["C_trans_down"]], cellMarkersN[["D_trans_down"]],
                                                        cellMarkersN[["E_trans_down"]], cellMarkersN[["F_trans_down"]])
cellMarkersAgg[["Principal_excitatory_neuron_trans_up"]] <- c(cellMarkersN[["C_trans_up"]], cellMarkersN[["D_trans_up"]],
                                                                cellMarkersN[["E_trans_up"]], cellMarkersN[["F_trans_up"]])
toKeep <- c("asdRiskGenes", "constrainedPLIScore", "developmentalDelay", "FMRP_targets",
            "PSD", "satterstromASD", "satterstromDDID", "sfari", "Glia_up", "Glia_down",
            "Glia_trans_up", "Glia_trans_down", "Inhibitory_interneuron_up",
            "Inhibitory_interneuron_down", "Inhibitory_interneuron_trans_up",
            "Inhibitory_interneuron_trans_down", "Principal_excitatory_neuron_up", 
            "Principal_excitatory_neuron_down", "Principal_excitatory_neuron_trans_up",
            "Principal_excitatory_neuron_trans_down")
cellMarkersAgg <- cellMarkersAgg[toKeep]

# Check for overlap.
cellMarkersAggVec <- do.call(c, cellMarkersAgg)
geneTable <- table(cellMarkersAggVec)
repeatedGenes <- names(geneTable)[which(geneTable > 1)]
cellMarkersAggUnique <- lapply(cellMarkersAgg, function(geneSet){
  return(setdiff(geneSet, repeatedGenes))
})
names(cellMarkersAggUnique) <- names(cellMarkersAgg)
cellMarkersAggUnique <- cellMarkersAggUnique[c(1:5, 8:20)]

# Read results.
outDirFinal <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/SFARI/profoundAutism/diffGeneExpression_ageBinned/"
iq1 <- read.csv(paste0(outDirFinal, "IQ1_nonverbalToBoth.csv"), row.names = 1)
iq2 <- read.csv(paste0(outDirFinal, "IQ2_mildToModerate.csv"), row.names = 1)
iq3 <- read.csv(paste0(outDirFinal, "IQ3_noToMild.csv"), row.names = 1)
iq4 <- read.csv(paste0(outDirFinal, "IQ5_giftedToNo.csv"), row.names = 1)
speech1 <- read.csv(paste0(outDirFinal, "Speech1_moderateToBoth.csv"), row.names = 1)
speech2 <- read.csv(paste0(outDirFinal, "Speech2_mildToNonverbal.csv"), row.names = 1)
speech3 <- read.csv(paste0(outDirFinal, "Speech3_noToNonverbal.csv"), row.names = 1)
speech4 <- read.csv(paste0(outDirFinal, "Speech4_giftedToNonverbal.csv"), row.names = 1)
intersectional1 <- read.csv(paste0(outDirFinal, "Intersectional1_mildToBoth.csv"), row.names = 1)
intersectional2 <- read.csv(paste0(outDirFinal, "Intersectional2_noToBoth.csv"), row.names = 1)
intersectional3 <- read.csv(paste0(outDirFinal, "Intersectional3_giftedToBoth.csv"), row.names = 1)

# Run pathway analysis.
sourceDirResults <- paste0(sourceDir, "/litmanGSEA_lm/")
dir.create(sourceDirResults)
library(fgsea)
set.seed(1)
RunPathwayAnalysis <- function(diffScores, pathwayFile, scoreType, pathways){
  
  # Run pathway analysis.
  ranks <- sort(diffScores)
  pathwayResult <- fgsea::fgsea(pathways = pathways, stats = diffScores,
                                scoreType = scoreType, sampleSize = 3)

  # Compile the "leading edge" vector into a list.
  leadingEdge <- unlist(lapply(1:length(pathwayResult$leadingEdge), function(i){
    return(paste(pathwayResult$leadingEdge[i][[1]], collapse = "; "))
  }))
  pathwayResultDf <- data.frame(pathway = pathwayResult$pathway,
                                pval = pathwayResult$pval,
                                padj = pathwayResult$padj,
                                pvalErrorSD = pathwayResult$log2err,
                                enrichmentScore = pathwayResult$ES,
                                normalizedEnrichmentScore = pathwayResult$NES,
                                remainingGeneCount = pathwayResult$size,
                                leadingGenesDrivingEnrichment = leadingEdge)
  if(scoreType == "neg"){
    write.csv(pathwayResultDf, paste0(pathwayFile, "_neg.csv"))
  }else{
    write.csv(pathwayResultDf, paste0(pathwayFile, "_pos.csv"))
  }
  
  # Return result.
  return(pathwayResultDf)
}
DE2PathwayAnalysisSex <- function(DEresults, pathwayFile, scoreType, pathways){
  
  # Extract score.
  score <- DEresults[,2] * -1 * log10(DEresults$pval) * DEresults$rsq
  names(score) <- unlist(lapply(DEresults$gene, function(gene){
    return(strsplit(gene, split = ".", fixed = TRUE)[[1]][1])
  }))
  
  # Add jitter to prevent ties.
  scoreTbl <- table(score)
  mode <- as.numeric(names(scoreTbl)[which.max(scoreTbl)])
  whichNotMode <- which(abs(score - mode) > 0.000000001)
  whichMode <- setdiff(1:length(score), whichNotMode)
  minDistFromMode <- abs(min(score[whichNotMode] - mode))
  jitter <- rnorm(length(whichMode), mean = mode, sd = minDistFromMode / 10)
  score[whichMode] <- jitter

  # Dedup.
  geneTable <- table(names(score))
  dupGenes <- unique(names(geneTable)[which(geneTable > 1)])
  toRemove <- unlist(lapply(dupGenes, function(dupGene){
    scoreIdx <- which(names(score) == dupGene)
    maxScoreIdx <- which.max(score[scoreIdx])
    rm <- setdiff(scoreIdx, scoreIdx[maxScoreIdx])
    return(rm)
  }))
  score <- score[setdiff(1:length(score), toRemove)]

  # Remove NA values.
  score <- score[which(!is.na(score))]
  pathways <- RunPathwayAnalysis(score, pathwayFile, scoreType, pathways)
  return(pathways)
}

# Run the pathway analysis with randomized jitter.
for(i in 1:10){
  # Compute the scores (Positive, full pathways).
  pwIQ1Pos <- DE2PathwayAnalysisSex(DEresults = iq1,pathwayFile = paste0(sourceDirResults, "iq1_fullpw_", i ), 
                                    scoreType = "pos", pathways = cellMarkersAgg)
  pwIQ2Pos <- DE2PathwayAnalysisSex(DEresults = iq2,pathwayFile = paste0(sourceDirResults, "iq2_fullpw_", i ), 
                                    scoreType = "pos", pathways = cellMarkersAgg)
  pwIQ3Pos <- DE2PathwayAnalysisSex(DEresults = iq3,pathwayFile = paste0(sourceDirResults, "iq3_fullpw_", i ), 
                                    scoreType = "pos", pathways = cellMarkersAgg)
  pwIQ4Pos <- DE2PathwayAnalysisSex(DEresults = iq4,pathwayFile = paste0(sourceDirResults, "iq4_fullpw_", i ), 
                                    scoreType = "pos", pathways = cellMarkersAgg)
  pwSpeech1Pos <- DE2PathwayAnalysisSex(DEresults = speech1,pathwayFile = paste0(sourceDirResults, "speech1_fullpw_", i ), 
                                        scoreType = "pos", pathways = cellMarkersAgg)
  pwSpeech2Pos <- DE2PathwayAnalysisSex(DEresults = speech2,pathwayFile = paste0(sourceDirResults, "speech2_fullpw_", i ), 
                                        scoreType = "pos", pathways = cellMarkersAgg)
  pwSpeech3Pos <- DE2PathwayAnalysisSex(DEresults = speech3,pathwayFile = paste0(sourceDirResults, "speech3_fullpw_", i ), 
                                        scoreType = "pos", pathways = cellMarkersAgg)
  pwSpeech4Pos <- DE2PathwayAnalysisSex(DEresults = speech4,pathwayFile = paste0(sourceDirResults, "speech4_fullpw_", i ), 
                                        scoreType = "pos", pathways = cellMarkersAgg)
  pwInt1Pos <- DE2PathwayAnalysisSex(DEresults = intersectional1,
                                     pathwayFile = paste0(sourceDirResults, "intersectional1_fullpw_", i ), 
                                     scoreType = "pos", pathways = cellMarkersAgg)
  pwInt2Pos <- DE2PathwayAnalysisSex(DEresults = intersectional2,
                                     pathwayFile = paste0(sourceDirResults, "intersectional2_fullpw_", i ), 
                                     scoreType = "pos", pathways = cellMarkersAgg)
  pwInt3Pos <- DE2PathwayAnalysisSex(DEresults = intersectional3,
                                     pathwayFile = paste0(sourceDirResults, "intersectional3_fullpw_", i ), 
                                     scoreType = "pos", pathways = cellMarkersAgg)
  
  # Positive, non-overlapping pathways
  pwIQ1PosU <- DE2PathwayAnalysisSex(DEresults = iq1,pathwayFile = paste0(sourceDirResults, "iq1_upw_", i ), 
                                     scoreType = "pos", pathways = cellMarkersAggUnique)
  pwIQ2PosU <- DE2PathwayAnalysisSex(DEresults = iq2,pathwayFile = paste0(sourceDirResults, "iq2_upw_", i ), 
                                     scoreType = "pos", pathways = cellMarkersAggUnique)
  pwIQ3PosU <- DE2PathwayAnalysisSex(DEresults = iq3,pathwayFile = paste0(sourceDirResults, "iq3_upw_", i ), 
                                     scoreType = "pos", pathways = cellMarkersAggUnique)
  pwIQ4PosU <- DE2PathwayAnalysisSex(DEresults = iq4,pathwayFile = paste0(sourceDirResults, "iq4_upw_", i ), 
                                     scoreType = "pos", pathways = cellMarkersAggUnique)
  pwSpeech1PosU <- DE2PathwayAnalysisSex(DEresults = speech1,pathwayFile = paste0(sourceDirResults, "speech1_upw_", i ), 
                                         scoreType = "pos", pathways = cellMarkersAggUnique)
  pwSpeech2PosU <- DE2PathwayAnalysisSex(DEresults = speech2,pathwayFile = paste0(sourceDirResults, "speech2_upw_", i ), 
                                         scoreType = "pos", pathways = cellMarkersAggUnique)
  pwSpeech3PosU <- DE2PathwayAnalysisSex(DEresults = speech3,pathwayFile = paste0(sourceDirResults, "speech3_upw_", i ), 
                                         scoreType = "pos", pathways = cellMarkersAggUnique)
  pwSpeech4PosU <- DE2PathwayAnalysisSex(DEresults = speech4,pathwayFile = paste0(sourceDirResults, "speech4_upw_", i ), 
                                         scoreType = "pos", pathways = cellMarkersAggUnique)
  pwInt1PosU <- DE2PathwayAnalysisSex(DEresults = intersectional1,
                                      pathwayFile = paste0(sourceDirResults, "intersectional1_upw_", i ), 
                                      scoreType = "pos", pathways = cellMarkersAggUnique)
  pwInt2PosU <- DE2PathwayAnalysisSex(DEresults = intersectional2,
                                      pathwayFile = paste0(sourceDirResults, "intersectional2_upw_", i ), 
                                      scoreType = "pos", pathways = cellMarkersAggUnique)
  pwInt3PosU <- DE2PathwayAnalysisSex(DEresults = intersectional3,
                                      pathwayFile = paste0(sourceDirResults, "intersectional3_upw_", i ), 
                                      scoreType = "pos", pathways = cellMarkersAggUnique)
  
  # Compute the scores (negative, full pathways)
  pwIQ1Neg <- DE2PathwayAnalysisSex(DEresults = iq1,pathwayFile = paste0(sourceDirResults, "iq1_fullpw_", i ), 
                                    scoreType = "neg", pathways = cellMarkersAgg)
  pwIQ2Neg <- DE2PathwayAnalysisSex(DEresults = iq2,pathwayFile = paste0(sourceDirResults, "iq2_fullpw_", i ), 
                                    scoreType = "neg", pathways = cellMarkersAgg)
  pwIQ3Neg <- DE2PathwayAnalysisSex(DEresults = iq3,pathwayFile = paste0(sourceDirResults, "iq3_fullpw_", i ), 
                                    scoreType = "neg", pathways = cellMarkersAgg)
  pwIQ4Neg <- DE2PathwayAnalysisSex(DEresults = iq4,pathwayFile = paste0(sourceDirResults, "iq4_fullpw_", i ), 
                                    scoreType = "neg", pathways = cellMarkersAgg)
  pwSpeech1Neg <- DE2PathwayAnalysisSex(DEresults = speech1,
                                        pathwayFile = paste0(sourceDirResults, "speech1_fullpw_", i ), 
                                        scoreType = "neg", pathways = cellMarkersAgg)
  pwSpeech2Neg <- DE2PathwayAnalysisSex(DEresults = speech2,
                                        pathwayFile = paste0(sourceDirResults, "speech2_fullpw_", i ), 
                                        scoreType = "neg", pathways = cellMarkersAgg)
  pwSpeech3Neg <- DE2PathwayAnalysisSex(DEresults = speech3,
                                        pathwayFile = paste0(sourceDirResults, "speech3_fullpw_", i ), 
                                        scoreType = "neg", pathways = cellMarkersAgg)
  pwSpeech4Neg <- DE2PathwayAnalysisSex(DEresults = speech4,
                                        pathwayFile = paste0(sourceDirResults, "speech4_fullpw_", i ), 
                                        scoreType = "neg", pathways = cellMarkersAgg)
  pwInt1Neg <- DE2PathwayAnalysisSex(DEresults = intersectional1,
                                     pathwayFile = paste0(sourceDirResults, "intersectional1_fullpw_", i ), 
                                     scoreType = "neg", pathways = cellMarkersAgg)
  pwInt2Neg <- DE2PathwayAnalysisSex(DEresults = intersectional2,
                                     pathwayFile = paste0(sourceDirResults, "intersectional2_fullpw_", i ), 
                                     scoreType = "neg", pathways = cellMarkersAgg)
  pwInt3Neg <- DE2PathwayAnalysisSex(DEresults = intersectional3,
                                     pathwayFile = paste0(sourceDirResults, "intersectional3_fullpw_", i ), 
                                     scoreType = "neg", pathways = cellMarkersAgg)
  
  # Negative, non-overlapping pathways
  pwIQ1NegU <- DE2PathwayAnalysisSex(DEresults = iq1,pathwayFile = paste0(sourceDirResults, "iq1_upw_", i ), 
                                     scoreType = "neg", pathways = cellMarkersAggUnique)
  pwIQ2NegU <- DE2PathwayAnalysisSex(DEresults = iq2,pathwayFile = paste0(sourceDirResults, "iq2_upw_", i ), 
                                     scoreType = "neg", pathways = cellMarkersAggUnique)
  pwIQ3NegU <- DE2PathwayAnalysisSex(DEresults = iq3,pathwayFile = paste0(sourceDirResults, "iq3_upw_", i ), 
                                     scoreType = "neg", pathways = cellMarkersAggUnique)
  pwIQ4NegU <- DE2PathwayAnalysisSex(DEresults = iq4,pathwayFile = paste0(sourceDirResults, "iq4_upw_", i ), 
                                     scoreType = "neg", pathways = cellMarkersAggUnique)
  pwSpeech1NegU <- DE2PathwayAnalysisSex(DEresults = speech1,
                                         pathwayFile = paste0(sourceDirResults, "speech1_upw_", i ), 
                                         scoreType = "neg", pathways = cellMarkersAggUnique)
  pwSpeech2NegU <- DE2PathwayAnalysisSex(DEresults = speech2,
                                         pathwayFile = paste0(sourceDirResults, "speech2_upw_", i ), 
                                         scoreType = "neg", pathways = cellMarkersAggUnique)
  pwSpeech3NegU <- DE2PathwayAnalysisSex(DEresults = speech3,
                                         pathwayFile = paste0(sourceDirResults, "speech3_upw_", i ), 
                                         scoreType = "neg", pathways = cellMarkersAggUnique)
  pwSpeech4NegU <- DE2PathwayAnalysisSex(DEresults = speech4,
                                         pathwayFile = paste0(sourceDirResults, "speech4_upw_", i ), 
                                         scoreType = "neg", pathways = cellMarkersAggUnique)
  pwInt1NegU <- DE2PathwayAnalysisSex(DEresults = intersectional1,
                                      pathwayFile = paste0(sourceDirResults, "intersectional1_upw_", i ), 
                                      scoreType = "neg", pathways = cellMarkersAggUnique)
  pwInt2NegU <- DE2PathwayAnalysisSex(DEresults = intersectional2,
                                      pathwayFile = paste0(sourceDirResults, "intersectional2_upw_", i ), 
                                      scoreType = "neg", pathways = cellMarkersAggUnique)
  pwInt3NegU <- DE2PathwayAnalysisSex(DEresults = intersectional3,
                                      pathwayFile = paste0(sourceDirResults, "intersectional3_upw_", i ), 
                                      scoreType = "neg", pathways = cellMarkersAggUnique)
}

# Paste together adjusted p-values.
logTransform <- function(dat){
  return(-1 * log10(dat$padj))
}
makeHeatmap <- function(dfPos, dfNeg){
  # Make negatives actually negative.
  dfNeg <- -1 * dfNeg
  
  # Combine data frames.
  matPos <- as.matrix(dfPos)
  matNeg <- as.matrix(dfNeg)
  whichNegGreater <- which(pmax(matPos, abs(matNeg)) > matPos)
  matPos[whichNegGreater] <- matNeg[whichNegGreater]
  dfAll <- as.data.frame(matPos)
  
  # Plot.
  library(gplots)
  # Set up color.
  thr <- 0#-log10(0.05)
  nbins <- 201 
  vmin <- min(dfAll, na.rm = TRUE)
  vmax <- max(dfAll, na.rm = TRUE)
  # split palette so white lands exactly at 'thr'
  n_below <- max(2, round(nbins * (thr - vmin) / (vmax - vmin)))
  n_above <- nbins - n_below
  cols <- c(
    colorRampPalette(c("#2166AC", "#FFFFFF"))(n_below),  # blue → white
    colorRampPalette(c("#FFFFFF", "#B2182B"))(n_above)   # white → red
  )
  brks <- c(
    seq(vmin, thr, length.out = n_below + 1),
    seq(thr, vmax, length.out = n_above + 1)[-1]
  )
  # Make heatmap.
  heatmap.2(t(as.matrix(dfAll)),
            Rowv = FALSE, Colv = FALSE, dendrogram = "none",
            trace = "none", scale = "none", breaks = brks,
            col = cols, key = TRUE, density.info = "none",
            
            
            # keep row/col names from being clipped
            margins = c(12, 14),   # c(bottom for columns, left for rows) in lines
            cexCol = 1.3, cexRow = 1.3,
            keysize = 0.5,
            # optional: make the key compact so labels aren't squeezed
            key.title = "Signed LogQ", key.xlab = "", key.par = list(mar = c(3,3,1,1))
  )
}
makeDataFrame <- function(pwIQ1, pwIQ2, pwIQ3, pwIQ4, pwSpeech1, pwSpeech2, pwSpeech3,
                        pwSpeech4, pwInt1, pwInt2, pwInt3, removeTwo){
  dfAll <- data.frame(IQ1 = logTransform(pwIQ1), 
                      IQ2 = logTransform(pwIQ2), IQ3 = logTransform(pwIQ3),
                      IQ4 = logTransform(pwIQ4), speech1 = logTransform(pwSpeech1),
                      speech2 = logTransform(pwSpeech2), speech3 = logTransform(pwSpeech3),
                      speech4 = logTransform(pwSpeech4),
                      intersectional1 = logTransform(pwInt1), intersectional2 = logTransform(pwInt2),
                      intersectional3 = logTransform(pwInt3), row.names = pwIQ1$pathway)
  
  colnames(dfAll) <- c("IQ1 - Nonspeaking",
                       "IQ2 - Mild to Mod-P ID",
                       "IQ3 - Avg IQ to Mild ID",
                       "IQ4 - Gifted to Avg IQ", 
                       "Speech1 - Mod-P ID",
                       "Speech2 - Mild ID",
                       "Speech3 - Avg IQ",
                       "Speech4 - Gifted",
                       "Comb1 - Mild ID",
                       "Comb2 - Avg IQ",
                       "Comb3 - Gifted")
  str(dfAll)
  if(removeTwo){
    rownames(dfAll) <- c("FMRP Targets", "Glial (Down)",
                           "Glial (Osc Down)", "Glial (Osc Up)",
                           "Glial (Up)", "II (Down)",
                           "II (Osc Down)", 
                           "II (Osc Up)", "II (Up)",
                           "PD",
                           "PEN (Down)", "PEN (Osc Down)",
                           "PEN (Osc Up)", "PEN (Up)",
                           "Autism (Sanders)", "HEC", "DD (Werling)",
                           "Autism (Satterstrom)", "DD/ID (Satterstrom)",
                           "Autism (SFARI)")
    dfAll <- dfAll[c(6:9, 11:14, 2:5, 10, 17, 1, 15, 18, 16),]
  }else{
    rownames(dfAll) <- c("FMRP Targets", "Glial (Down)",
                         "Glial (Osc Down)", "Glial (Osc Up)",
                         "Glial (Up)", "II (Down)",
                         "II (Osc Down)", 
                         "II (Osc Up)", "II (Up)",
                         "PD",
                         "PEN (Down)", "PEN (Osc Down)",
                         "PEN (Osc Up)", "PEN (Up)",
                         "Autism (Sanders)", "HEC", "DD (Werling)",
                         "Autism (Satterstrom)", "DD/ID (Satterstrom)",
                         "Autism (SFARI)")
    dfAll <- dfAll[c(6:9, 11:14, 2:5, 10, 17, 19, 1, 15, 18, 20, 16),]
  }
  
  return(dfAll)
}
saveMeanAndSD <- function(fullOrPartialStr, posOrNeg, removeTwo){
  
  # Check mean and SD for data frame.
  allResultsList <- lapply(1:10, function(i){
    return(makeDataFrame(read.csv(paste0(sourceDirResults, "iq1_", fullOrPartialStr, "_", i , "_", posOrNeg, ".csv"), row.names = 1),
                         read.csv(paste0(sourceDirResults, "iq2_", fullOrPartialStr, "_", i , "_", posOrNeg, ".csv"), row.names = 1),
                         read.csv(paste0(sourceDirResults, "iq3_", fullOrPartialStr, "_", i , "_", posOrNeg, ".csv"), row.names = 1),
                         read.csv(paste0(sourceDirResults, "iq4_", fullOrPartialStr, "_", i , "_", posOrNeg, ".csv"), row.names = 1),
                         read.csv(paste0(sourceDirResults, "speech1_", fullOrPartialStr, "_", i , "_", posOrNeg, ".csv"), row.names = 1),
                         read.csv(paste0(sourceDirResults, "speech2_", fullOrPartialStr, "_", i , "_", posOrNeg, ".csv"), row.names = 1), 
                         read.csv(paste0(sourceDirResults, "speech3_", fullOrPartialStr, "_", i , "_", posOrNeg, ".csv"), row.names = 1),
                         read.csv(paste0(sourceDirResults, "speech4_", fullOrPartialStr, "_", i , "_", posOrNeg, ".csv"), row.names = 1),
                         read.csv(paste0(sourceDirResults, "intersectional1_", fullOrPartialStr, "_", i , "_", posOrNeg, ".csv"), row.names = 1),
                         read.csv(paste0(sourceDirResults, "intersectional2_", fullOrPartialStr, "_", i , "_", posOrNeg, ".csv"), row.names = 1),
                         read.csv(paste0(sourceDirResults, "intersectional3_", fullOrPartialStr, "_", i , "_", posOrNeg, ".csv"), row.names = 1),
                         removeTwo  = removeTwo))
  })

  # Calculate mean.
  runningSum <- as.matrix(allResultsList[[1]])
  for(i in 2:length(allResultsList)){
    runningSum <- runningSum + as.matrix(allResultsList[[i]])
  }
  meanVal <- runningSum / length(allResultsList)
  write.csv(meanVal, paste0(sourceDirResults, "mean_geneSet_pvals", "_", 
                            fullOrPartialStr, "_", posOrNeg, ".csv"))
  
  # Calculate standard deviation.
  firstMatrix <- as.matrix(allResultsList[[1]]) - meanVal
  runningSumOfSquares <- firstMatrix * firstMatrix
  for(i in 2:length(allResultsList)){
    firstMatrix <- as.matrix(allResultsList[[i]]) - meanVal
    runningSumOfSquares <- runningSumOfSquares + (firstMatrix * firstMatrix)
  }
  sdVal <- sqrt(runningSumOfSquares / length(allResultsList))
  write.csv(sdVal, paste0(sourceDirResults, "sd_geneSet_pvals", "_", 
                          fullOrPartialStr, "_", posOrNeg, ".csv"))
}
saveMeanAndSD("fullpw", "pos", removeTwo = FALSE)
saveMeanAndSD("upw", "pos", removeTwo = TRUE)
saveMeanAndSD("fullpw", "neg", removeTwo = FALSE)
saveMeanAndSD("upw", "neg", removeTwo = TRUE)

makeHeatmap(read.csv(paste0(sourceDirResults, "mean_geneSet_pvals_fullpw_pos.csv"), row.names = 1,
                     check.names = FALSE),
            read.csv(paste0(sourceDirResults, "mean_geneSet_pvals_fullpw_neg.csv"), row.names = 1,
                     check.names = FALSE))
makeHeatmap(read.csv(paste0(sourceDirResults, "mean_geneSet_pvals_upw_pos.csv"), row.names = 1,
                     check.names = FALSE),
            read.csv(paste0(sourceDirResults, "mean_geneSet_pvals_upw_neg.csv"), row.names = 1,
                     check.names = FALSE))
