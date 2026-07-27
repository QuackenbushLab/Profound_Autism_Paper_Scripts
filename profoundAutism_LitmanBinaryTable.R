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

#############################
# Put in a binary table.
############################
# List genes.
allGenes <- sort(unique(unlist(cellMarkersAgg)))
# Remove blanks.
allGenes <- allGenes[-1]
# Make rows.
geneMembershipList <- lapply(allGenes, function(gene){
  myDF <- t(data.frame(matrix(rep(0, 20))))
  colnames(myDF) <- toKeep
  for(i in 1:length(toKeep)){
    if(gene %in% cellMarkersAgg[[i]]){
      myDF[,i] <- 1
    }
  }
  str(myDF)
  return(myDF)
})
geneMembershipBin <- do.call(rbind, geneMembershipList)
rownames(geneMembershipBin) <- allGenes
# Save.
write.csv(geneMembershipBin, paste0(sourceDir, "membershipBin.csv"))