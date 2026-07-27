library("ComplexHeatmap")

# Read eQTL data.
eQTLDir <- "../eQTL/"
eQTLRawCombCis <- readRDS(paste0(eQTLDir, "eQTLMatCis.csv_comb.RDS"))
eQTLDiffCis <- readRDS(paste0(eQTLDir, "eQTLDiffMatCis.csv_comb.RDS"))
eQTLRawCombTrans <- readRDS(paste0(eQTLDir, "eQTLMatTrans.csv_comb.RDS"))
eQTLDiffTrans <- readRDS(paste0(eQTLDir, "eQTLDiffMatTrans.csv_comb.RDS"))

# Read eQTLs persistent across random data.
randLogCPMCis <- readRDS(paste0(eQTLDir, "overlappingPairsLogCPMCis.RDS"))
randLogCPMTrans <- readRDS(paste0(eQTLDir, "overlappingPairsLogCPMTrans.RDS"))
randDiffCis <- readRDS(paste0(eQTLDir, "overlappingPairsDiffCis.RDS"))
randDiffTrans <- readRDS(paste0(eQTLDir, "overlappingPairsDiffTrans.RDS"))

# Look at overlap between the pairs.
paNames <- c("Profound, Nonverbal and ID", "Profound, Moderate ID Only")
asdNames <- c("ASD, Gifted", "ASD, No ID", "ASD, Mild ID")
findNoOverlap <- function(matDat, rand, fileName){
  mat <- attr(matDat, "data")
  randShared <- intersect(rand, rownames(mat))
  str(randShared)
  paSub <- mat[randShared, paNames]
  asdSub <- mat[randShared, asdNames]
  whichPaSpec <- intersect(which(rowSums(paSub) > 0), which(rowSums(paSub) == 0))
  str(which(rowSums(paSub) < 2))
  str(which(rowSums(asdSub) < 3))
  print(randShared[whichPaSpec])
  saveRDS(randShared[whichPaSpec], paste0(eQTLDir, fileName))
}
findNoOverlap(eQTLRawCombCis, randLogCPMCis, "paSpecLogCPMCis.RDS")
findNoOverlap(eQTLDiffCis, randDiffCis, "paSpecDiffCis.RDS")
findNoOverlap(eQTLRawCombTrans, randLogCPMTrans, "paSpecLogCPMTrans.RDS")
findNoOverlap(eQTLDiffTrans, randDiffTrans, "paSpecDiffTrans.RDS")