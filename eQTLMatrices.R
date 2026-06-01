# Read eQTL data.
eQTLDir <- "../eQTL/"
eQTLRaw <- readRDS(paste0(eQTLDir, "alleQTLResults.RDS"))
eQTLDiff <- readRDS(paste0(eQTLDir, "alleQTLDiffResults.RDS"))

# Create an UpSet plot of pairs.
library(ComplexHeatmap)
makeMatrix <- function(eQTL, matFile, type){
  
  eQTLPairList <- unique(unlist(eQTL))
  eQTLMat <- matrix(data = rep(0, length(eQTLPairList) * 7), ncol = 7)
  eQTLMat <- eQTLMat[,3:7]
  colnames(eQTLMat) <- c("Profound, Nonverbal and ID", 
                         "Profound, Moderate ID Only", "ASD, Gifted", "ASD, No ID", "ASD, Mild ID")
  rownames(eQTLMat) <- eQTLPairList
  if(type == "cis"){
    eQTLMat[eQTL[["pairsProfoundBothCis"]], "Profound, Nonverbal and ID"] <- 1
    eQTLMat[eQTL[["pairsProfoundModerateIDCis"]], "Profound, Moderate ID Only"] <- 1
    eQTLMat[eQTL[["pairsGiftedCis"]], "ASD, Gifted"] <- 1
    eQTLMat[eQTL[["pairsNoIDCis"]], "ASD, No ID"] <- 1
    eQTLMat[eQTL[["pairsMildIDCis"]], "ASD, Mild ID"] <- 1
  }else if(type == "trans"){
    eQTLMat[eQTL[["pairsProfoundBothTransSig"]], "Profound, Nonverbal and ID"] <- 1
    eQTLMat[eQTL[["pairsProfoundModerateIDTransSig"]], "Profound, Moderate ID Only"] <- 1
    eQTLMat[eQTL[["pairsGiftedTransSig"]], "ASD, Gifted"] <- 1
    eQTLMat[eQTL[["pairsNoIDTransSig"]], "ASD, No ID"] <- 1
    eQTLMat[eQTL[["pairsMildIDTransSig"]], "ASD, Mild ID"] <- 1
  }
  write.csv(eQTLMat, matFile)
  
  # Make combination matrix.
  eQTLMatComb <- make_comb_mat(eQTLMat, mode = "distinct")
  saveRDS(eQTLMatComb, paste0(matFile, "_comb.RDS"))
}
makeMatrix(eQTL = eQTLRaw, matFile = paste0(eQTLDir, "eQTLMatCis.csv"), type = "cis")
makeMatrix(eQTL = eQTLDiff, matFile = paste0(eQTLDir, "eQTLDiffMatCis.csv"), type = "cis")
makeMatrix(eQTL = eQTLRaw, matFile = paste0(eQTLDir, "eQTLMatTrans.csv"), type = "trans")
makeMatrix(eQTL = eQTLDiff, matFile = paste0(eQTLDir, "eQTLDiffMatTrans.csv"), type = "trans")