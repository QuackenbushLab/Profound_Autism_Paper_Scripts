eQTLDir <- "../eQTL/"
logCPM <- readRDS(paste0(eQTLDir, "/alleQTLResultsRand.RDS"))
diff <- readRDS(paste0(eQTLDir, "/alleQTLResultsRandDiff.RDS"))

doJaccard <- function(pairs){
  jaccard <- unlist(lapply(1:10, function(i){
    cis1 <- pairs[[i]][[1]]
    cis2 <- pairs[[i]][[2]]
    union <- unique(c(cis1, cis2))
    intersection <- intersect(cis1, cis2)
    jaccard <- length(intersection) / length(union)
    str(jaccard)
    return(jaccard)
  }))
  names(jaccard) <- 1:10
  str(jaccard)
  return(jaccard)
}
jaccardCis <- doJaccard(logCPM$pairsCis)
jaccardTrans <- doJaccard(logCPM$pairsTrans)
jaccardCisDiff <- doJaccard(diff$pairsCis)
jaccardTransDiff <- doJaccard(diff$pairsTrans)

overlappingPairs <- function(pairs){
  overlap <- lapply(1:10, function(i){
    cis1 <- pairs[[i]][[1]]
    cis2 <- pairs[[i]][[2]]
    return(intersect(cis1, cis2))
  })
  str(overlap)
  names(overlap) <- 1:10
  return(overlap)
}
overlappingPairsCis <- overlappingPairs(logCPM$pairsCis)
overlappingPairsTrans <- overlappingPairs(logCPM$pairsTrans)
overlappingPairsCisDiff <- overlappingPairs(diff$pairsCis)
overlappingPairsTransDiff <- overlappingPairs(diff$pairsTrans)

saveRDS(list(jaccardCis = jaccardCis,
             jaccardTrans = jaccardTrans,
             overlappingPairsCis = overlappingPairsCis,
             overlappingPairsTrans = overlappingPairsTrans), 
        paste0(eQTLDir, "/jaccardRandLogCPM.RDS"))
saveRDS(list(jaccardCis = jaccardCisDiff,
             jaccardTrans = jaccardTransDiff,
             overlappingPairsCis = overlappingPairsCisDiff,
             overlappingPairsTrans = overlappingPairsTransDiff), 
        paste0(eQTLDir, "/jaccardRandDiff.RDS"))