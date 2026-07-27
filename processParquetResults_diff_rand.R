if(!require("arrow")){
  install.packages("arrow")
}
library("arrow")

# Read files.
eQTLDir <- NULL
cutoff <- 0.05
readAllCis <- function(fpath){
  chroms <- c(as.character(1:22), "X", "Y")
  allRes <- do.call(rbind, lapply(chroms, function(chrom){
    parquet <- arrow::read_parquet(paste0(fpath, ".", chrom, ".parquet"))
    return(parquet)
  }))
  padj <- p.adjust(allRes$pval_nominal, method = "fdr")
  parquetSig <- allRes[which(padj < cutoff),]
  cat(".")
  return(parquetSig)
}
cisResults <- lapply(1:10, function(i){
  cisResultsSplit <- lapply(1:2, function(j){
    return(readAllCis(paste0(eQTLDir, "/profoundAutismBoth_rand_", i, "_", j, "_diff.cis_qtl_pairs")))
  })
  names(cisResultsSplit) <- 1:2
  return(cisResultsSplit)
})
names(cisResults) <- 1:10

# We use the number of phenotypes x number of variants as the upper bound for M (the number of tests).
# Technically it is this number - the number of cis pairs filtered out, but that number
# is negligible in comparison to the number of tests.
M <- 34070 * 2440283
getAboveCutoff <- function(file){
  ranks <- rank(file$pval, ties.method = "first")
  padj <- file$pval * (M / ranks)
  sig <- file[which(padj < cutoff),]
  return(sig)
}
transResults <- lapply(1:10, function(i){
  transResultsSplit <- lapply(1:2, function(j){
    cat("*")
    return(getAboveCutoff(file = arrow::read_parquet(paste0(eQTLDir, "/profoundAutismBoth_rand_", i, "_", j, "_diff_trans.trans_qtl_pairs.parquet"))))
  })
  names(transResultsSplit) <- 1:2
  return(transResultsSplit)
})
names(transResults) <- 1:10

# Look at overlap in pairs.
getPairs <- function(result){
  return(paste(result$variant_id, result$phenotype_id, sep = "__"))
}
pairsCis <- lapply(1:10, function(i){
  pairsSplit <- lapply(1:2, function(j){
    return(getPairs(cisResults[[i]][[j]]))
  })
  names(pairsSplit) <- 1:2
  return(pairsSplit)
})
names(pairsCis) <- 1:10
pairsTrans <- lapply(1:10, function(i){
  pairsSplit <- lapply(1:2, function(j){
    return(getPairs(transResults[[i]][[j]]))
  })
  names(pairsSplit) <- 1:2
  return(pairsSplit)
})
names(pairsTrans) <- 1:10

saveRDS(list(pairsCis = pairsCis,
             pairsTrans = pairsTrans), paste0(eQTLDir, "/alleQTLResultsRandDiff.RDS"))