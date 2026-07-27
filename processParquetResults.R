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
  str(allRes)
  padj <- p.adjust(allRes$pval_nominal, method = "fdr")
  parquetSig <- allRes[which(padj < cutoff),]
  str(parquetSig)
  return(parquetSig)
}
profoundBothCis <- readAllCis(paste0(eQTLDir, "/profoundBoth.cis_qtl_pairs"))
profoundModerateIDCis <- readAllCis(paste0(eQTLDir, "/profoundModerateID.cis_qtl_pairs"))
verbalGiftedCis <- readAllCis(paste0(eQTLDir, "/verbalGifted.cis_qtl_pairs"))
verbalNoIDCis <- readAllCis(paste0(eQTLDir, "/verbalNoID.cis_qtl_pairs"))
verbalMildCis <- readAllCis(paste0(eQTLDir, "/verbalMildID.cis_qtl_pairs"))

# We use the number of phenotypes x number of variants as the upper bound for M (the number of tests).
# Technically it is this number - the number of cis pairs filtered out, but that number
# is negligible in comparison to the number of tests.
M <- 34070 * 2440283
getAboveCutoff <- function(file){
  ranks <- rank(file$pval, ties.method = "first")
  padj <- file$pval * (M / ranks)
  sig <- file[which(padj < cutoff),]
  str(sig)
  return(sig)
}
profoundBothSig <- getAboveCutoff(file = arrow::read_parquet(paste0(eQTLDir, "/profoundBoth.trans_qtl_pairs.parquet")))
profoundModerateIDOnlySig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/profoundModerateID.trans_qtl_pairs.parquet")))
verbalGiftedTransSig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/verbalGifted.trans_qtl_pairs.parquet")))
verbalNoIDTransSig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/verbalNoID.trans_qtl_pairs.parquet")))
verbalMildIDTransSig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/verbalMildID.trans_qtl_pairs.parquet")))

# Look at overlap in pairs.
getPairs <- function(result){
  return(paste(result$variant_id, result$phenotype_id, sep = "__"))
}

pairsProfoundBothCis <- getPairs(profoundBothCis)
pairsProfoundModerateIDCis <- getPairs(profoundModerateIDCis)
pairsGiftedCis <- getPairs(verbalGiftedCis)
pairsNoIDCis <- getPairs(verbalNoIDCis)
pairsMildIDCis <- getPairs(verbalMildCis)

pairsProfoundBothTransSig <- getPairs(profoundBothSig)
pairsProfoundModerateIDTransSig <- getPairs(profoundModerateIDOnlySig)
pairsGiftedTransSig <- getPairs(verbalGiftedTransSig)
pairsNoIDTransSig <- getPairs(verbalNoIDTransSig)
pairsMildIDTransSig <- getPairs(verbalMildIDTransSig)

saveRDS(list(pairsProfoundBothCis = pairsProfoundBothCis,
             pairsProfoundModerateIDCis = pairsProfoundModerateIDCis,
             pairsGiftedCis = pairsGiftedCis,
             pairsNoIDCis = pairsNoIDCis,
             pairsMildIDCis = pairsMildIDCis,
             pairsProfoundBothTransSig = pairsProfoundBothTransSig,
             pairsProfoundModerateIDTransSig = pairsProfoundModerateIDTransSig,
             pairsGiftedTransSig = pairsGiftedTransSig,
             pairsNoIDTransSig = pairsNoIDTransSig,
           pairsMildIDTransSig = pairsMildIDTransSig), paste0(eQTLDir, "/alleQTLResults.RDS"))