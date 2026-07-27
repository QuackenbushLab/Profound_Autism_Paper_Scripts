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
profoundBothCis <- readAllCis(paste0(eQTLDir, "/profoundBothDiffResults.cis_qtl_pairs"))
profoundModerateIDCis <- readAllCis(paste0(eQTLDir, "/profoundModerateIDDiffResults.cis_qtl_pairs"))
verbalGiftedCis <- readAllCis(paste0(eQTLDir, "/verbalGiftedDiffResults.cis_qtl_pairs"))
verbalNoIDCis <- readAllCis(paste0(eQTLDir, "/verbalNoIDDiffResults.cis_qtl_pairs"))
verbalMildCis <- readAllCis(paste0(eQTLDir, "/verbalMildIDDiffResults.cis_qtl_pairs"))

M <- 34070 * 2440283
getAboveCutoff <- function(file){
  ranks <- rank(file$pval, ties.method = "first")
  padj <- file$pval * (M / ranks)
  sig <- file[which(padj < cutoff),]
  str(sig)
  return(sig)
}
profoundBothSig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/profoundBothDiffTransResults.trans_qtl_pairs.parquet")))
profoundModerateIDOnlySig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/profoundModerateIDDiffTransResults.trans_qtl_pairs.parquet")))
verbalGiftedTransSig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/verbalGiftedDiffTransResults.trans_qtl_pairs.parquet")))
verbalNoIDTransSig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/verbalNoIDDiffTransResults.trans_qtl_pairs.parquet")))
verbalMildIDTransSig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/verbalMildIDDiffTransResults.trans_qtl_pairs.parquet")))

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
           pairsMildIDTransSig = pairsMildIDTransSig), paste0(eQTLDir, "/alleQTLDiffResults.RDS"))