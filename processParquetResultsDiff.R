if(!require("arrow")){
  install.packages("arrow")
}
library("arrow")

# Read files.
eQTLDir <- "../eQTL"
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
#profoundCis <- readAllCis(paste0(eQTLDir, "/profoundAllResults.cis_qtl_pairs"))
#notProfoundCis <- readAllCis(paste0(eQTLDir, "/notProfoundAllResults.cis_qtl_pairs"))
profoundBothCis <- readAllCis(paste0(eQTLDir, "/profoundBothDiffResults.cis_qtl_pairs"))
profoundModerateIDCis <- readAllCis(paste0(eQTLDir, "/profoundModerateIDDiffResults.cis_qtl_pairs"))
profoundNonverbalCis <- readAllCis(paste0(eQTLDir, "/profoundNonverbalDiffResults.cis_qtl_pairs"))
verbalGiftedCis <- readAllCis(paste0(eQTLDir, "/verbalGiftedDiffResults.cis_qtl_pairs"))
verbalNoIDCis <- readAllCis(paste0(eQTLDir, "/verbalNoIDDiffResults.cis_qtl_pairs"))
verbalMildCis <- readAllCis(paste0(eQTLDir, "/verbalMildIDDiffResults.cis_qtl_pairs"))

getAboveCutoff <- function(file){
  padj <- p.adjust(file$pval, method = "fdr")
  return(file[which(padj < cutoff),])
}
#profoundTransSig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/profoundAllResultsTrans.trans_qtl_pairs.parquet")))
profoundBothSig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/profoundBothDiffTransResults.trans_qtl_pairs.parquet")))
profoundNonverbalOnlySig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/profoundNonverbalDiffTransResults.trans_qtl_pairs.parquet")))
profoundModerateIDOnlySig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/profoundModerateIDDiffTransResults.trans_qtl_pairs.parquet")))
#notProfoundTransSig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/notProfoundAllResultsTrans.trans_qtl_pairs.parquet")))
verbalGiftedTransSig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/verbalGiftedDiffTransResults.trans_qtl_pairs.parquet")))
verbalNoIDTransSig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/verbalNoIDDiffTransResults.trans_qtl_pairs.parquet")))
verbalMildIDTransSig <- getAboveCutoff(arrow::read_parquet(paste0(eQTLDir, "/verbalMildIDDiffTransResults.trans_qtl_pairs.parquet")))

# Look at overlap in pairs.
getPairs <- function(result){
  return(paste(result$variant_id, result$phenotype_id, sep = "__"))
}
#pairsProfoundCis <- getPairs(profoundCis)
#pairsNotProfoundCis <- getPairs(notProfoundCis)
#pairsProfoundTrans <- getPairs(profoundTransSig)
#pairsNotProfoundTrans <- getPairs(notProfoundTransSig)

pairsProfoundBothCis <- getPairs(profoundBothCis)
pairsProfoundNonverbalCis <- getPairs(profoundNonverbalCis)
pairsProfoundModerateIDCis <- getPairs(profoundModerateIDCis)
pairsGiftedCis <- getPairs(verbalGiftedCis)
pairsNoIDCis <- getPairs(verbalNoIDCis)
pairsMildIDCis <- getPairs(verbalMildCis)

pairsProfoundBothTransSig <- getPairs(profoundBothSig)
pairsProfoundNonverbalTransSig <- getPairs(profoundNonverbalOnlySig)
pairsProfoundModerateIDTransSig <- getPairs(profoundModerateIDOnlySig)
pairsGiftedTransSig <- getPairs(verbalGiftedTransSig)
pairsNoIDTransSig <- getPairs(verbalNoIDTransSig)
pairsMildIDTransSig <- getPairs(verbalMildIDTransSig)

saveRDS(list(#pairsProfoundCis = pairsProfoundCis,
             #pairsNotProfoundCis = pairsNotProfoundCis,
             #pairsProfoundTrans = pairsProfoundTrans,
             #pairsNotProfoundTrans = pairsNotProfoundTrans,
             pairsProfoundBothCis = pairsProfoundBothCis,
             pairsProfoundNonverbalCis = pairsProfoundNonverbalCis,
             pairsProfoundModerateIDCis = pairsProfoundModerateIDCis,
             pairsGiftedCis = pairsGiftedCis,
             pairsNoIDCis = pairsNoIDCis,
             pairsMildIDCis = pairsMildIDCis,
             pairsProfoundBothTransSig = pairsProfoundBothTransSig,
             pairsProfoundNonverbalTransSig = pairsProfoundNonverbalTransSig,
             pairsProfoundModerateIDTransSig = pairsProfoundModerateIDTransSig,
             pairsGiftedTransSig = pairsGiftedTransSig,
             pairsNoIDTransSig = pairsNoIDTransSig,
           pairsMildIDTransSig = pairsMildIDTransSig), paste0(eQTLDir, "/alleQTLDiffResults.RDS"))