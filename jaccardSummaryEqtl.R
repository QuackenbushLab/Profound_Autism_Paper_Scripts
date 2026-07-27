library(biomaRt)

eqtlDir <- NULL
jaccardLogCPM <- readRDS(paste0(eqtlDir, "jaccardRandLogCPM.RDS"))
jaccardDiff <- readRDS(paste0(eqtlDir, "jaccardRandDiff.RDS"))

# Get summary stats.
print(paste("Cis, logCPM:", mean(jaccardLogCPM$jaccardCis), sd(jaccardLogCPM$jaccardCis)))
print(paste("Trans, logCPM:", mean(jaccardLogCPM$jaccardTrans), sd(jaccardLogCPM$jaccardTrans)))
print(paste("Cis, diff:", mean(jaccardDiff$jaccardCis), sd(jaccardDiff$jaccardCis)))
print(paste("Trans, diff:", mean(jaccardDiff$jaccardTrans), sd(jaccardDiff$jaccardTrans)))

# Find shared eQTLs.
findShared <- function(overlappingPairs){
  longList <- unlist(overlappingPairs)
  counts <- table(longList)
  shared <- names(counts)[which(counts == 10)]
  return(shared)
}
sharedLogCPMCis <- findShared(jaccardLogCPM$overlappingPairsCis)
sharedLogCPMTrans <- findShared(jaccardLogCPM$overlappingPairsTrans)
sharedDiffCis <- findShared(jaccardDiff$overlappingPairsCis)
sharedDiffTrans <- findShared(jaccardDiff$overlappingPairsTrans)

# Read in the reference.
refFile <- NULL
ref <- as.data.frame(fread(refFile))
snpFile <- NULL
snpPvar <- as.data.frame(fread(snpFile))
rownames(snpPvar) <- snpPvar$ID

# Map each SNP.
#snpMart <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
#attributes <- c("chr_name", "chrom_start", "clinical_significance", "associated_gene",
#                "distance_to_transcript", "allele")

mapSNPs <- function(snps){
  #snpInfo <- getBM(filters = "snp_filter", values = snps, mart = snpMart,
  #                 attributes = attributes)
  snpInfo <- snpPvar[snps,]
  colnames(snpInfo)[1] <- "CHROM"
  snpGene <- do.call(rbind, lapply(1:nrow(snpInfo), function(i){

    # Subset chromosome.
    chrName <- paste0("chr", as.character(snpInfo[i,"CHROM"]))
    if(snpInfo[i,"CHROM"] == 24){
      chrName <- "chrX"
    }else if(snpInfo[i,"CHROM"] == 25){
      chrName <- "chrY"
    }
    refChrom <- ref[which(ref$V1 == chrName),]

    # Find starting points where SNP position is greater.
    refGreater <- refChrom[which(refChrom$V2 <= snpInfo[i,"POS"]),]

    # Find ending points where SNP position is less.
    refLess <- refGreater[which(refGreater$V3 >= snpInfo[i,"POS"]),]

    # Order of importance.
    refImportant <- refLess
    ooi <- c("CDS", "start_codon", "stop_codon", "Selenocysteine", "exon", "gene", "transcript",
             "UTR")
    if(nrow(refImportant) > 0){
      found <- FALSE
      j = 1
      while(found == FALSE && j < length(ooi)){
        if(ooi[j] %in% refLess$V8){
          refImportant <- refLess[which(refLess$V8 == ooi[j]),]
          found <- TRUE
        }else{
          j <- j + 1
        }
      }
    }else{
      # If the variant doesn't map to a gene, then map it to its closest gene.
      closestStart <- refChrom[which.min(abs(snpInfo[i,"POS"] - refChrom$V2)),]
      closestEnd <- refChrom[which.min(abs(refChrom$V3 - snpInfo[i,"POS"])),]
      refImportant <- closestStart
      if(abs(snpInfo[i,"POS"] - closestStart$V2) < abs(closestEnd$V3 - snpInfo[i,"POS"])){
        refImportant <- closestEnd
      }
    }
    snpDetails <- refImportant[,ncol(refImportant)]
    if(nrow(refImportant) > 1){
      refImportant <- refImportant[1,]
    }

    refImportant <- refImportant[,c(2:4, 6, 8)]
    colnames(refImportant) <- c("geneStart", "geneEnd", "name", "strand", "region")
  
    # Paste together the SNP and gene information.
    snpGeneInfo <- cbind(snpInfo[i,], refImportant)
    splitDetails <- strsplit(snpDetails, "; ")[[1]]
    geneType <- splitDetails[startsWith(splitDetails, "gene_type")]
    geneTypeSplit <- strsplit(geneType, split = "\"", fixed = TRUE)[[1]][2]
    geneName <- splitDetails[startsWith(splitDetails, "gene_name")]
    geneNameSplit <- strsplit(geneName, split = "\"", fixed = TRUE)[[1]][2]
    snpGeneInfo$geneName <- geneNameSplit
    snpGeneInfo$geneType <- geneTypeSplit
    return(snpGeneInfo)
  }))
  return(snpGene)
}
mapGenes <- function(genes){
  refGenesList <- lapply(genes, function(gene){
    refGene <- ref[which(ref$V4 == gene)[1], c(1:4, 6, 8)]
    colnames(refGene) <- c("expChrome", "expStart", "expEnd", "expId", "expStrand", "expRegion")
    geneDetails <- ref[which(ref$V4 == gene)[1], 10]
    splitDetails <- strsplit(geneDetails, "; ")[[1]]
    geneType <- splitDetails[startsWith(splitDetails, "gene_type")]
    geneTypeSplit <- strsplit(geneType, split = "\"", fixed = TRUE)[[1]][2]
    geneName <- splitDetails[startsWith(splitDetails, "gene_name")]
    geneNameSplit <- strsplit(geneName, split = "\"", fixed = TRUE)[[1]][2]
    refGene$expGeneName <- geneNameSplit
    refGene$expGeneType <- geneTypeSplit
    return(refGene)
  })
  refGenes <- do.call(rbind, refGenesList)
  str(refGenes)
  return(refGenes)
}
makePairTable <- function(pairs, fileName){
  snps <- unlist(lapply(pairs, function(eqtl){return(strsplit(eqtl, "__")[[1]][1])}))
  genes <- unlist(lapply(pairs, function(eqtl){return(strsplit(eqtl, "__")[[1]][2])}))
  mapping <- mapSNPs(snps)
  geneExtra <- mapGenes(genes)
  eqtlInfo <- cbind(mapping, geneExtra)
  write.csv(eqtlInfo, paste0(eqtlDir, fileName))
  return(eqtlInfo)
}
logCPMCisTable <- makePairTable(sharedLogCPMCis, "overlappingPairsLogCPMCis.csv")
logCPMTransTable <- makePairTable(sharedLogCPMTrans, "overlappingPairsLogCPMTrans.csv")
diffTransTable <- makePairTable(sharedDiffTrans, "overlappingPairsDiffTrans.csv")
fullTable <- data.frame(expressionType = c(rep("logCPM", nrow(logCPMCisTable) + nrow(logCPMTransTable)),
                                           rep("diff", nrow(diffTransTable))),
                        eqtlType = c(rep("cis", nrow(logCPMCisTable)),
                                     rep("trans", nrow(logCPMTransTable) + nrow(diffTransTable))))
fullTableAdd <- do.call(rbind, list(logCPMCisTable, logCPMTransTable, diffTransTable))
fullTable <- cbind(fullTable, fullTableAdd)
write.csv(fullTable, paste0(eQTLDir, "eqtlTable.csv"))