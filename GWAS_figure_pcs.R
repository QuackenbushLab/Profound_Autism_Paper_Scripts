sourceDir <- NULL
refFile <- NULL
# Read in the reference.
library(data.table)
ref <- fread(refFile)

# Read in the GWAS results.
profoundAll <- read.csv(paste0(sourceDir, "profoundAll_noSex_pcs.PHENO1.glm.logistic.hybrid.sig"), row.names = 1)
profoundBoth <- read.csv(paste0(sourceDir, "profoundBoth_noSex_pcs.PHENO1.glm.logistic.hybrid.sig"), row.names = 1)
profoundBothID <- read.csv(paste0(sourceDir, "profoundBothID_noSex_pcs.PHENO1.glm.logistic.hybrid.sig"), row.names = 1)
profoundBothNonverbal <- read.csv(paste0(sourceDir, "profoundBothNonverbal_noSex_pcs.PHENO1.glm.logistic.hybrid.sig"), row.names = 1)
verbalAll <- read.csv(paste0(sourceDir, "notProfoundAll_noSex_pcs.PHENO1.glm.logistic.hybrid.sig"), row.names = 1)
verbalNoGifted <- read.csv(paste0(sourceDir, "notProfoundNoGifted_noSex_pcs.PHENO1.glm.logistic.hybrid.sig"), row.names = 1)
verbalNoID <- read.csv(paste0(sourceDir, "notProfoundNoID_noSex_pcs.PHENO1.glm.logistic.hybrid.sig"), row.names = 1)

# Make a consolidated table of SNPs.
snps <- unique(c(verbalAll$ID, verbalNoGifted$ID, verbalNoID$ID,
                 profoundAll$ID, profoundBothNonverbal$ID, profoundBothID$ID))
snpMapping <- do.call(rbind, lapply(snps, function(snp){
  
  # Initialize.
  whereSig <- data.frame(nonProfoundAll = 0, nonProfoundExclGifted = 0, nonProfoundExclID = 0,
                         profoundAll = 0, profoundExclIDOnly = 0, profoundExclNonverbalOnly = 0)
  
  foundOnce <- FALSE
  dat <- NULL
  
  # Build data frame with info.
  if(snp %in% verbalAll$ID){
    whereSig[1,"nonProfoundAll"] <- 1
    dat <- verbalAll[which(verbalAll$ID == snp),1:4]
    foundOnce <- TRUE
  }
  if(snp %in% verbalNoGifted$ID){
    whereSig[1,"nonProfoundExclGifted"] <- 1
    if(foundOnce == FALSE){
      dat <- verbalNoGifted[which(verbalNoGifted$ID == snp),1:4]
      foundOnce <- TRUE
    }
  }
  if(snp %in% verbalNoID$ID){
    whereSig[1,"nonProfoundExclID"] <- 1
    if(foundOnce == FALSE){
      dat <- verbalNoID[which(verbalNoID$ID == snp),1:4]
      foundOnce <- TRUE
    }
  }
  if(snp %in% profoundAll$ID){
    whereSig[1,"profoundAll"] <- 1
    if(foundOnce == FALSE){
      dat <- profoundAll[which(profoundAll$ID == snp),1:4]
      foundOnce <- TRUE
    }
  }
  if(snp %in% profoundBothID$ID){
    whereSig[1,"profoundExclIDOnly"] <- 1
    if(foundOnce == FALSE){
      dat <- profoundBothID[which(profoundBothID$ID == snp),1:4]
      foundOnce <- TRUE
    }
  }
  if(snp %in% profoundBothNonverbal$ID){
    whereSig[1,"profoundExclNonverbalOnly"] <- 1
    if(foundOnce == FALSE){
      dat <- profoundBothNonverbal[which(profoundBothNonverbal$ID == snp),1:4]
      foundOnce <- TRUE
    }
  }
  
  # Add list of tests.
  finalDat <- cbind(dat, whereSig)
  return(finalDat)
}))

# Map each SNP.
mapSNPs <- do.call(rbind, lapply(1:nrow(snpMapping), function(i){
  
  # Subset chromosome.
  chrName <- paste0("chr", as.character(snpMapping[i,"X.CHROM"]))
  chrNewName <- paste0("chr", as.character(snpMapping[i,"X.CHROM"]))
  if(snpMapping[i,"X.CHROM"] == "XY"){
    chrName <- "chrX"
    chrNewName <- "chrXY"
  }
  snpMapping[i,"X.CHROM"] <- chrNewName
  refChrom <- ref[which(ref$V1 == chrName),]
  
  # Find starting points where SNP position is greater.
  refGreater <- refChrom[which(refChrom$V2 <= snpMapping[i,"POS"])]

  # Find ending points where SNP position is less.
  refLess <- refGreater[which(refGreater$V3 >= snpMapping[i,"POS"])]

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
    str(snpMapping[i,"POS"])
    str(refChrom)
    closestStart <- refChrom[which.min(abs(snpMapping[i,"POS"] - refChrom$V2)),]
    closestEnd <- refChrom[which.min(abs(refChrom$V3 - snpMapping[i,"POS"])),]
    refImportant <- closestStart
    if(abs(snpMapping[i,"POS"] - closestStart$V2) < abs(closestEnd$V3 - snpMapping[i,"POS"])){
      refImportant <- closestEnd
    }
  }
  # There are two cases where this happens. In each case, we want the first one.
  if(nrow(refImportant) == 2){
    refImportant <- refImportant[1,]
  }
  
  # Paste together the SNP and gene information.
  snpGeneInfo <- cbind(snpMapping[i,], refImportant)
  snpDetails <- snpGeneInfo[,ncol(snpGeneInfo)]
  snpGeneInfo <- snpGeneInfo[,c(1:10, 12:14, 16, 18)]
  colnames(snpGeneInfo)[1] <- "chrom"
  colnames(snpGeneInfo)[11:15] <- c("chromStart", "chromEnd", "name",
                                   "strand", "region")
  splitDetails <- strsplit(snpDetails, "; ")[[1]]
  geneType <- splitDetails[startsWith(splitDetails, "gene_type")]
  geneTypeSplit <- strsplit(geneType, split = "\"", fixed = TRUE)[[1]][2]
  geneName <- splitDetails[startsWith(splitDetails, "gene_name")]
  geneNameSplit <- strsplit(geneName, split = "\"", fixed = TRUE)[[1]][2]
  snpGeneInfo$geneName <- geneNameSplit
  snpGeneInfo$geneType <- geneTypeSplit
  return(snpGeneInfo)
}))
write.table(mapSNPs, paste0(sourceDir, "snpsAnnotated_noSex_pcs.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

# Make an UpSet plot.
library(ComplexHeatmap)
mapSNPsMat <- mapSNPs[,5:10]
colnames(mapSNPsMat) <- c("ASD Composite", "ASD Not Gifted",
                          "ASD No ID", "PA Composite",
                          "PA Nonspeaking", "PA Moderate-Profound ID")
mapSNPComb <- make_comb_mat(mapSNPsMat, mode = "distinct")

tanno = upset_top_annotation(mapSNPComb, 
                                      add_numbers = TRUE)
UpSet(mapSNPComb, top_annotation = tanno,
      set_order = c("ASD Composite", "ASD Not Gifted",
                                "ASD No ID", "PA Composite",
                                "PA Nonspeaking", 
                                "PA Moderate-Profound ID")) -> ht
draw(ht, padding = unit(c(5, 10, 5, 5), "mm"))