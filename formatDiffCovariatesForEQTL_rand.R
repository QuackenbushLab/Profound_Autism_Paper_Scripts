for(i in 1:10){
  for(j in 1:2){
    expressionBed <- paste0("../eQTL/expression.pheno.profoundBoth_", i, "_", j, "_diff.bed")
    
    # Filter to include only the samples also in the BED file.
    bed <- read.table(expressionBed, sep = "\t",
                      comment.char = "", header = TRUE, check.names = FALSE)
    covarNew <- read.table(paste0("../eQTL/covarNewProfoundBoth_", i, "_", j, ".txt"), 
                           sep = "\t", header = FALSE, row.names = 1)
    print(rownames(covarNew))
    colnames(covarNew) <- unname(unlist(covarNew[1,]))
    shared <- intersect(colnames(covarNew), colnames(bed))
    str(covarNew)
    covarNew <- covarNew[,shared]
    bed <- bed[,c("#chr", "start", "end", "phenotype_id", shared)]
    bed[,"#chr"] <- unlist(lapply(bed[,"#chr"], function(chr){
      return(strsplit(chr, "chr")[[1]][2])
    }))
    print(rownames(covarNew))
    str(covarNew)
    write.table(covarNew, paste0("../eQTL/covarNewProfoundBoth_", i, "_", j, "_diff.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
    write.table(bed, paste0("../eQTL/expression.pheno.profoundBoth_", i, "_", j, "_diff_filt.bed"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}