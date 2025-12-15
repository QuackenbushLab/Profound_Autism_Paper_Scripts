# Read the data.
sourceDir <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/SFARI/"
large <- read.csv(paste0(sourceDir, "geneExpressionLogCPMLarge.csv"),
                  row.names = 1)

# Split data into families.
familyInfo <- do.call(rbind, lapply(colnames(large), function(samp){
  
  # Get the family label. Remove the "X" added by R.
  sampSplitDot <- strsplit(samp, ".", fixed = TRUE)[[1]]
  family <- substr(sampSplitDot[1], 2, nchar(sampSplitDot[1]))

  # Get the position within the family.
  positionStr <- substr(sampSplitDot[2], 1, 3)
  
  # Construct a data frame.
  retval <- data.frame(family = family, position = positionStr)
  return(retval)
}))

# Extract only probands and siblings.
whichChildren <- union(which(familyInfo$position == "s1"),
                       which(familyInfo$position == "p1"))
children <- large[,whichChildren]
childrenInfo <- familyInfo[whichChildren,]
write.csv(children, paste0(sourceDir, "childrenGeneExpression.csv"))
write.csv(childrenInfo, paste0(sourceDir, "childrenLabels.csv"), row.names = FALSE)

# z-scale each gene.
geneMeans <- rowMeans(children)
geneSds <- unlist(lapply(rownames(children), function(gene){
    return(sd(children[gene,]))
}))
childrenScaled <- (children - geneMeans) / geneSds

# For each family, obtain differences.
# Note that we had to exclude 2 families because there was no sibling data.
families <- sort(unique(childrenInfo$family))
diffExpressionList <- lapply(families, function(family){
  
  # Initialize differential expression.
  diffExpression <- NULL
  
  # Extract proband and sibling data.
  whichProband <- intersect(which(childrenInfo$family == family),
                            which(childrenInfo$position == "p1"))
  whichSibling <- intersect(which(childrenInfo$family == family),
                            which(childrenInfo$position == "s1"))
  if(length(whichSibling) == 0){
    print(paste("No sibling for family", family, "-- skipping"))
  }else if(length(whichProband) == 0){
    print(paste("No proband for family", family, "-- skipping"))
  }else{
    proband <- childrenScaled[,whichProband]
    sibling <- childrenScaled[,whichSibling]
    names(proband) <- rownames(childrenScaled)
    names(sibling) <- rownames(childrenScaled)
    
    # Compute differential expression.
    diffExpression <- data.frame(v1 = proband - sibling)
    rownames(diffExpression) <- names(proband)
    colnames(diffExpression) <- family
  }
  
  return(diffExpression)
})
isNull <- unlist(lapply(diffExpressionList,is.null))
diffExpressionList <- diffExpressionList[which(isNull == FALSE)]
diffExpression <- do.call(cbind, diffExpressionList)

# Save the differential expression matrix.
diffExpressionMat <- as.matrix(diffExpression)
write.csv(diffExpressionMat, paste0(sourceDir, "diffExpression.csv"))