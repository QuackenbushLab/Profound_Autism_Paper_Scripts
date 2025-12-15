# Calculate and save percent variance.
inDir <- "/n/holylabs/LABS/quackenbush_lab/Lab/teicher/PLINK/resultFiles"
eigenvals <- read.csv(paste0(inDir, "/omniPCA.eigenval"), header = FALSE)
eigenvals$percentVariance <- (eigenvals$V1 / sum(eigenvals$V1)) * 100
write.csv(eigenvals, paste0(inDir, "/omniPCA.eigenval.withVariance.csv"))

# Plot the percentages.
pdf(paste0(inDir, "/omniPCA.eigenval.percentVariance.pdf"))
plot(x = 1:50, y = eigenvals$percentVariance[1:50], xlab = "Principal Component",
     ylab = "Percent Variance", pch = 19, type = "b")
dev.off()