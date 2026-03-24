# Read in the profound autism (both) file.
paDir <- "/n/holylabs/LABS/quackenbush_lab/Lab/teicher/profoundAutism/"
eqtlDir <- "/n/holylabs/LABS/quackenbush_lab/Lab/teicher/eQTL/"
paData <- read.csv(paste0(paDir, "profoundAutismBoth_above8.csv"), row.names = 1)
paSamps <- row.names(paData)
merged <- read.table(paste0(eqtlDir, "merged.txt"), sep = "\t", header = TRUE)

# Generate random subsets and save.
set.seed(1)
for(i in 1:10){
  
  # Save the phenotypic data.
  subset1 <- sample(paSamps, floor(length(paSamps) / 2))
  subset2 <- setdiff(paSamps, subset1)
  write.csv(paData[subset1,], paste0(paDir, "profoundAutismBoth_rand_", i, "_1.csv"))
  write.csv(paData[subset2,], paste0(paDir, "profoundAutismBoth_rand_", i, "_2.csv"))
  
  # Save sample data for PLINK.
  whichHasSNP1 <- which(merged$individual_id %in% subset1)
  whichHasSNP2 <- which(merged$individual_id %in% subset2)
  subset1PLINK <- merged[whichHasSNP1, c("family_id", "array_id")]
  subset2PLINK <- merged[whichHasSNP2, c("family_id", "array_id")]
  str(subset1PLINK)
  str(subset2PLINK)
  write.table(subset1PLINK, paste0(paDir, "profoundAutismBoth_rand_samps_", i, "_1.csv"),
              quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")
  write.table(subset2PLINK, paste0(paDir, "profoundAutismBoth_rand_samps_", i, "_2.csv"),
              quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")
}