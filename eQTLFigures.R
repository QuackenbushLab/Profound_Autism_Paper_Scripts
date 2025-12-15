library("ComplexHeatmap")

# Read eQTL data.
eQTLDir <- "../eQTL/"
eQTLRawComb <- readRDS(paste0(eQTLDir, "eQTLMat.csv_comb.RDS"))
eQTLDiff <- readRDS(paste0(eQTLDir, "eQTLDiffMat.csv_comb.RDS"))
combs <- normalize_comb_mat(list(eQTLRawComb, eQTLDiff))

# Subset.
subset1 = combs[[1]][, 17:26]
subset2 = combs[[2]][, 17:26]

upsetRaw <- UpSet(subset1,
               top_annotation = upset_top_annotation(subset1, add_numbers = TRUE),
               right_annotation = upset_right_annotation(combs[[1]], add_numbers = TRUE))
upsetRaw@name = "rawExpression"
upsetDiff <- UpSet(subset2,
                   top_annotation = upset_top_annotation(subset2, add_numbers = TRUE),
                   right_annotation = upset_right_annotation(combs[[2]], add_numbers = TRUE))
upsetDiff@row_order <- upsetRaw@row_order
upsetDiff@column_order <- upsetRaw@column_order
upsetDiff@name = "diffExpression"
str(upsetRaw)
str(upsetDiff)
par(mfrow = c(2, 1))
pdf(paste0(eQTLDir, "eQTLCombinedUpsets.pdf"), width = 8, height = 6)
draw(upsetRaw %v% upsetDiff, padding = unit(c(0, 0, 0, 1), "cm"))
dev.off()