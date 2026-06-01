library("ComplexHeatmap")

# Read eQTL data.
eQTLDir <- "../eQTL/"
eQTLRawComb <- readRDS(paste0(eQTLDir, "eQTLMatCis.csv_comb.RDS"))
eQTLDiff <- readRDS(paste0(eQTLDir, "eQTLDiffMatCis.csv_comb.RDS"))
combs <- normalize_comb_mat(list(eQTLRawComb, eQTLDiff))

# Subset.
subset1 = combs[[1]][, 13:22]
subset2 = combs[[2]][, 13:22]

upsetRaw <- UpSet(subset1,
                  top_annotation = upset_top_annotation(subset1, add_numbers = TRUE),
                  right_annotation = upset_right_annotation(combs[[1]], add_numbers = TRUE))
upsetRaw@name = "logCPM"
upsetDiff <- UpSet(subset2,
                   top_annotation = upset_top_annotation(subset2, add_numbers = TRUE),
                   right_annotation = upset_right_annotation(combs[[2]], add_numbers = TRUE))
upsetDiff@row_order <- upsetRaw@row_order
upsetDiff@column_order <- upsetRaw@column_order
upsetDiff@name = "diffExpression"
str(upsetRaw)
str(upsetDiff)
par(mfrow = c(2, 1))
pdf(paste0(eQTLDir, "eQTLCombinedUpsetsCis2.pdf"), width = 8, height = 6)
draw(upsetRaw %v% upsetDiff, padding = unit(c(0, 0, 0, 1), "cm"))
dev.off()

# Do trans data.
eQTLRawComb <- readRDS(paste0(eQTLDir, "eQTLMatTrans.csv_comb.RDS"))
eQTLDiff <- readRDS(paste0(eQTLDir, "eQTLDiffMatTrans.csv_comb.RDS"))
combs <- normalize_comb_mat(list(eQTLRawComb, eQTLDiff))

# Subset.
subset1 = combs[[1]][, 17:26]
subset2 = combs[[2]][, 17:26]

upsetRaw <- UpSet(subset1,
                  top_annotation = upset_top_annotation(subset1, add_numbers = TRUE),
                  right_annotation = upset_right_annotation(combs[[1]], add_numbers = TRUE))
upsetRaw@name = "logCPM"
upsetDiff <- UpSet(subset2,
                   top_annotation = upset_top_annotation(subset2, add_numbers = TRUE),
                   right_annotation = upset_right_annotation(combs[[2]], add_numbers = TRUE))
upsetDiff@row_order <- upsetRaw@row_order
upsetDiff@column_order <- upsetRaw@column_order
upsetDiff@name = "diffExpression"
str(upsetRaw)
str(upsetDiff)
par(mfrow = c(2, 1))
pdf(paste0(eQTLDir, "eQTLCombinedUpsetsTrans.pdf"), width = 8, height = 6)
draw(upsetRaw %v% upsetDiff, padding = unit(c(0, 0, 0, 1), "cm"))
dev.off()
