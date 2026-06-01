# ---- user inputs ----
expr_csv   <- NULL  # genes x samples
gtf_path   <- NULL                            # e.g., GENCODE/GTF
geno_prefix <- NULL
use_tss    <- FALSE   # TRUE = 1bp TSS; FALSE = full gene body
keep_chr_prefix <- FALSE  # set TRUE if your VCF/PLINK uses "chr1" etc.

# ---- packages ----
suppressPackageStartupMessages({
  library(rtracklayer)  # to read GTF
})

# ---- 1) read expression (genes in rows, samples in columns) ----
expr <- read.csv(expr_csv, row.names = 1, check.names = FALSE)
# ensure numeric
expr[] <- lapply(expr, function(x) suppressWarnings(as.numeric(as.character(x))))
expr <- as.data.frame(expr)
print("Read")

# Map samples.
map <- read.table("../eQTL/merged.txt", sep = "\t", header = TRUE)
str(map)
iid <- unlist(lapply(colnames(expr), function(samp){
  return(paste0(samp, ".p1"))
}))
aid <- unlist(lapply(iid, function(id){
  retval <- NA
  which_iid <- which(map$individual_id == id)
  aid <- map[which_iid, "array_id"]
  if(length(which_iid) == 1){
          retval <-  aid
  }
  return(retval)
}))
str(expr)
expr <- expr[,which(!is.na(aid))]
str(expr)
aid <- aid[which(!is.na(aid))]
colnames(expr) <- aid
str(expr)

# drop genes with all-NA
expr <- expr[rowSums(!is.na(expr)) > 0, , drop = FALSE]
print("Dropped")

# ---- 2) gene coordinates from GTF ----
gr <- import(gtf_path)
gr_gene <- gr[gr$type == "gene"]
# gene_id and (optionally) gene_name
gid <- if (!is.null(mcols(gr_gene)$gene_id)) mcols(gr_gene)$gene_id else mcols(gr_gene)$ID
gname <- mcols(gr_gene)$gene_name
print("coords")

# Retain only the identifiers that have matching IDs.
expr <- expr[which(rownames(expr) %in% gid),]

# choose which gene identifier to match by:
# If your expr rownames are Ensembl IDs -> use gene_id
# If your expr rownames are gene symbols -> use gene_name
if (all(rownames(expr) %in% gid)) {
  key <- gid
} else if (!is.null(gname) && all(rownames(expr) %in% gname)) {
  key <- gname
} else {
  stop("Row names of expr do not match GTF gene_id or gene_name. ",
       "Use the matching identifier or provide your own gene coord table.")
}
print("Chose ID")

df_anno <- data.frame(
  gene_key = key,
  chr  = as.character(seqnames(gr_gene)),
  start = start(gr_gene),
  end   = end(gr_gene),
  strand= as.character(strand(gr_gene)),
  stringsAsFactors = FALSE
)
print("Data frame")

# keep only genes present in expression
df_anno <- df_anno[!duplicated(df_anno$gene_key), ]
common <- intersect(rownames(expr), df_anno$gene_key)
df_anno <- df_anno[match(common, df_anno$gene_key), ]
expr    <- expr[common, , drop = FALSE]
print("Kept")

# region to use
if (use_tss) {
  # 1-bp TSS interval (BED requires start<end; use [TSS-1, TSS])
  tss <- ifelse(df_anno$strand == "-", df_anno$end, df_anno$start)
  bed_start <- pmax(0L, tss - 1L)  # 0-based
  bed_end   <- tss
} else {
  # full gene body (0-based start, 1-based end)
  bed_start <- pmax(0L, df_anno$start - 1L)
  bed_end   <- df_anno$end
}
print("Region")

# chromosome naming to match genotypes
if (!keep_chr_prefix) {
  df_anno$chr <- sub("^chr", "", df_anno$chr, ignore.case = TRUE)
}
print("Naming")

# phenotype_id = the same IDs you matched by
phenotype_id <- df_anno$gene_key

# ---- 3) reorder expression columns to match genotype sample order ----
psam_path <- paste0(geno_prefix, ".psam")
fam_path  <- paste0(geno_prefix, ".fam")
if (file.exists(psam_path)) {
  psam <- read.table(psam_path, header = TRUE, sep = "\t", check.names = FALSE,
                     stringsAsFactors = FALSE)
  str("Found PSAM")
  str(psam)
  sample_ids <- psam[,1]
} else if (file.exists(fam_path)) {
  fam <- read.table(fam_path, header = FALSE, sep = "", stringsAsFactors = FALSE)
  sample_ids <- fam[[2]]  # IID is column 2
} else {
  stop("Could not find .psam or .fam for 'geno_prefix' = ", geno_prefix)
}
print("Read FAM")

# ensure all samples exist
shared <- intersect(sample_ids, colnames(expr))
str(shared)
expr <- expr[, shared, drop = FALSE]
str(expr)
print("Exists")
# ---- 4) assemble phenotype BED and sort ----
bed <- data.frame(
  `#chr` = df_anno$chr,
  start = as.integer(bed_start),
  end   = as.integer(bed_end),
  phenotype_id = phenotype_id,
  check.names = FALSE
)
bed <- cbind(bed, expr)
str(bed)

# sort by chrom and start (natural human order if possible)
chr_order <- function(x) {
  # try to place 1..22, X, Y, MT at the end if present
  x2 <- toupper(x)
  map <- setNames(seq_along(x2), x2)
  # numeric part
  num <- suppressWarnings(as.integer(x2))
  is_num <- !is.na(num)
  ord <- order(
    ifelse(is_num, 0L, ifelse(x2 %in% c("X","Y","MT","M"), 1L, 2L)),
    ifelse(is_num, num,
           match(x2, c("X","Y","MT","M"), nomatch = 999L)),
    bed$start
  )
  ord
}
bed <- bed[chr_order(bed$`#chr`), ]
print("sorted by chrom")
str(bed)
# ---- 5) write, bgzip, tabix ----
out_tsv <- NULL
write.table(bed, file = out_tsv, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
print("wrote")
# compress & index (requires bgzip/tabix in PATH)
system(sprintf("bgzip -f %s", shQuote(out_tsv)))
system(sprintf("tabix -f -p bed %s.gz", shQuote(out_tsv)))

message("Wrote: ", out_tsv, ".gz and ", out_tsv, ".gz.tbi")

