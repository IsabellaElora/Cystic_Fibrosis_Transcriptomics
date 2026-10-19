# Cystic fibrosis transcriptomics: CF versus non-CF in UNC medium
# Cleaned portfolio workflow based on the original exploratory R script.
# Input: unmodified raw integer gene-count table, genes in rows, samples in columns.
# Run from the repository root: Rscript scripts/01_deseq2_analysis.R

required <- c("DESeq2", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))

input_file <- "data/raw_counts.tsv"
if (!file.exists(input_file)) stop("Missing ", input_file, "; see README for input requirements.")
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

raw <- read.delim(input_file, row.names = 1, check.names = FALSE)
if (anyNA(raw) || anyDuplicated(rownames(raw))) stop("Missing values or duplicated gene IDs.")
if (!all(vapply(raw, is.numeric, logical(1)))) stop("Count columns must be numeric.")
mat <- as.matrix(raw)
if (any(!is.finite(mat)) || any(mat < 0) || any(mat != round(mat))) {
  stop("DESeq2 requires finite, non-negative integer raw counts.")
}
storage.mode(mat) <- "integer"

# Metadata are inferred from the sample names used in the original project.
# Example: 9_NCF_2_UNC. Check against the actual sample sheet before use.
pattern <- "^[0-9]+_(CF|NCF)_([0-9]+)_(UNC|SC)$"
if (!all(grepl(pattern, colnames(mat)))) stop("Unexpected sample names; check metadata parsing.")
parts <- regmatches(colnames(mat), regexec(pattern, colnames(mat)))
meta <- data.frame(
  Sample = colnames(mat),
  Status = factor(vapply(parts, `[`, character(1), 2), levels = c("NCF", "CF")),
  Donor = vapply(parts, `[`, character(1), 3),
  Medium = factor(vapply(parts, `[`, character(1), 4), levels = c("UNC", "SC")),
  row.names = colnames(mat)
)
write.csv(meta, "results/sample_metadata.csv")
unc <- meta$Medium == "UNC"
unc_mat <- mat[, unc, drop = FALSE]
unc_meta <- droplevels(meta[unc, , drop = FALSE])
stopifnot(identical(colnames(unc_mat), rownames(unc_meta)))
if (any(table(unc_meta$Status) < 2)) stop("At least two samples per status are required.")

run_de <- function(count_mat, sample_meta) {
  stopifnot(identical(colnames(count_mat), rownames(sample_meta)))
  dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = sample_meta, design = ~ Status)
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  if (!nrow(dds)) stop("No genes passed the low-count filter.")
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("Status", "CF", "NCF"), alpha = 0.05)
  list(dds = dds, res = res)
}

primary <- run_de(unc_mat, unc_meta)
result_df <- as.data.frame(primary$res)
result_df$GeneID <- rownames(result_df)
write.csv(result_df, "results/DESeq2_UNC_CF_vs_NCF_results.csv", row.names = FALSE)
write.csv(subset(result_df, !is.na(padj) & padj < 0.05),
          "results/significant_genes.csv", row.names = FALSE)

# MA plot: plotMA uses DESeq2's result object directly.
png("figures/MA_plot_UNC_CF_vs_NCF.png", width = 2400, height = 1800, res = 300)
plotMA(primary$res, main = "CF versus non-CF (UNC medium)", ylim = c(-5, 5), alpha = 0.05)
dev.off()

# PCA on variance-stabilised counts (exploratory visualisation).
vsd <- vst(primary$dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "Status", returnData = TRUE)
pct <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, colour = Status)) + geom_point(size = 3) +
  labs(title = "PCA: CF versus non-CF (UNC medium)",
       x = paste0("PC1 (", pct[1], "% variance)"),
       y = paste0("PC2 (", pct[2], "% variance)")) + theme_minimal()
ggsave("figures/PCA_UNC_CF_vs_NCF.png", p, width = 8, height = 6, dpi = 300)

# Volcano plot. Non-significant genes are grey; significant genes are highlighted.
vol <- result_df[!is.na(result_df$padj) & !is.na(result_df$log2FoldChange), ]
vol$significant <- vol$padj < 0.05
vol$neg_log10_padj <- -log10(pmax(vol$padj, .Machine$double.xmin))
v <- ggplot(vol, aes(log2FoldChange, neg_log10_padj, colour = significant)) +
  geom_point(alpha = 0.55, size = 0.7) +
  scale_colour_manual(values = c("FALSE" = "grey65", "TRUE" = "steelblue")) +
  labs(title = "Differential expression: CF versus non-CF", x = "log2 fold change",
       y = "-log10 adjusted p-value", colour = "FDR < 0.05") + theme_minimal()
ggsave("figures/Volcano_UNC_CF_vs_NCF.png", v, width = 8, height = 6, dpi = 300)

# Sensitivity analysis: exclude the originally investigated sample for comparison.
# This is NOT evidence that the sample should be removed from the primary analysis.
sample_to_check <- "9_NCF_2_UNC"
if (sample_to_check %in% colnames(unc_mat)) {
  keep_sample <- colnames(unc_mat) != sample_to_check
  smat <- unc_mat[, keep_sample, drop = FALSE]
  smeta <- droplevels(unc_meta[keep_sample, , drop = FALSE])
  sensitivity <- run_de(smat, smeta)
  shared <- intersect(rownames(primary$res), rownames(sensitivity$res))
  comparison <- data.frame(
    GeneID = shared,
    log2FC_primary = primary$res[shared, "log2FoldChange"],
    log2FC_sensitivity = sensitivity$res[shared, "log2FoldChange"]
  )
  write.csv(comparison, "results/sensitivity_fold_changes.csv", row.names = FALSE)
  diag <- ggplot(comparison, aes(log2FC_primary, log2FC_sensitivity)) +
    geom_point(alpha = 0.25, size = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = "Fold-change comparison between analyses",
         subtitle = "UNC: all samples versus excluding one non-CF sample",
         x = "Log2 fold change: primary analysis",
         y = "Log2 fold change: sensitivity analysis") + theme_minimal()
  ggsave("figures/Sensitivity_fold_change_comparison.png", diag,
         width = 9, height = 7, dpi = 300)
  sig1 <- rownames(primary$res)[!is.na(primary$res$padj) & primary$res$padj < 0.05]
  sig2 <- rownames(sensitivity$res)[!is.na(sensitivity$res$padj) & sensitivity$res$padj < 0.05]
  union_n <- length(union(sig1, sig2))
  metrics <- data.frame(
    primary_significant = length(sig1), sensitivity_significant = length(sig2),
    jaccard = if (union_n) length(intersect(sig1, sig2)) / union_n else NA_real_,
    spearman_log2FC = cor(comparison$log2FC_primary,
                          comparison$log2FC_sensitivity, method = "spearman", use = "complete.obs")
  )
  write.csv(metrics, "results/sensitivity_summary.csv", row.names = FALSE)
} else {
  message("Sensitivity analysis skipped: sample not present: ", sample_to_check)
}
message("Completed DESeq2 workflow. Inspect results/ and figures/.")
