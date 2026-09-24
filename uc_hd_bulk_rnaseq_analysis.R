## UC vs HD bulk RNA-seq analysis
## A group = UC; B group = HD
## Input file: E:/肠上皮-乳酸化ACSS2/转录组测序/20250503-1.csv

## =========================
## 1. Packages
## =========================
cran_pkgs <- c("data.table", "ggplot2", "pheatmap", "RColorBrewer")
bioc_pkgs <- c("edgeR", "limma")

install_if_missing_cran <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) install.packages(miss, repos = "https://cloud.r-project.org")
}

install_if_missing_bioc <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(miss, ask = FALSE, update = FALSE)
  }
}

install_if_missing_cran(cran_pkgs)
install_if_missing_bioc(bioc_pkgs)

library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

## =========================
## 2. Parameters
## =========================
input_file <- "E:/肠上皮-乳酸化ACSS2/转录组测序/20250503-1.csv"
outdir <- "E:/肠上皮-乳酸化ACSS2/转录组测序/UC_HD_R_analysis"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

gene_col <- "gene_name"
min_cpm <- 1
min_samples <- 2
adj_p_cutoff <- 0.05
logfc_cutoff <- 1

## Genes shown in your example heatmap: glycolysis + ferroptosis related genes
target_genes <- c(
  "ALOX5AP", "PKM", "TFRC", "ACSL4", "HK2",
  "LOX", "HK1", "PFKM", "HK3", "ALDOC",
  "PFKL", "GAPDH", "GPX4", "NFE2L2", "PFKP"
)

group_colors <- c(UC = "#E41A1C", HD = "#4DA3D9")

## =========================
## 3. Read expression matrix
## =========================
raw_dt <- fread(input_file, data.table = FALSE, check.names = FALSE)

if (!gene_col %in% colnames(raw_dt)) {
  stop("Cannot find gene column: ", gene_col)
}

sample_cols <- setdiff(colnames(raw_dt), gene_col)

metadata <- data.frame(
  sample_id = sample_cols,
  group = ifelse(grepl("^A", sample_cols, ignore.case = TRUE), "UC",
                 ifelse(grepl("^B", sample_cols, ignore.case = TRUE), "HD", NA)),
  stringsAsFactors = FALSE
)

if (any(is.na(metadata$group))) {
  stop("Some samples cannot be assigned to UC/HD by A*/B* naming: ",
       paste(metadata$sample_id[is.na(metadata$group)], collapse = ", "))
}

metadata$group <- factor(metadata$group, levels = c("HD", "UC"))
metadata <- metadata[order(metadata$group, metadata$sample_id), ]
metadata$plot_sample <- ave(as.character(metadata$group), metadata$group,
                            FUN = function(x) paste0(x, seq_along(x)))

expr <- raw_dt[, sample_cols, drop = FALSE]
expr[] <- lapply(expr, function(x) as.numeric(as.character(x)))
gene_names <- raw_dt[[gene_col]]

keep_gene <- !is.na(gene_names) & gene_names != ""
expr <- expr[keep_gene, , drop = FALSE]
gene_names <- gene_names[keep_gene]

## If duplicate gene symbols exist, sum counts by gene.
expr_mat <- rowsum(as.matrix(expr), group = gene_names, reorder = FALSE)
expr_mat <- expr_mat[, metadata$sample_id, drop = FALSE]

write.csv(metadata, file.path(outdir, "sample_metadata.csv"), row.names = FALSE)
write.csv(expr_mat, file.path(outdir, "raw_count_matrix_by_gene.csv"))

## =========================
## 4. edgeR filtering + TMM + voom
## =========================
dge <- DGEList(counts = expr_mat, group = metadata$group)
keep <- filterByExpr(dge, group = metadata$group, min.count = 10)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

design <- model.matrix(~ 0 + group, data = metadata)
colnames(design) <- gsub("^group", "", colnames(design))

pdf(file.path(outdir, "voom_mean_variance.pdf"), width = 6, height = 5)
v <- voom(dge, design, plot = TRUE)
dev.off()

logcpm <- v$E
colnames(logcpm) <- metadata$plot_sample
write.csv(logcpm, file.path(outdir, "TMM_voom_logCPM.csv"))

## =========================
## 5. PCA plot
## =========================
var_order <- order(apply(logcpm, 1, var), decreasing = TRUE)
pca_mat <- logcpm[var_order[seq_len(min(5000, length(var_order)))], , drop = FALSE]
pca <- prcomp(t(pca_mat), scale. = FALSE)
pca_df <- data.frame(
  sample = metadata$plot_sample,
  original_sample = metadata$sample_id,
  group = metadata$group,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)

pc_var <- round(100 * summary(pca)$importance[2, 1:2], 2)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group, fill = group)) +
  stat_ellipse(geom = "polygon", type = "norm", alpha = 0.15, linewidth = 0.7) +
  stat_ellipse(type = "norm", linewidth = 0.8) +
  geom_point(size = 3.2) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(x = paste0("PC1 (", pc_var[1], "%)"),
       y = paste0("PC2 (", pc_var[2], "%)"),
       color = NULL, fill = NULL) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = c(0.92, 0.15),
    legend.background = element_rect(color = "black", fill = "white"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black", face = "bold")
  )

ggsave(file.path(outdir, "PCA_UC_vs_HD.pdf"), p_pca, width = 6, height = 5)
ggsave(file.path(outdir, "PCA_UC_vs_HD.png"), p_pca, width = 6, height = 5, dpi = 300)

## =========================
## 6. Differential expression: UC vs HD
## Positive logFC means higher in UC.
## =========================
fit <- lmFit(v, design)
contrast_mat <- makeContrasts(UC_vs_HD = UC - HD, levels = design)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

deg_all <- topTable(fit2, coef = "UC_vs_HD", number = Inf, sort.by = "P")
deg_all$gene_name <- rownames(deg_all)
deg_all$change <- "NS"
deg_all$change[deg_all$adj.P.Val < adj_p_cutoff & deg_all$logFC >= logfc_cutoff] <- "Up_in_UC"
deg_all$change[deg_all$adj.P.Val < adj_p_cutoff & deg_all$logFC <= -logfc_cutoff] <- "Down_in_UC"

deg_all <- deg_all[, c("gene_name", setdiff(colnames(deg_all), "gene_name"))]
deg_sig <- deg_all[deg_all$change != "NS", ]

write.csv(deg_all, file.path(outdir, "DEG_all_UC_vs_HD.csv"), row.names = FALSE)
write.csv(deg_sig, file.path(outdir, "DEG_sig_adjP0.05_logFC1_UC_vs_HD.csv"), row.names = FALSE)

## =========================
## 7. Volcano plot
## =========================
p_volcano <- ggplot(deg_all, aes(logFC, -log10(adj.P.Val), color = change)) +
  geom_point(alpha = 0.75, size = 1.6) +
  scale_color_manual(values = c(Up_in_UC = "#E41A1C", Down_in_UC = "#377EB8", NS = "grey70")) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(adj_p_cutoff), linetype = "dashed", color = "grey40") +
  labs(x = "log2 Fold Change (UC / HD)", y = "-log10 adjusted P value", color = NULL) +
  theme_classic(base_size = 13) +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"))

ggsave(file.path(outdir, "Volcano_UC_vs_HD.pdf"), p_volcano, width = 6, height = 5)
ggsave(file.path(outdir, "Volcano_UC_vs_HD.png"), p_volcano, width = 6, height = 5, dpi = 300)

## =========================
## 8. Target-gene heatmap
## =========================
rownames_upper <- toupper(rownames(logcpm))
target_upper <- toupper(target_genes)
idx <- match(target_upper, rownames_upper)
found_genes <- rownames(logcpm)[idx[!is.na(idx)]]
missing_genes <- target_genes[is.na(idx)]

writeLines(missing_genes, file.path(outdir, "missing_target_genes.txt"))

if (length(found_genes) < 2) {
  stop("Too few target genes found in expression matrix. Check gene symbols.")
}

heat_mat <- logcpm[found_genes, metadata$plot_sample, drop = FALSE]
heat_z <- t(scale(t(heat_mat)))
heat_z[is.na(heat_z)] <- 0
heat_z[heat_z > 2] <- 2
heat_z[heat_z < -2] <- -2

ann_col <- data.frame(Group = metadata$group)
rownames(ann_col) <- metadata$plot_sample
ann_colors <- list(Group = group_colors)

pdf(file.path(outdir, "Heatmap_glycolysis_ferroptosis.pdf"), width = 7, height = 7)
pheatmap(
  heat_z,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  breaks = seq(-2, 2, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  border_color = NA,
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 9,
  main = "Glycolysis & Ferroptosis"
)
dev.off()

png(file.path(outdir, "Heatmap_glycolysis_ferroptosis.png"), width = 2100, height = 2100, res = 300)
pheatmap(
  heat_z,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  breaks = seq(-2, 2, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  border_color = NA,
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 9,
  main = "Glycolysis & Ferroptosis"
)
dev.off()

write.csv(heat_mat, file.path(outdir, "target_gene_logCPM.csv"))
write.csv(heat_z, file.path(outdir, "target_gene_zscore_for_heatmap.csv"))

## =========================
## 9. Brief summary
## =========================
summary_txt <- c(
  paste0("Input file: ", input_file),
  paste0("Output directory: ", outdir),
  paste0("Samples: ", paste(metadata$sample_id, metadata$group, sep = "=", collapse = "; ")),
  paste0("Genes before filtering: ", nrow(expr_mat)),
  paste0("Genes after filtering: ", nrow(logcpm)),
  paste0("Significant DEGs: ", nrow(deg_sig)),
  paste0("Up in UC: ", sum(deg_all$change == "Up_in_UC")),
  paste0("Down in UC: ", sum(deg_all$change == "Down_in_UC")),
  paste0("Missing target genes: ", ifelse(length(missing_genes) == 0, "None", paste(missing_genes, collapse = ", ")))
)
writeLines(summary_txt, file.path(outdir, "analysis_summary.txt"))

sink(file.path(outdir, "sessionInfo.txt"))
sessionInfo()
sink()

message("Analysis finished. Results saved to: ", outdir)
