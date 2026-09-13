# ggplot2 版 MA 图 + 火山图（读取 results/deseq2_results.txt，不重跑差异分析）
# 用法：Rscript scripts/plot_ggplot2.R

suppressMessages(library(ggplot2))

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    return(normalizePath(".", mustWork = TRUE))
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
}

PROJECT_ROOT <- normalizePath(file.path(dirname(get_script_path()), ".."), mustWork = TRUE)
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
FIGURES_DIR <- file.path(PROJECT_ROOT, "figures")

dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

res <- read.delim(file.path(RESULTS_DIR, "deseq2_results.txt"), row.names = 1)

# 基因名映射，用于火山图标注
if (file.exists(file.path(RESULTS_DIR, "id2name.tsv"))) {
  nm <- read.delim(file.path(RESULTS_DIR, "id2name.tsv"),
                   header = FALSE, col.names = c("id", "name"))
  res$name <- nm$name[match(rownames(res), nm$id)]
} else {
  res$name <- rownames(res)
}

# 分组着色：上调 = padj<0.05 且 log2FC>1；下调 = padj<0.05 且 log2FC<-1
res$dir <- ifelse(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1, "up",
          ifelse(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1, "down", "ns"))
res$dir <- factor(res$dir, levels = c("up", "down", "ns"))
col_manual <- c(up = "#C0392B", down = "#2471A3", ns = "grey75")

# ---------- MA 图 ----------
p1 <- ggplot(res, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = dir), size = 0.7, alpha = 0.55) +
  scale_x_log10() +
  scale_color_manual(values = col_manual, labels = c("显著上调", "显著下调", "不显著")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  coord_cartesian(ylim = c(-4, 4)) +
  labs(x = "mean of normalized counts (log10)", y = "log2 fold change",
       title = "MA Plot - Dexamethasone vs Control (airway)",
       subtitle = "red: padj<0.05 & log2FC>1; blue: padj<0.05 & log2FC<-1",
       color = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIGURES_DIR, "MAplot.pdf"), p1, width = 7.5, height = 5.5)
ggsave(file.path(FIGURES_DIR, "MAplot.png"), p1, width = 7.5, height = 5.5, dpi = 300)

# ---------- 火山图 ----------
top <- head(res[order(res$padj, na.last = TRUE), ], 10)

p2 <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = dir), size = 0.7, alpha = 0.55) +
  scale_color_manual(values = col_manual, labels = c("显著上调", "显著下调", "不显著")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  ggrepel::geom_text_repel(data = top, aes(label = name), size = 3.2,
                           max.overlaps = 20, seed = 42, color = "grey20") +
  coord_cartesian(xlim = c(-8, 8)) +
  labs(x = "log2 fold change", y = "-log10(adjusted p-value)",
       title = "Volcano Plot - Dexamethasone vs Control (airway)",
       subtitle = "dashed: |log2FC|=1, padj=0.05",
       color = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIGURES_DIR, "volcano.pdf"), p2, width = 8, height = 6)
ggsave(file.path(FIGURES_DIR, "volcano.png"), p2, width = 8, height = 6, dpi = 300)

cat("完成：figures/MAplot.*、figures/volcano.*\n")