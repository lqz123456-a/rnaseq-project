# GO/KEGG 富集分析：把显著基因放到功能/通路层面解读
# 输入：results/deseq2_results.txt；输出：results/ 与 figures/

suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
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

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

res <- read.delim(file.path(RESULTS_DIR, "deseq2_results.txt"), row.names = 1)

# 显著基因：padj < 0.05 且 |log2FC| > 1
sig <- rownames(res[!is.na(res$padj) & !is.na(res$log2FoldChange) &
                    res$padj < 0.05 & abs(res$log2FoldChange) > 1, ])
cat("显著基因数：", length(sig), "\n")

# ---- GO 富集（BP：生物过程）----
ego <- enrichGO(gene = sig, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff = 0.05, qvalueCutoff = 0.05)

if (nrow(as.data.frame(ego)) > 0) {
  write.table(as.data.frame(ego), file.path(RESULTS_DIR, "GO_BP_enrichment.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  p1 <- dotplot(ego, showCategory = 15) + ggtitle("GO BP Enrichment")
  ggsave(file.path(FIGURES_DIR, "GO_dotplot.pdf"), p1, width = 8, height = 6)
  ggsave(file.path(FIGURES_DIR, "GO_dotplot.png"), p1, width = 8, height = 6, dpi = 300)
} else {
  cat("GO 无显著富集项\n")
}

# ---- KEGG 通路富集（需联网访问 kegg.jp，网络不通时跳过，不中断脚本）----
kegg_table_path <- file.path(RESULTS_DIR, "KEGG_enrichment.txt")
kegg_pdf_path <- file.path(FIGURES_DIR, "KEGG_dotplot.pdf")
kegg_png_path <- file.path(FIGURES_DIR, "KEGG_dotplot.png")
unlink(c(kegg_table_path, kegg_pdf_path, kegg_png_path))

gene_entrez <- bitr(sig, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
kegg <- tryCatch(
  enrichKEGG(gene = gene_entrez$ENTREZID, organism = "hsa",
             pvalueCutoff = 0.05, qvalueCutoff = 0.05),
  error = function(e) {
    cat("KEGG 连接失败，已跳过：", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
  write.table(as.data.frame(kegg), kegg_table_path,
              sep = "\t", quote = FALSE, row.names = FALSE)
  p2 <- dotplot(kegg, showCategory = 15) + ggtitle("KEGG Pathway Enrichment")
  ggsave(kegg_pdf_path, p2, width = 8, height = 6)
  ggsave(kegg_png_path, p2, width = 8, height = 6, dpi = 300)
  cat("完成：results/KEGG_enrichment.txt；figures/KEGG_dotplot.*\n")
} else {
  cat("KEGG 无显著富集项或服务不可达（GO 已独立保存，网络恢复后可重跑）\n")
}

cat("完成：results/GO_BP_enrichment.txt；figures/GO_dotplot.*\n")
