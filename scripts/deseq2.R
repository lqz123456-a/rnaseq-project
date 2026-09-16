# 真实 RNA-seq 差异分析（airway 8 样本，salmon + tximport + DESeq2）
# 原理：salmon 把 reads 拟比对/定量到转录本；tximport 按 tx2gene 汇总到基因级；
#       DESeq2 用负二项分布建模，对 treated vs untreated 做 Wald 检验 + BH 校正（padj）。
# 画图统一用 ggplot2。

suppressMessages(library(tximport))
suppressMessages(library(DESeq2))
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
QUANT_DIR <- file.path(PROJECT_ROOT, "quant")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

# 1) 转录本 -> 基因 映射（由 make_tx2gene.py 生成）
tx2gene <- read.delim(file.path(RESULTS_DIR, "tx2gene.tsv"),
                      header = FALSE, col.names = c("tx", "gene"))

# 2) 8 个样本的 salmon 定量文件
samples <- c("SRR1039508", "SRR1039509", "SRR1039512", "SRR1039513",
             "SRR1039516", "SRR1039517", "SRR1039520", "SRR1039521")
files <- file.path(QUANT_DIR, samples, "quant.sf")
names(files) <- samples

# 3) 导入并汇总到基因级
#    dropInfReps=TRUE   ：跳过 salmon 推断重复信息
#    ignoreTxVersion=TRUE：忽略转录本 ID 的版本后缀
txi <- tximport(files, type = "salmon", tx2gene = tx2gene,
                dropInfReps = TRUE, ignoreTxVersion = TRUE)

# 4) 读入样本信息表
coldata <- read.table(file.path(RESULTS_DIR, "coldata.txt"),
                      header = TRUE, row.names = 1)
coldata <- coldata[colnames(txi$counts), , drop = FALSE]

# 5) 构建 DESeqDataSet 并跑完整差异分析
#    design = ~cell + condition：配对模型，扣除 4 个细胞系的基线差异后再估计处理效应
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~cell + condition)
dds <- DESeq(dds)

# 6) 提取结果：treated vs untreated，按 padj 排序
res <- results(dds, contrast = c("condition", "treated", "untreated"))
res <- res[order(res$padj), ]

# 7) 输出结果表 + 摘要
write.table(as.data.frame(res), file.path(RESULTS_DIR, "deseq2_results.txt"),
            sep = "\t", quote = FALSE)
summary(res)

# 8) 转成数据框，并加上基因名
res <- as.data.frame(res)
nm <- read.delim(file.path(RESULTS_DIR, "id2name.tsv"),
                 header = FALSE, col.names = c("id", "name"))
res$name <- nm$name[match(rownames(res), nm$id)]

# 9) 分组着色：上调 = padj<0.05 且 log2FC>1；下调 = padj<0.05 且 log2FC<-1
res$dir <- ifelse(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1, "up",
          ifelse(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1, "down", "ns"))
res$dir <- factor(res$dir, levels = c("up", "down", "ns"))
cols <- c(up = "#C0392B", down = "#2471A3", ns = "grey75")

# 10) MA 图
p1 <- ggplot(res, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = dir), size = 0.7, alpha = 0.55) +
  scale_x_log10() +
  scale_color_manual(values = cols, labels = c("显著上调", "显著下调", "不显著")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  coord_cartesian(ylim = c(-4, 4)) +
  labs(x = "mean of normalized counts (log10)",
       y = "log2 fold change",
       title = "MA Plot - treated vs untreated",
       color = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

ggsave(file.path(FIGURES_DIR, "MAplot.pdf"), p1, width = 7.5, height = 5.5)
ggsave(file.path(FIGURES_DIR, "MAplot.png"), p1, width = 7.5, height = 5.5, dpi = 300)

# 11) 火山图，标注 top 10 基因
top <- head(res[order(res$padj, na.last = TRUE), ], 10)

p2 <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = dir), size = 0.7, alpha = 0.55) +
  scale_color_manual(values = cols, labels = c("显著上调", "显著下调", "不显著")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  ggrepel::geom_text_repel(data = top, aes(label = name), size = 3.2,
                           max.overlaps = 20, seed = 42, color = "grey20") +
  coord_cartesian(xlim = c(-8, 8)) +
  labs(x = "log2 fold change", y = "-log10(adjusted p-value)",
       title = "Volcano Plot - treated vs untreated",
       color = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

ggsave(file.path(FIGURES_DIR, "volcano.pdf"), p2, width = 8, height = 6)
ggsave(file.path(FIGURES_DIR, "volcano.png"), p2, width = 8, height = 6, dpi = 300)

cat("完成：results/deseq2_results.txt；figures/MAplot.*、figures/volcano.*\n")