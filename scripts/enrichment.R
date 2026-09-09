# GO/KEGG 富集分析：把显著基因放到功能/通路层面解读
# 输入：第 8 步的 deseq2_results.txt；输出：富集结果表 + 气泡图
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(ggplot2))

res <- read.delim("deseq2_results.txt", row.names=1)

# 显著基因（padj<0.05；想更严格可加 & abs(log2FoldChange)>1）
sig <- rownames(res[!is.na(res$padj) & res$padj < 0.05, ])
cat("显著基因数：", length(sig), "\n")

# ---- GO 富集（BP：生物过程）----
ego <- enrichGO(gene=sig, OrgDb=org.Hs.eg.db, keyType="ENSEMBL",
                ont="BP", pAdjustMethod="BH",
                pvalueCutoff=0.05, qvalueCutoff=0.05)

if (nrow(as.data.frame(ego)) > 0) {
  write.table(as.data.frame(ego), "GO_BP_enrichment.txt",
              sep="\t", quote=FALSE, row.names=FALSE)
  p1 <- dotplot(ego, showCategory=15) + ggtitle("GO BP Enrichment")
  ggsave("GO_dotplot.pdf", p1, width=8, height=6)
  ggsave("GO_dotplot.png", p1, width=8, height=6, dpi=300)
} else { cat("GO 无显著富集项\n") }

# ---- KEGG 通路富集（需联网访问 kegg.jp，你网络已实测可达）----
gene_entrez <- bitr(sig, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
kegg <- enrichKEGG(gene=gene_entrez$ENTREZID, organism="hsa",
                   pvalueCutoff=0.05, qvalueCutoff=0.05)

if (nrow(as.data.frame(kegg)) > 0) {
  write.table(as.data.frame(kegg), "KEGG_enrichment.txt",
              sep="\t", quote=FALSE, row.names=FALSE)
  p2 <- dotplot(kegg, showCategory=15) + ggtitle("KEGG Pathway Enrichment")
  ggsave("KEGG_dotplot.pdf", p2, width=8, height=6)
  ggsave("KEGG_dotplot.png", p2, width=8, height=6, dpi=300)
} else { cat("KEGG 无显著富集项（网络不通可跳过，GO 已够用）\n") }

cat("完成！GO/KEGG 结果表 + 气泡图已生成\n")
