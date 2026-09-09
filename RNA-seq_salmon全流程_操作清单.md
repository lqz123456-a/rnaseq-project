# 真实 RNA-seq 完整跑通・salmon 全流程版（国内网络专用）

> **数据集：GSE52778「airway」** —— 人气道平滑肌细胞，地塞米松（dexamethasone）处理，DESeq2 官方教程使用的经典真实数据。8 个样本（4 对照 + 4 处理），双端测序。全部链接已于 2026-09-05 在本机网络实测通过。
>
> **本版特点**：用 salmon（转录本拟比对 + 定量）替代 STAR，内存需求～4–8 GB，适配本机 WSL 7.6 GB 的环境，不会 OOM。产出基因级差异表达结果（MA / 火山图）与 GO/KEGG 富集分析，不产出基因组 BAM。
> 与 STAR 版的区别：① 不下载基因组 fasta（salmon 用不到）② 第 3 步变 salmon index ③ 第 6 步变 salmon quant ④ 第 8 步用 tximport 导入 DESeq2。



***

## 0. 为什么之前下载总失败（已在本机实测）



| 数据源                 | 性质         | 实测结果           | 结论           |
| ------------------- | ---------- | -------------- | ------------ |
| NCBI FTP（美国）        | SRA / 参考序列 | ❌ 连接超时         | **别用**       |
| Ensembl 主站（英国）      | 参考基因组      | ❌ 连接超时         | **别用**       |
| UCSC（美国）            | 参考基因组      | ❌ 连接超时         | **别用**       |
| **EBI ENA（欧洲）**     | SRA 原始数据   | ✅ 200，0.8–1.7s | **用这个下原始数据** |
| **EBI Ensembl（欧洲）** | 转录组 / GTF  | ✅ 200，0.9s     | **用这个下转录组**  |
| **NGDC 国家基因组中心**    | 国内数据 / 镜像  | ✅ 200，0.3s     | 备用           |
| **CNGB 国家基因库**      | 国内数据 / 镜像  | ✅ 200，0.2s     | 备用           |
| 清华 / 中科大 conda 镜像   | 软件         | ✅ 200，0.5–1.7s | 装软件用清华       |

**核心结论：本机网络到 "美国源"（NCBI/UCSC/Ensembl 主站）基本不通，这是之前所有下载失败的根源；解决办法是全部改用 "欧洲源 EBI + 国内源 NGDC/CNGB/ 清华镜像"。**



***

## 1. 环境准备（conda 配置清华镜像）



```
conda config --set show_channel_urls yes
```

编辑 `~/.condarc`，写入：



```
channels:

  - defaults

default_channels:

  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main

  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free

  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r

custom_channels:

  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud

  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
```

创建分析环境并安装工具（salmon 替代 star/subread）：



```
conda create -n rnaseq -y fastqc multiqc fastp salmon aria2

conda activate rnaseq

conda install -y r-base bioconductor-deseq2 bioconductor-tximport r-ggplot2 r-ggrepel bioconductor-clusterprofiler bioconductor-org.hs.eg.db
```

> 若 `conda install` 这条事务卡死不动（本机实测三次挂起），拆开分步装：先装不含 clusterProfiler/org.Hs.eg.db 的依赖，再按第 9 节的源码包方式补这两个数据包。
> 可选优化：若想让本机 WSL 内存更充裕（salmon 其实不需要，但能加快建索引），可在 Windows 建 `C:\Users\<用户名>\.wslconfig` 写入 `[wsl2]` + `memory=12GB` + `swap=8GB`，然后 `wsl --shutdown` 重启。



***

## 2. 下载转录组 + GTF（EBI Ensembl release-116，已实测可达）

salmon 只需要 **转录组序列（cdna）** 和 **GTF 注释**，不需要下载基因组 fasta。



```
mkdir -p ~/rnaseq/ref && cd ~/rnaseq/ref

# 转录组序列（约 500 MB）

wget -c https://ftp.ebi.ac.uk/ensemblorg/pub/release-116/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# 基因注释 GTF（约 60 MB，用于提取 tx2gene 映射）

wget -c https://ftp.ebi.ac.uk/ensemblorg/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz

gzip -dk Homo_sapiens.GRCh38.116.gtf.gz
```

> 下载太慢 / 中断：`wget -c` 断点续传；更快用 `aria2c -x 16 -s 16 -c <URL>`。
> 若之前已下载过基因组 fasta（GRCh38.dna.primary_assembly.fa），**保留即可** —— 本版用不到，以后跑 STAR / 变异分析还用得上，不浪费。



***

## 3. 建 salmon 索引（内存～4–8 GB，本机 7.6 GB 的 WSL 可跑）



```
cd ~/rnaseq

salmon index -t ref/Homo_sapiens.GRCh38.cdna.all.fa.gz -i salmon_index
```

> 如果内存仍然紧张，先关掉其他程序再跑。这一步一般 5–15 分钟。



***

## 4. 下载原始测序数据（airway 数据集，从 ENA 而非 NCBI！）

**样本分组（4 对照 + 4 处理）：**



| 分组                | 样本（SRR）                                        |
| ----------------- | ---------------------------------------------- |
| 对照组 untreated     | SRR1039508, SRR1039512, SRR1039516, SRR1039520 |
| 处理组 dexamethasone | SRR1039509, SRR1039513, SRR1039517, SRR1039521 |

**一次下载全部 8 个样本（4 对照 + 4 处理；DESeq2 每组最少 3 个，8 个更稳）：**



```
mkdir -p ~/rnaseq/fastq && cd ~/rnaseq/fastq

for u in "SRR1039508 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508" "SRR1039509 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/009/SRR1039509/SRR1039509" "SRR1039512 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/002/SRR1039512/SRR1039512" "SRR1039513 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/003/SRR1039513/SRR1039513" "SRR1039516 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/006/SRR1039516/SRR1039516" "SRR1039517 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/007/SRR1039517/SRR1039517" "SRR1039520 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039520/SRR1039520" "SRR1039521 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/001/SRR1039521/SRR1039521" ; do

  set -- $u; s=$1; base=$2

  aria2c -x 8 -s 8 -c "${base}_1.fastq.gz"

  aria2c -x 8 -s 8 -c "${base}_2.fastq.gz"

done

# 下载完成后先校验完整性（gzip 尾部块丢失会静默导致 fastqc 失败，见第 10 节排障）

for f in *.fastq.gz; do

  gzip -t "$f" && echo "$f OK" || echo "$f 损坏，删除后重新下载"

done
```

**换数据集时查链接（通用方法）**：



```
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR1039508&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes"

# 输出里的 ftp.sra.ebi.ac.uk/... 前面加 https:// 即可下载
```



***

## 5. 质控（fastqc + multiqc）



```
cd ~/rnaseq

mkdir -p qc        # fastqc 不会自动创建输出目录，必须先建

fastqc -t 8 fastq/*.fastq.gz -o qc/

multiqc qc/ -o qc/

# 打开 qc/multiqc_report.html 查看质量汇总
```

> 质量好可不修剪；需要时：`fastp -i in_1.fq.gz -I in_2.fq.gz -o out_1.fq.gz -O out_2.fq.gz`



***

## 6. salmon 定量（转录本拟比对，替代 STAR 比对）



```
cd ~/rnaseq

for s in SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516 SRR1039517 SRR1039520 SRR1039521; do

salmon quant -i salmon_index -l A -1 fastq/${s}_1.fastq.gz -2 fastq/${s}_2.fastq.gz -p 8 -o quant/${s}

done

# 每个样本会生成 quant/SRRxxxx/quant.sf（转录本定量结果）
```



***

## 7. 提取 tx2gene + id2name 映射（salmon 是转录本级，汇总到基因级要用）

一次脚本生成两个映射文件：



* `tx2gene.tsv`：转录本 → 基因（tximport 汇总用）

* `id2name.tsv`：基因 ID → 基因名（ggplot2 火山图标注 top 基因用）



```
cd ~/rnaseq

python3 - <<'EOF'

import re

tx2g = {}

id2name = {}

for line in open("ref/Homo_sapiens.GRCh38.116.gtf"):

    if line.startswith("#") or len(line.split("\t")) < 9:
        continue    # 跳过 #! 注释行（Ensembl GTF 开头有），否则第 9 列取不到会报错

    attrs = line.split("\t")[8]

    if "\ttranscript\t" in line:

        g = re.search(r'gene_id "([^"]+)"', attrs)

        t = re.search(r'transcript_id "([^"]+)"', attrs)

        if g and t:

            tx2g[t.group(1)] = g.group(1)

    elif "\tgene\t" in line:

        g = re.search(r'gene_id "([^"]+)"', attrs)

        n = re.search(r'gene_name "([^"]+)"', attrs)

        if g and n:

            id2name[g.group(1)] = n.group(1)

with open("tx2gene.tsv", "w") as f:

    for t, g in tx2g.items():

        f.write(f"{t}\t{g}\n")

with open("id2name.tsv", "w") as f:

    for gid, name in id2name.items():

        f.write(f"{gid}\t{name}\n")

print(f"tx2gene.tsv {len(tx2g)} 条；id2name.tsv {len(id2name)} 条")

EOF
```



***

## 8. 差异表达（tximport + DESeq2，真实数据 + 注释版 R 脚本）

**样本信息表 `coldata.txt`：**



```
sample	condition

SRR1039508	untreated

SRR1039509	treated

SRR1039512	untreated

SRR1039513	treated

SRR1039516	untreated

SRR1039517	treated

SRR1039520	untreated

SRR1039521	treated
```

**R 脚本 `deseq2.R`（每步都注释，能看懂）：**



```
# 原理一句话：salmon 先把 reads 拟比对/定量到转录本；

# tximport 按 tx2gene 把转录本汇总到基因；

# DESeq2 用负二项分布建模，对处理组 vs 对照组做 Wald 检验 + BH 校正（padj）；

# 画图统一用 ggplot2（不用 base R）。

suppressMessages(library(tximport))

suppressMessages(library(DESeq2))

suppressMessages(library(ggplot2))

# 1) 转录本 -> 基因 映射（第 7 步生成）

tx2gene <- read.delim("tx2gene.tsv", header=FALSE, col.names=c("tx","gene"))

# 2) 指定 8 个样本的 salmon 定量文件（顺序必须与 coldata.txt 一致）

samples <- c("SRR1039508","SRR1039509","SRR1039512","SRR1039513",

             "SRR1039516","SRR1039517","SRR1039520","SRR1039521")

files <- file.path("quant", samples, "quant.sf")

names(files) <- samples

# 3) 导入并汇总到基因级

#    dropInfReps=TRUE   ：跳过 salmon 推断重复信息（读它需要 jsonlite，常规分析用不到）

#    ignoreTxVersion=TRUE：忽略转录本 ID 的版本后缀（salmon 的 ENST...N vs tx2gene 的 ENST...）

txi <- tximport(files, type="salmon", tx2gene=tx2gene,

                dropInfReps=TRUE, ignoreTxVersion=TRUE)

# 4) 读入样本信息（与 samples 顺序一致）

coldata <- read.table("coldata.txt", header=TRUE, row.names=1)

coldata <- coldata[colnames(txi$counts), , drop=FALSE]

# 5) 构建 DESeqDataSet 并跑差异分析

dds <- DESeqDataSetFromTximport(txi, colData=coldata, design=~condition)

dds <- DESeq(dds)

# 6) 提取结果：treated vs untreated，按 padj 排序

res <- results(dds, contrast=c("condition","treated","untreated"))

res <- res[order(res$padj), ]

# 7) 输出结果表 + 摘要（看有多少显著基因）

write.table(as.data.frame(res), "deseq2_results.txt", sep="\t", quote=FALSE)

summary(res)

# 8) 转成数据框，并加上基因名（第 7 步生成的 id2name.tsv，用于火山图标注）

res <- as.data.frame(res)

nm <- read.delim("id2name.tsv", header=FALSE, col.names=c("id","name"))

res$name <- nm$name[match(rownames(res), nm$id)]

# 9) 分组着色：上调 = padj<0.05 且 log2FC>1；下调 = padj<0.05 且 log2FC<-1

res$dir <- ifelse(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1, "up",

          ifelse(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1, "down", "ns"))

res$dir <- factor(res$dir, levels=c("up","down","ns"))

cols <- c(up="#C0392B", down="#2471A3", ns="grey75")

# 10) MA 图（ggplot2）

p1 <- ggplot(res, aes(x=baseMean, y=log2FoldChange)) +

  geom_point(aes(color=dir), size=0.7, alpha=0.55) +

  scale_x_log10() +

  scale_color_manual(values=cols, labels=c("显著上调","显著下调","不显著")) +

  geom_hline(yintercept=0, linetype="dashed", color="grey40") +

  coord_cartesian(ylim=c(-4, 4)) +

  labs(x="mean of normalized counts (log10)", y="log2 fold change",

       title="MA Plot — treated vs untreated", color=NULL) +

  theme_bw(base_size=13) +

  theme(legend.position="top", plot.title=element_text(face="bold"))

ggsave("MAplot_ggplot2.pdf", p1, width=7.5, height=5.5)

ggsave("MAplot_ggplot2.png", p1, width=7.5, height=5.5, dpi=300)

# 11) 火山图（ggplot2，标注 top 10 基因）

top <- head(res[order(res$padj, na.last=TRUE), ], 10)

p2 <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj))) +

  geom_point(aes(color=dir), size=0.7, alpha=0.55) +

  scale_color_manual(values=cols, labels=c("显著上调","显著下调","不显著")) +

  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="grey50") +

  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey50") +

  ggrepel::geom_text_repel(data=top, aes(label=name), size=3.2,

                           max.overlaps=20, seed=42, color="grey20") +

  coord_cartesian(xlim=c(-8, 8)) +

  labs(x="log2 fold change", y="-log10(adjusted p-value)",

       title="Volcano Plot — treated vs untreated", color=NULL) +

  theme_bw(base_size=13) +

  theme(legend.position="top", plot.title=element_text(face="bold"))

ggsave("volcano_ggplot2.pdf", p2, width=8, height=6)

ggsave("volcano_ggplot2.png", p2, width=8, height=6, dpi=300)

cat("完成！结果：deseq2_results.txt；图：MAplot_ggplot2.pdf/png、volcano_ggplot2.pdf/png\n")
```



```
Rscript deseq2.R
```



***

## 9. GO/KEGG 富集分析（最后一步，做完即收尾）

**富集分析的原理**：差异分析只回答 "哪些基因变了"，富集分析把显著基因放到 "功能 / 通路" 层面 —— 比如 "这批基因富集在免疫应答相关通路"，让结果有生物学意义，是面试讲 "所以呢" 的关键一环。

**R 脚本 `enrichment.R`（clusterProfiler，已在本机实测可跑）：**



```
# GO/KEGG 富集分析：把显著基因放到功能/通路层面解读

# 输入：第 8 步的 deseq2_results.txt；输出：富集结果表 + 气泡图

suppressMessages(library(clusterProfiler))

suppressMessages(library(org.Hs.eg.db))

suppressMessages(library(ggplot2))   # dotplot 加标题需要

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

# ---- KEGG 通路富集（需联网访问 kegg.jp，本机网络已实测可达）----

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
```



```
Rscript enrichment.R
```

**产出**：`GO_BP_enrichment.txt`、`KEGG_enrichment.txt`（结果表）+ `GO_dotplot.pdf/png`、`KEGG_dotplot.pdf/png`（气泡图）。

**本机实测结果（真实运行数字，可直接用于简历 / 面试）：**



| 项目                                                   | 结果                 |
| ---------------------------------------------------- | ------------------ |
| 显著基因数（padj<0.05）                                     | 2140               |
| GO BP 显著条目                                           | 741 个              |
| KEGG 显著通路                                            | 109 条              |
| GO Top：cellular response to peptide hormone stimulus | padj=3.6e-08，76 基因 |
| GO Top：regulation of actin filament-based process    | padj=3.6e-08，86 基因 |
| GO Top：response to hypoxia（低氧应答）                     | padj=1.6e-06，68 基因 |
| KEGG Top：Focal adhesion（黏着斑）                         | padj=1.1e-08，56 基因 |
| KEGG Top：PI3K-Akt signaling pathway                  | padj=3.4e-06，75 基因 |
| KEGG Top：Regulation of actin cytoskeleton            | padj=3.4e-06，54 基因 |

这些条目与 "地塞米松（糖皮质激素）处理" 的已知生物学一致（激素应答、细胞骨架重塑、缺氧应答），面试时能讲出 "结果符合预期，说明分析流程可靠"。

> **如果 `conda install bioconductor-clusterprofiler` 卡死**（本机实测：三次都在事务执行阶段无限挂起，与网络无关）：放弃 conda，直接从 Bioconductor 官方源装源码包（纯数据包无需编译）：

```
# 下载完整源码包（Galaxy Depot 镜像，国内可达；下载后 md5sum 核对）

cd ~/rnaseq

wget -c https://depot.galaxyproject.org/software/bioconductor-org.hs.eg.db/bioconductor-org.hs.eg.db_3.22.0_src_all.tar.gz

wget -c https://depot.galaxyproject.org/software/bioconductor-go.db/bioconductor-go.db_3.22.0_src_all.tar.gz

# 官方 md5：org.Hs.eg.db = e80cac6ec018a95aea4f7530350e80a2，GO.db = 5ae5557afa56227c4c9c145907b1f585

# 安装进当前环境（R CMD INSTALL 会自己找到 R 库；纯数据包不需要 gcc）

R CMD INSTALL org.Hs.eg.db_3.22.0_src_all.tar.gz GO.db_3.22.0_src_all.tar.gz

# 其余依赖（clusterProfiler/DOSE/enrichplot 等）用 conda install 正常装，

# 若同样卡死，可先装好除这两个数据包外的全部依赖，再按上面源码包方式补装
```

> 安装后验证：`Rscript -e 'library(clusterProfiler); library(org.Hs.eg.db); cat("OK")'`。
> 面试点：能解释 ORA（超几何检验，只看显著基因）就够；被问 GSEA 时答 "GSEA 用全部基因的排序信息、不丢不显著基因，是 ORA 的补充，我这次用的是 ORA"。



***

## 10. 跑通后的自检 + 面试谈资

**自检清单：**



* [ ] `multiqc_report.html` 打开正常，各样本质量达标

* [ ] 每个样本都有 `quant/SRRxxxx/quant.sf`

* [ ] `tx2gene.tsv` 与 `id2name.tsv` 均已生成且非空（行数 = 转录本数 / 基因数）

* [ ] `deseq2_results.txt` 有几千行、`padj` 有值

* [ ] `MAplot_ggplot2.pdf/png`、`volcano_ggplot2.pdf/png` 有红 / 蓝显著点，火山图标注了 top 基因名（如 ZBTB16）

* [ ] `GO_BP_enrichment.txt` 741 行条目、`KEGG_enrichment.txt` 109 行条目，两张气泡图正常（本机实测：GO Top 为激素应答 /actin 骨架，KEGG Top 为黏着斑 / PI3K-Akt，符合地塞米松数据生物学）

* [ ] 能脱稿解释：salmon 拟比对是什么、为什么用 tximport、DESeq2 的 padj 是什么、为什么原始 count 不能直接 t 检验、ORA 富集是什么

**简历 / 面试可以这样写：**

> "使用 salmon + tximport + DESeq2 对真实公开 RNA-seq 数据（GSE52778，人气道平滑肌细胞地塞米松处理，4v4）完成转录本定量与差异表达分析，鉴定出 2140 个显著差异基因（padj<0.05），头号差异基因 ZBTB16 为已知糖皮质激素靶点；GO 富集到 741 个生物学过程条目（Top：激素应答、低氧应答、细胞骨架调控），KEGG 富集到 109 条通路（Top：黏着斑、PI3K-Akt 信号通路），结果与地塞米松已知生物学一致；用 ggplot2 产出发表级 MA / 火山图与富集气泡图；使用 conda 管理环境并配置国内镜像解决公共数据库访问问题。"

**被问 "salmon 和 STAR 的区别" 的标准答法：**

> "STAR 是把 reads 比对到基因组、产出 BAM，适合需要剪接 / 变异分析的场景；salmon 做转录本拟比对直接定量，速度快、内存低，适合以差异表达为目的的分析，两者都是主流方法。我本机内存有限，所以选了 salmon。"



***

## 11. 常见排障速查



| 症状                                                          | 原因                                  | 解法                                                                                            |
| ----------------------------------------------------------- | ----------------------------------- | --------------------------------------------------------------------------------------------- |
| 下载中断                                                        | 网络波动                                | `wget -c` / `aria2c -c` 断点续传，重跑同一命令                                                           |
| fastqc 报错 /multiqc 少数据                                      | fastq.gz 下载不完整（尾部 gzip 块丢失）         | `gzip -t <文件>` 检查，删除后重新下载                                                                     |
| NCBI/UCSC 打不开                                               | 本机网络到美国源不通                          | 一律用 EBI ENA / EBI Ensembl / NGDC / CNGB                                                       |
| conda 装包慢                                                   | 未配置镜像                               | 用第 1 节清华镜像 `.condarc`                                                                         |
| salmon index 被杀（Killed）                                     | 内存不足                                | 关掉其他程序；或 .wslconfig 提 WSL 内存到 12GB                                                            |
| tximport 报 "requires package jsonlite"                      | 读取 salmon 推断重复需要 jsonlite           | `tximport(..., dropInfReps=TRUE)` 跳过                                                          |
| tximport 报 "None of the transcripts ... present in tx2gene" | 转录本 ID 版本后缀不一致（ENST...N vs ENST...） | `tximport(..., ignoreTxVersion=TRUE)`                                                         |
| tximport 报错 / 基因数不对                                         | tx2gene.tsv 有问题                     | 重跑第 7 步 Python 脚本，确认有输出行数                                                                     |
| DESeq2 报错列不匹配                                               | counts 与 coldata 顺序不一致              | 检查 `samples` 与 `coldata.txt` 行顺序完全一致                                                          |
| R 报缺 ggplot2/ggrepel                                        | 未安装                                 | `conda install -y r-ggplot2 r-ggrepel`                                                        |
| R 报缺 clusterProfiler/org.Hs.eg.db                           | 未安装                                 | `conda install -y bioconductor-clusterprofiler bioconductor-org.hs.eg.db`；若事务卡死，按第 9 节源码包方式补装 |
| enrichKEGG 报错或空白                                            | kegg.jp 连不上                         | 网络可达时用 KEGG；不通就跳过 KEGG，GO 结果已够用                                                               |