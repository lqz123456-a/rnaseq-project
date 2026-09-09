# RNA-seq 差异表达分析项目（salmon + DESeq2）

真实 RNA-seq 全流程分析项目：8 个人支气管上皮样本（GSE52778 / airway 数据集，地塞米松处理 vs 未处理），从公开数据下载、质控、定量到差异分析与 GO/KEGG 富集，完整跑通。

## 技术栈

| 环节 | 工具 | 说明 |
| --- | --- | --- |
| 数据下载 | wget / ENA（EBI） | 8 个样本单端 50bp FASTQ，`gzip -t` 完整性校验 |
| 质控 | FastQC + MultiQC | 数据质量良好，未修剪 |
| 定量 | salmon | 转录本定量，参考基因组 GRCh38（Ensembl release-116） |
| 差异分析 | tximport + DESeq2 | 负二项分布建模，`padj < 0.05` 判定显著 |
| 可视化 | ggplot2 | MA 图、火山图（含 top 基因标注） |
| 富集分析 | clusterProfiler | GO BP / KEGG 通路富集（ORA） |

## 目录结构

```
rnaseq-project/
├── RNA-seq_salmon全流程_操作清单.md   # 全流程可照做复现文档（含排障表）
├── scripts/
│   ├── make_tx2gene.py                # 从 GTF 生成 tx2gene.tsv / id2name.tsv
│   ├── deseq2.R                       # 差异分析（tximport + DESeq2 + 出 MA/火山图）
│   ├── plot_ggplot2.R                 # 独立重画图：读 deseq2_results.txt 重画 MA/火山图（不重跑差异分析）
│   └── enrichment.R                   # GO BP / KEGG 富集 + 气泡图
├── results/
│   ├── coldata.txt                    # 样本分组表（8 样本）
│   ├── tx2gene.tsv                    # 转录本 → 基因 映射（646,577 条）
│   ├── id2name.tsv                    # 基因 → 基因名 映射（43,458 条）
│   ├── deseq2_results.txt             # 全部基因差异分析结果（34,712 基因）
│   ├── GO_BP_enrichment.txt           # GO BP 富集显著条目
│   └── KEGG_enrichment.txt            # KEGG 通路富集显著条目
└── figures/
    ├── MAplot.pdf / .png      # MA 图
    ├── volcano.pdf / .png     # 火山图（标注 top 基因）
    ├── GO_dotplot.pdf / .png          # GO BP 富集气泡图
    └── KEGG_dotplot.pdf / .png        # KEGG 富集气泡图
```

> fastq / ref / salmon_index / quant / qc 等原始数据与中间文件体积过大，不入库（详见 .gitignore）。

## 关键结果（真实运行产出）

- 定量与差异：检测 34,712 个基因；地塞米松处理显著差异（`padj < 0.05`）**2,140** 个，其中上调 383、下调 325
- 头号基因：**ZBTB16**（log2FC = +5.61），与 airway 数据集已知结论一致
- GO BP 富集：741 个显著条目（Top：cellular response to peptide hormone stimulus、regulation of actin filament-based process、response to hypoxia）
- KEGG 通路：109 条显著通路（Top：Focal adhesion、PI3K-Akt signaling pathway）

## 复现

1. 按 `RNA-seq_salmon全流程_操作清单.md` 第 1–6 节准备环境、下载数据、质控、salmon 定量；
2. 运行 `scripts/make_tx2gene.py` 生成基因映射表；
3. 运行 `scripts/deseq2.R`（差异分析，输出结果表与 MA / 火山图）→ `scripts/enrichment.R`（GO/KEGG 富集）；
4. （可选）`scripts/plot_ggplot2.R`：跳过差异分析，直接读 `deseq2_results.txt` 重画 MA / 火山图，适合单独调整图样式；
5. 关键结果写入 `results/`，图输出到 `figures/`。

环境：WSL Ubuntu + conda（rnaseq1 环境），R 4.x + Bioconductor（DESeq2 / tximport / clusterProfiler），salmon 1.x。
