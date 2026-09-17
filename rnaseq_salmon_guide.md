# RNA-seq 完整复现指南（salmon + DESeq2）

本指南用于复现 GSE52778「airway」数据集的处理流程：从 ENA 下载 paired-end 63 bp FASTQ，使用 FastQC/MultiQC 质控、salmon 定量、tximport 汇总和 DESeq2 差异分析，最后完成 ggplot2 可视化与 GO/KEGG 富集。

> 数据集包含 8 个气道平滑肌细胞样本，即 4 对细胞系的处理与对照。`deseq2.R` 使用 `design=~cell + condition` 配对模型：先扣除 4 个细胞系（N61311、N052611、N080611、N061011）的基线表达差异，再估计地塞米松处理效应。
>
> 本指南在 conda 命令中显式使用 `conda-forge` 和 `bioconda`。国内网络下载缓慢时可额外配置清华镜像，但安装命令的渠道声明保持不变。

## 目录

- [0. 工作目录约定](#0-工作目录约定)
- [1. 环境准备](#1-环境准备)
- [2. 下载参考转录组和注释](#2-下载参考转录组和注释)
- [3. 构建 salmon 索引](#3-构建-salmon-索引)
- [4. 下载 GSE52778 原始数据](#4-下载-gse52778-原始数据)
- [5. 质控](#5-质控)
- [6. salmon 定量](#6-salmon-定量)
- [7. 生成 tx2gene、id2name 和样本表](#7-生成-tx2geneid2name-和样本表)
- [8. 差异表达分析](#8-差异表达分析)
- [9. GO 和 KEGG 富集](#9-go-和-kegg-富集)
- [10. 可选：单独重画 MA 图和火山图](#10-可选单独重画-ma-图和火山图)
- [11. 自检](#11-自检)
- [12. 常见问题](#12-常见问题)
- [13. 数据来源](#13-数据来源)

## 0. 工作目录约定

以下命令默认从项目根目录执行。示例使用 `~/rnaseq`：

```bash
export PROJECT_ROOT="$HOME/rnaseq"
mkdir -p "$PROJECT_ROOT"
cd "$PROJECT_ROOT"
```

克隆仓库后包含脚本和文档；运行流程后会额外生成：

```text
docs/images/  # README 展示用 PNG，已纳入 Git
scripts/
results/    # 本地结果，不纳入 Git
figures/    # 本地图片，不纳入 Git
ref/        # 不纳入 Git
fastq/      # 不纳入 Git
quant/      # 不纳入 Git
```

所有脚本根据自身位置定位项目根目录。参考文件从 `ref/` 读取，salmon 定量从 `quant/` 读取，样本表从 `results/coldata.txt` 读取；结果表写入 `results/`，图片写入 `figures/`。`results/` 和 `figures/` 仅保留在本地，不会上传到 GitHub；仓库仅保留 `docs/images/` 中用于 README 展示的 PNG 副本。

## 1. 环境准备

创建 conda 环境并安装主要工具：

```bash
conda create -n rnaseq -y -c conda-forge -c bioconda python=3.14 fastqc multiqc fastp salmon aria2
conda activate rnaseq
conda install -y -c conda-forge -c bioconda r-base bioconductor-deseq2 bioconductor-tximport r-ggplot2 r-ggrepel bioconductor-clusterprofiler bioconductor-org.hs.eg.db
```

国内网络下载 conda 包较慢时，可选配置清华镜像。没有网络限制时无需设置：

```bash
conda config --set show_channel_urls yes
```

```yaml
channels:
  - defaults
  - conda-forge
  - bioconda
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
```

如果 `clusterProfiler` 或注释包的 conda 安装长时间无响应，可先安装其余依赖，再按第 9 节的源码包方案补装。

### 本次运行版本快照

| 工具 / 包 | 版本 | 用途 |
| --- | --- | --- |
| FastQC | 0.12.1 | 原始数据质控 |
| MultiQC | 1.35 | 质控报告汇总 |
| fastp | 1.3.6 | 可选修剪，本流程未执行 |
| salmon | 2.7.0 | 转录本定量 |
| aria2 | 1.37.0 | 多线程下载 |
| Python | 3.14.7 | 生成转录本和基因映射 |
| R | 4.5.3 | 差异分析、绘图和富集 |
| DESeq2 | 1.50.2 | 差异表达分析 |
| tximport | 1.38.2 | 转录本汇总到基因 |
| ggplot2 | 4.0.3 | MA 图、火山图和气泡图 |
| ggrepel | 0.9.8 | 火山图标签避让 |
| clusterProfiler | 4.18.4 | GO / KEGG 富集 |
| org.Hs.eg.db | 3.22.0 | 人类基因注释 |
| GO.db | 3.22.0 | GO 注释数据 |

如需固定版本，可在安装命令中显式指定，例如 `salmon=2.7.0` 或 `bioconductor-deseq2=1.50.2`。当前仓库没有附带 `environment.yml` 或 `renv.lock`，以上表格只是实测快照，不代表后续安装会自动得到完全相同的版本；在线资源和富集结果可能随时间变化。

## 2. 下载参考转录组和注释

salmon 只需要 cDNA 序列和 GTF，不需要基因组 FASTA：

```bash
mkdir -p ref
wget -c https://ftp.ebi.ac.uk/ensemblorg/pub/release-116/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -P ref
wget -c https://ftp.ebi.ac.uk/ensemblorg/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz -P ref
gzip -dk ref/Homo_sapiens.GRCh38.116.gtf.gz
```

断点续传也可以使用：

```bash
aria2c -x 16 -s 16 -c -d ref "https://ftp.ebi.ac.uk/ensemblorg/pub/release-116/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
```

参考文件大小：cDNA 压缩包约 **176 MB**；GTF 压缩包约 **135 MB**（解压后约 4.4 GB）。Ensembl release 和下载地址可能调整，重跑前应确认链接仍然有效。

## 3. 构建 salmon 索引

```bash
salmon index \
  -t ref/Homo_sapiens.GRCh38.cdna.all.fa.gz \
  -i salmon_index
```

索引构建通常需要约 4–8 GB 内存，在 7.6 GB 内存的 WSL 环境实测可运行。若进程被系统终止，先关闭其他程序或增加 WSL 可用内存。

## 4. 下载 GSE52778 原始数据

样本为 4 对匹配样本，每个细胞系各有一个 untreated 和一个 treated：

| 配对 | untreated | treated |
| --- | --- | --- |
| 1 | SRR1039508 | SRR1039509 |
| 2 | SRR1039512 | SRR1039513 |
| 3 | SRR1039516 | SRR1039517 |
| 4 | SRR1039520 | SRR1039521 |

从 ENA 下载全部 paired-end FASTQ：

```bash
mkdir -p fastq
cd fastq

for u in \
  "SRR1039508 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508" \
  "SRR1039509 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/009/SRR1039509/SRR1039509" \
  "SRR1039512 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/002/SRR1039512/SRR1039512" \
  "SRR1039513 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/003/SRR1039513/SRR1039513" \
  "SRR1039516 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/006/SRR1039516/SRR1039516" \
  "SRR1039517 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/007/SRR1039517/SRR1039517" \
  "SRR1039520 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/000/SRR1039520/SRR1039520" \
  "SRR1039521 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/001/SRR1039521/SRR1039521"
do
  set -- $u
  sample=$1
  base=$2
  aria2c -x 8 -s 8 -c "${base}_1.fastq.gz"
  aria2c -x 8 -s 8 -c "${base}_2.fastq.gz"
done

for f in *.fastq.gz; do
  gzip -t "$f" && echo "$f OK" || echo "$f 损坏，请删除后重新下载"
done

cd "$PROJECT_ROOT"
```

如果链接需要重新查询：

```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR1039508&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes"
```

## 5. 质控

```bash
mkdir -p qc
fastqc -t 8 fastq/*.fastq.gz -o qc/
multiqc qc/ -o qc/
```

打开 `qc/multiqc_report.html` 查看每个样本的质量指标。本项目的已记录流程未执行修剪；如数据质量需要处理，可使用 fastp，但应重新执行后续定量步骤。

## 6. salmon 定量

```bash
mkdir -p quant

for sample in \
  SRR1039508 SRR1039509 SRR1039512 SRR1039513 \
  SRR1039516 SRR1039517 SRR1039520 SRR1039521
do
  salmon quant \
    -i salmon_index \
    -l A \
    -1 "fastq/${sample}_1.fastq.gz" \
    -2 "fastq/${sample}_2.fastq.gz" \
    -p 8 \
    -o "quant/${sample}"
done
```

每个样本应生成 `quant/SRRxxxx/quant.sf`。

## 7. 生成 tx2gene、id2name 和样本表

仓库脚本从 `ref/Homo_sapiens.GRCh38.116.gtf` 读取注释，并将结果写入 `results/`：

- `tx2gene.tsv`：转录本到基因映射，供 tximport 汇总。
- `id2name.tsv`：Ensembl 基因 ID 到基因名映射，供火山图标注。

```bash
mkdir -p results
cat > results/coldata.txt <<'EOF'
sample condition cell
SRR1039508 untreated N61311
SRR1039509 treated N61311
SRR1039512 untreated N052611
SRR1039513 treated N052611
SRR1039516 untreated N080611
SRR1039517 treated N080611
SRR1039520 untreated N061011
SRR1039521 treated N061011
EOF

python3 scripts/make_tx2gene.py
```

本次运行得到 `tx2gene.tsv` 646,577 条、`id2name.tsv` 43,458 条。

## 8. 差异表达分析

```bash
Rscript scripts/deseq2.R
```

`deseq2.R` 执行以下步骤：

1. 使用 `tximport` 按 `tx2gene.tsv` 将转录本定量汇总到基因级。
2. 从 `results/coldata.txt` 读取样本分组（含细胞系与处理信息）。
3. 使用 `DESeqDataSetFromTximport()` 和 `DESeq()` 完成差异分析。
4. 通过 `results(..., contrast=c("condition", "treated", "untreated"))` 提取处理组相对对照组的结果。
5. 将 `results/deseq2_results.txt`、`figures/MAplot.*` 和 `figures/volcano.*` 写入本地输出目录。

当前脚本使用 `design=~cell + condition` 配对模型。GSE52778 的 4 个细胞系各含一对处理/对照样本，模型中先扣除细胞系基线表达差异，再估计处理效应，符合配对实验设计。

`results/` 和 `figures/` 已加入 `.gitignore`。运行这些脚本不会改变 Git 跟踪状态，也不会把结果上传到 GitHub。

## 9. GO 和 KEGG 富集

```bash
Rscript scripts/enrichment.R
```

脚本从 `deseq2_results.txt` 中提取同时满足 `padj < 0.05` 和 `|log2FC| > 1` 的基因，使用 clusterProfiler 完成 GO BP 和 KEGG ORA，并输出：

- `results/GO_BP_enrichment.txt`
- `results/KEGG_enrichment.txt`
- `figures/GO_dotplot.pdf` / `.png`
- `figures/KEGG_dotplot.pdf` / `.png`

KEGG 分析依赖在线数据库和服务状态，结果可能随时间变化。如果 KEGG 暂时不可用，可以先保留 GO 结果。

### 可选：安装注释包的备用方案

若 conda 安装 `org.Hs.eg.db` 或 `GO.db` 长时间无响应，可从 Bioconductor 源码包安装：

```bash
wget -c https://depot.galaxyproject.org/software/bioconductor-org.hs.eg.db/bioconductor-org.hs.eg.db_3.22.0_src_all.tar.gz
wget -c https://depot.galaxyproject.org/software/bioconductor-go.db/bioconductor-go.db_3.22.0_src_all.tar.gz
R CMD INSTALL org.Hs.eg.db_3.22.0_src_all.tar.gz GO.db_3.22.0_src_all.tar.gz
Rscript -e 'library(clusterProfiler); library(org.Hs.eg.db); cat("OK\n")'
```

下载后应核对官方校验值再安装。

## 10. 可选：单独重画 MA 图和火山图

若只需要调整图形样式，不必重新运行差异分析：

```bash
Rscript scripts/plot_ggplot2.R
```

该脚本读取 `results/deseq2_results.txt` 和 `results/id2name.tsv`，重新生成 `figures/MAplot.*` 和 `figures/volcano.*`。

## 11. 自检

- `qc/multiqc_report.html` 可以正常打开，样本质量指标合理。
- `quant/SRRxxxx/quant.sf` 对 8 个样本均存在。
- `results/tx2gene.tsv` 和 `results/id2name.tsv` 非空。
- `results/deseq2_results.txt` 包含 34,712 个基因的差异分析结果。
- 按 `padj < 0.05` 统计得到 3,387 个显著基因，其中上调 1,843 个、下调 1,544 个。
- 同时要求 `|log2FC| > 1` 时，得到上调 441 个、下调 403 个。
- GO/KEGG 富集分析使用同时满足 `padj < 0.05` 和 `|log2FC| > 1` 的 844 个基因。
- 当前结果中 ZBTB16 的 `log2FC = +5.68`、`padj = 3.79e-130`，按 `padj` 排序为第 1；靠前基因还包括 DUSP1、NEXN 和 SAMHD1。
- GO BP 和 KEGG 结果表及气泡图正常生成。本机运行快照为 GO 611 条、KEGG 29 条，重新运行可能变化。

可使用以下命令检查样本表、映射表、结果规模和主要统计量：

```bash
test -s results/coldata.txt
test -s results/tx2gene.tsv
test -s results/id2name.tsv
wc -l results/tx2gene.tsv results/id2name.tsv results/deseq2_results.txt
Rscript -e 'res <- read.delim("results/deseq2_results.txt", row.names = 1); sig <- !is.na(res$padj) & res$padj < 0.05; effect <- sig & abs(res$log2FoldChange) > 1; cat("genes:", nrow(res), "significant:", sum(sig), "effect:", sum(effect), "up:", sum(sig & res$log2FoldChange > 0), "down:", sum(sig & res$log2FoldChange < 0), "\n")'
```

注意：当前 `deseq2_results.txt` 第一列是基因 ID，但表头没有单独的 `gene_id` 列。项目内 R 脚本按行名读取；使用 pandas 等其他工具时需要通过 `index_col=0` 或等价方式处理。后续可考虑在脚本中显式输出 `gene_id` 列。

## 12. 常见问题

| 症状 | 可能原因 | 处理方式 |
| --- | --- | --- |
| 下载中断 | 网络波动 | 使用 `wget -c` 或 `aria2c -c` 断点续传 |
| FastQC 报错或 MultiQC 样本缺失 | FASTQ gzip 文件不完整 | 使用 `gzip -t` 检查并重新下载 |
| NCBI 或 UCSC 访问不稳定 | 网络路径问题 | 使用 ENA、Ensembl 或本地镜像 |
| conda 安装缓慢 | `conda-forge` 或 `bioconda` 网络较慢 | 可选配置清华镜像，或分步安装依赖 |
| salmon index 被终止 | 内存不足 | 关闭其他程序或增加 WSL 内存 |
| tximport 缺少 jsonlite | 读取 salmon 推断重复信息 | 当前脚本已设置 `dropInfReps=TRUE` |
| tximport 找不到转录本 | 转录本版本后缀不一致 | 当前脚本已设置 `ignoreTxVersion=TRUE` |
| DESeq2 列不匹配 | 定量文件和分组表顺序不一致 | 检查 `coldata.txt` 与 `quant/` 中的样本名 |
| KEGG 为空或报错 | 在线服务不可达 | 保留 GO 结果，网络恢复后重跑 KEGG |
| 脚本找不到参考或定量文件 | `ref/` 或 `quant/` 不完整 | 检查参考文件、salmon 索引和每个样本的 `quant.sf` |

## 13. 数据来源

- GEO：[GSE52778 airway](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778)
- 原始测序数据：EBI ENA，accessions 见第 4 节
- 关联论文：Himes et al. (2014)，[DOI: 10.1371/journal.pone.0099625](https://doi.org/10.1371/journal.pone.0099625)
- 参考转录组和注释：[Ensembl GRCh38 release-116](https://ftp.ensembl.org/pub/release-116/)
- 下载链接和在线数据库状态会变化；本指南记录的是 2026-09 验证时的可用路径。
