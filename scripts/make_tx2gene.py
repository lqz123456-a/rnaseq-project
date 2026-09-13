#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 从 Ensembl GTF 生成 tx2gene.tsv（转录本->基因）与 id2name.tsv（基因->基因名）
# 输出位置固定为项目根目录下的 results/

import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
GTF_PATH = PROJECT_ROOT / "ref" / "Homo_sapiens.GRCh38.116.gtf"
RESULTS_DIR = PROJECT_ROOT / "results"
TX2GENE_PATH = RESULTS_DIR / "tx2gene.tsv"
ID2NAME_PATH = RESULTS_DIR / "id2name.tsv"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

tx2g = {}
id2name = {}

with GTF_PATH.open(encoding="utf-8") as handle:
    for line in handle:
        if line.startswith("#") or len(line.split("\t")) < 9:
            continue

        attrs = line.split("\t")[8]

        if "\ttranscript\t" in line:
            gene = re.search(r'gene_id "([^"]+)"', attrs)
            transcript = re.search(r'transcript_id "([^"]+)"', attrs)
            if gene and transcript:
                tx2g[transcript.group(1)] = gene.group(1)
        elif "\tgene\t" in line:
            gene = re.search(r'gene_id "([^"]+)"', attrs)
            name = re.search(r'gene_name "([^"]+)"', attrs)
            if gene and name:
                id2name[gene.group(1)] = name.group(1)

with TX2GENE_PATH.open("w", encoding="utf-8", newline="") as handle:
    for transcript, gene in tx2g.items():
        handle.write(f"{transcript}\t{gene}\n")

with ID2NAME_PATH.open("w", encoding="utf-8", newline="") as handle:
    for gene_id, name in id2name.items():
        handle.write(f"{gene_id}\t{name}\n")

print(f"tx2gene.tsv {len(tx2g)} 条 -> {TX2GENE_PATH}")
print(f"id2name.tsv {len(id2name)} 条 -> {ID2NAME_PATH}")