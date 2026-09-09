#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 从 Ensembl GTF 生成 tx2gene.tsv（转录本→基因）与 id2name.tsv（基因→基因名）
# 用法: cd rnaseq 项目目录后运行  python3 make_tx2gene.py
# 输入: ref/Homo_sapiens.GRCh38.116.gtf  输出: tx2gene.tsv / id2name.tsv

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

