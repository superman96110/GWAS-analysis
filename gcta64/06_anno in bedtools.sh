#按照第四列排序，将小于阈值的SNP全部提出到文件中
awk 'NR==1 || $4 < 1e-7' 825_gcta_sex_pc1-5_loco_result_sorted.txt > sig_p_less_1e7.txt




#anno.sh
#./anno.sh sig_p_less_1e7.txt horse3-gene_nochr.bed


#!/bin/bash

# 使用方法：
#   ./anno.sh GWAS_FILE GENE_FILE
#
# 举例：
#   ./anno.sh 825_gcta_sex_pc1-5_loco_result_sorted.txt genes.bed
#
# 说明：
#   1) 对 $GWAS_FILE 做一次注释，输出：final_annotated_results_all.txt
#   2) 如果当前目录存在 sig_p_less_1e7.txt，再对该文件做一次注释：
#        输入：sig_p_less_1e7.txt
#        输出：final_annotated_results_sig_p_less_1e7.txt

GWAS_FILE="$1"  # 第一个参数是GWAS文件（全量）
GENE_FILE="$2"  # 第二个参数是基因Bed文件

# ---- 基本检查 ----
if [[ -z "$GWAS_FILE" || -z "$GENE_FILE" ]]; then
    echo "用法: $0 GWAS_FILE GENE_FILE"
    exit 1
fi

if [[ ! -f "$GWAS_FILE" ]]; then
    echo "错误: 找不到GWAS文件: $GWAS_FILE"
    exit 1
fi

if [[ ! -f "$GENE_FILE" ]]; then
    echo "错误: 找不到基因文件: $GENE_FILE"
    exit 1
fi

# ---- 定义一个函数：对任意一个 GWAS 文件做注释 ----
annotate_gwas () {
    local GWAS_IN="$1"    # 输入GWAS文件
    local LABEL="$2"      # 用于区分输出文件名的标签

    echo "👉 开始注释 GWAS 文件: $GWAS_IN （标签: $LABEL）"

    # 步骤 1: 准备 GWAS SNP Bed 文件，转换为 Bed 格式并扩展 50KB
    # 假设输入文件格式类似：
    # SNP Chr bp p
    # 6:82664891 6 82664891 2.90322e-14
    awk '
    NR > 1 {
        # 提取 Chr 和 bp：从第一列 SNP（形如 6:82664891）拆分
        split($1, a, ":");
        chr = a[1];   # Chr 是前半部分
        bp  = a[2];   # bp 是后半部分
        snp_id = $1;

        # p 值在第 4 列
        p_value = $4;

        # 计算 50KB 范围，注意: Bed 格式是 0-based 半开区间
        start = bp - 50000;
        end   = bp + 50000;

        # 确保 start 不小于 1
        if (start < 1) {
            start = 1;
        }

        # 输出为 Bed 格式: Chr, Start-1, End, SNP_ID, p
        print chr "\t" (start - 1) "\t" end "\t" snp_id "\t" p_value;
    }
    ' "$GWAS_IN" > "snp_50kb_${LABEL}.bed"

    # 步骤 2: 使用 bedtools intersect 进行重叠分析
    bedtools intersect -a "snp_50kb_${LABEL}.bed" -b "$GENE_FILE" -wao | \
    awk '
    BEGIN {
        OFS = "\t";
    }
    {
        # snp_50kb_*.bed 的列: $1=Chr, $2=Start-1, $3=End, $4=SNP_ID, $5=p
        snp_id     = $4;
        p_value    = $5;

        # 假设基因Bed文件的第 9 列是 GeneName，第 10 列是 Overlap_Length（来自 -wao）
        gene_name  = $9;
        overlap_len = $10;

        # 存储原始数据和所有重叠的基因
        if (! (snp_id in data)) {
            # 存储原始数据 (Chr, bp, p)
            split(snp_id, a, ":");
            chr = a[1];
            bp  = a[2];
            data[snp_id] = chr OFS bp OFS p_value;
            genes[snp_id] = "";
        }

        # 聚合基因名（如果有重叠）
        if (overlap_len > 0) {
            if (genes[snp_id] == "") {
                genes[snp_id] = gene_name;
            } else {
                # 防止重复拼接同一个基因名（简单处理）
                if (index(genes[snp_id], gene_name) == 0) {
                    genes[snp_id] = genes[snp_id] ";" gene_name;
                }
            }
        }
    }
    END {
        print "SNP","Chr","bp","p","Gene_Annotation";
        for (id in data) {
            ann = genes[id];
            if (ann == "") {
                ann = "NA";
            }
            print id, data[id], ann;
        }
    }
    ' > "final_annotated_results_${LABEL}.txt"

    echo "✅ 注释完成！输出文件: final_annotated_results_${LABEL}.txt"
}

# ---- 1) 先对全量 GWAS 文件做注释 ----
annotate_gwas "$GWAS_FILE" "all"

# ---- 2) 如果存在 sig_p_less_1e7.txt，则也做一遍注释 ----
if [[ -f "sig_p_less_1e7.txt" ]]; then
    annotate_gwas "sig_p_less_1e7.txt" "sig_p_less_1e7"
else
    echo "⚠️ 未找到 sig_p_less_1e7.txt，跳过该文件的注释。"
fi
