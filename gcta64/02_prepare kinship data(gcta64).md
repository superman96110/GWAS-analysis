#使用plink计算IBD，看是否有超高的相似的个体
plink --bfile 546_filter --chr-set 31 --keep-allele-order --genome --out 546_filter
#提取并排序
sort -k10,10n 546_filter.genome > 546_filter_sorted_by_IBD.genome

#将bim文件中第二列的“.”替换为“chr:position”,其余列保持不变
awk -v OFS="\t" '{ $2 = $1 ":" $4; print }' 825_filter_maf001_geno01_mind01.bim > 825_renamed.bim

#从更大的数据中，根据fam文件中提取sex
awk -v OFS='\t' 'NR==FNR{ids[$1]=1; next} ($1 in ids)' \
  825_filter_maf001_geno01_mind01.fam \
  /data/supeng/bodysize/horse/gwas/1160/1160_sex.txt \
  > 825_sex.txt

# 从更大的数据中，根据fam文件中提取pheno
awk -v OFS='\t' 'NR==FNR{ids[$1]=1; next} ($1 in ids)' \
  825_filter_maf001_geno01_mind01.fam \
  /data/supeng/bodysize/horse/gwas/1160/1160_pheno.txt \
  > 825_pheno.txt


## 进行GWAS分析前，需要准备K矩阵和PCA协变量等。
##这里使用的是gcta64计算grm矩阵


(base) [supeng@jianglin 546]$ /home/jianglin/software/gcta64 --bfile 546_filter --autosome --autosome-num 31 --make-grm --thread-num 16 --out 546_filter_grm



*******************************************************************
* Genome-wide Complex Trait Analysis (GCTA)
* version 1.24.7
* (C) 2010-2013 Jian Yang, Hong Lee, Michael Goddard and Peter Visscher
* The University of Queensland
* MIT License
*******************************************************************
Analysis started: Tue Nov 11 12:46:29 2025

Options:
--bfile 546_filter
--autosome
--autosome-num 26
--make-grm
--thread-num 16
--out 546_filter_grm

Note: the program will be running on 16 threads.

Reading PLINK FAM file from [546_filter.fam].
546 individuals to be included from [546_filter.fam].
Reading PLINK BIM file from [546_filter.bim].
11117484 SNPs to be included from [546_filter.bim].
10158075 SNPs from chromosome 1 to chromosome 26 are included in the analysis.
Reading PLINK BED file from [546_filter.bed] in SNP-major format ...
Genotype data for 546 individuals and 10158075 SNPs to be included from [546_filter.bed].
Calculating allele frequencies ...
Recoding genotypes (individual major mode) ...

Calculating the genetic relationship matrix (GRM) ... (Note: default speed-optimized mode, may use huge RAM)

Summary of the GRM:
Mean of diagonals = 1.10836
Variance of diagonals = 0.00567863
Mean of off-diagonals = -0.00203366
Variance of off-diagonals = 0.00398875
GRM of 546 individuals has been saved in the file [546_filter_grm.grm.bin] (in binary format).
Number of SNPs to calcuate the genetic relationship between each pair of individuals has been saved in the file [546_filter_grm.grm.N.bin] (in binary format).
IDs for the GRM file [546_filter_grm.grm.bin] have been saved in the file [546_filter_grm.grm.id].

Analysis finished: Tue Nov 11 12:51:45 2025
Computational time: 0:5:16
