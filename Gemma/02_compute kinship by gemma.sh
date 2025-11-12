#gk参数解释
#-gk 0	标准化的相关系数矩阵（legacy）	最早版本方法，不常用	⚠️ 已较少使用
#-gk 1	中心化（centered）关系矩阵	对基因型数据中心化，不进行标准化	⭐ 最常用、GCTA兼容
#-gk 2	标准化（standardized）关系矩阵	对基因型中心化并按方差标准化	适用于混合群体或等位基因频率差异大
#-gk 3	非标准化关系矩阵（原始）	不进行中心化和标准化	仅用于特殊分析，不推荐

(base) [supeng@jianglin 546]$ gemma -bfile 546_filter -gk 1 -o 546_filter


Reading Files ...
## number of total individuals = 546
## number of analyzed individuals = 546
## number of covariates = 1
## number of phenotypes = 1
## number of total SNPs = 11117484
## number of analyzed SNPs = 11113829
Calculating Relatedness Matrix ...
Reading SNPs  ==================================================100.00%
