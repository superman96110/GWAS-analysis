#-lmm参数解释
#-lmm 1	Wald test	使用基于最大似然估计的固定效应检验	⭐ 速度最快，适合大样本 GWAS
#-lmm 2	Likelihood Ratio Test (LRT)	基于似然比的严格检验	计算更耗时，结果更保守（P值略大）
#-lmm 3	Score test	不需要估计 β 值，基于分数统计量	速度快、稳健，对小样本有优势
#-lmm 4	All three (Wald + LRT + Score)	同时输出三种检验结果	🧪 最全面但最耗时；常用于比较验证


(base) [supeng@jianglin 546]$ gemma -bfile 546_filter -k /data/supeng/bodysize/horse/gwas/546/output/546_filter.cXX.txt -c 546_gemma_c.txt -p 546_gemma_pheno.txt -lmm 1 -o 546_filter_gemma_pc1-3_sex


Reading Files ...
no intecept term is found in the cvt file. a column of 1s is added.
## number of total individuals = 546
## number of analyzed individuals = 546
## number of covariates = 5
## number of phenotypes = 1
## number of total SNPs = 11117484
## number of analyzed SNPs = 11113829
Start Eigen-Decomposition...
pve estimate =0.871644
se(pve) =0.0310082
loading SNPs                                                    0.90%
