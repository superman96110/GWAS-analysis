(base) [supeng@jianglin 546]$ plink --bfile 546_filter --chr-set 31 --pheno 546_pheno.txt --covar 546_covar_plink.txt --linear hide-covar --threads 16 --out 546_plink_linear_pc1to3_sex
PLINK v1.90b3y 64-bit (4 Nov 2015)         https://www.cog-genomics.org/plink2
(C) 2005-2015 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to 546_plink_linear_pc1to3_sex.log.
Options in effect:
  --bfile 546_filter #plink生成的bfile文件
  --chr-set 31 #染色体数为31，马物种
  --covar 546_covar_plink.txt #协变量
  --linear hide-covar #使用一般线性模型，结果不显示协变量add结果
  --out 546_plink_linear_pc1to3_sex 
  --pheno 546_pheno.txt #表型数据
  --threads 16 #服务器线程为16

1033908 MB RAM detected; reserving 516954 MB for main workspace.
11117484 variants loaded from .bim file.
546 samples (211 males, 335 females) loaded from .fam.
546 phenotype values present after --pheno.
Using 1 thread.
Warning: This run includes BLAS/LAPACK linear algebra operations which
currently disregard the --threads limit.  If this is problematic, you may want
to recompile against single-threaded BLAS/LAPACK.
--covar: 3 covariates loaded.
Before main variant filters, 546 founders and 0 nonfounders present.
Calculating allele frequencies... done.
11117484 variants and 546 samples pass filters and QC.
Phenotype data is quantitative.
Writing linear model association results to
546_plink_linear_pc1to3_sex.assoc.linear ... 32%
