(base) [supeng@jianglin 546]$ plink --bfile 546_filter --pca 10 --chr-set 31 --out 546_filter_pca


PLINK v1.90b3y 64-bit (4 Nov 2015)         https://www.cog-genomics.org/plink2
(C) 2005-2015 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to 546_filter_pca.log.
Options in effect:
  --bfile 546_filter
  --chr-set 31
  --out 546_filter_pca
  --pca 10

1033908 MB RAM detected; reserving 516954 MB for main workspace.
11117484 variants loaded from .bim file.
546 samples (0 males, 0 females, 546 ambiguous) loaded from .fam.
Ambiguous sex IDs written to 546_filter_pca.nosex .
546 phenotype values loaded from .fam.
Warning: Ignoring phenotypes of missing-sex samples.  If you don't want those
phenotypes to be ignored, use the --allow-no-sex flag.
Using up to 143 threads (change this with --threads).
Warning: This run includes BLAS/LAPACK linear algebra operations which
currently disregard the --threads limit.  If this is problematic, you may want
to recompile against single-threaded BLAS/LAPACK.
Before main variant filters, 546 founders and 0 nonfounders present.
Calculating allele frequencies... done.
11117484 variants and 546 samples pass filters and QC.
Note: No phenotypes present.
Relationship matrix calculation complete.
--pca: Results saved to 546_filter_pca.eigenval and 546_filter_pca.eigenvec .
