#准备GWAS P值文件
awk 'NR==1{print "SNP","P"; next} {print $2, $11}' /data/supeng/bodysize/horse/gwas/852/825_maf005/output/825_filter_maf005_geno01_mind01_gemma_pc1-5_sex_results.assoc.txt > gwas_for_magma.pval

#准备gene_loc.txt文件，把bed转为MAGMA需要的格式
awk '{print $4, $1, $2, $3}' /data/supeng/bodysize/horse/gwas/852/825_maf005/output/horse3-gene_nochr.bed > gene_loc.txt


#去除gene_loc.txt的重复选项，重新生成文件
awk '
{
  gene=$1; chr=$2; start=$3; end=$4;
  if(!(gene in minS) || start < minS[gene]) minS[gene]=start;
  if(!(gene in maxE) || end   > maxE[gene]) maxE[gene]=end;
  if(!(gene in CHR)) CHR[gene]=chr;
  else if(CHR[gene]!=chr) CHR[gene]="MULTI";
}
END{
  for(g in CHR){
    if(CHR[g]!="MULTI") print g, CHR[g], minS[g], maxE[g];
  }
}
' gene_loc.txt | sort -k2,2n -k3,3n > gene_loc.unique.txt




#准备snp_loc.txt文件，把bim转为MAGMA需要的格式
awk '{print $2, $1, $4}' /data/supeng/bodysize/horse/gwas/852/825_maf005/825_filter_maf005_geno01_mind01.bim > snp_loc_chr.txt

#将SNP与gene映射
/data/supeng/software/MAGMA/magma --annotate nonhuman  --snp-loc snp_loc.txt   --gene-loc gene_loc.unique.txt  --out horse_height


#进行gene-based GWAS计算
/data/supeng/software/MAGMA/magma --bfile /data/supeng/bodysize/horse/gwas/852/825_maf005/825_filter_maf005_geno01_mind01 --pval gwas_for_magma.pval N=825 --gene-annot horse_height.genes.annot --out horse_height_gene
