#进行SusieR的包进行fine mapping
#需要准备的文件有给定区域的mlma的结果文件，给定区域的bim文件喝给定区域的ld结果文件

#提取区域的bim文件
plink --bfile 825_filter_maf001_geno01_mind01 --chr 6 --from-bp 82513000 --to-bp 82664950 --make-bed --out 825_filter_maf001_geno01_mind01_chr6_82513_82665 --chr-set 31 --keep-allele-order

#计算LD区域
plink --bfile 825_filter_maf001_geno01_mind01_chr6_82513_82665 --r square --chr-set 31 --keep-allele-order --out 825_filter_maf001_geno01_mind01_chr6_82513_82665

#提取区域的mlma的结果文件
awk 'BEGIN{FS=OFS="\t"} NR==1 || ($1 == 6 && $3 >= 82513000 && $3 <= 82664950)' 825_filter_maf001_geno01_mind01_gcta_pc1-5_sex_loco_result.loco.mlma > chr6_82513_82665_mlma.txt
awk '{print $1,$3,$9}' chr6_8248_82665_mlma.txt > chr6_8248_82665_mlma1.txt


library(susieR)
library(data.table)
setwd("F:/caas/毕业课题/第四章_GWAS/horse/852/")
LD <- as.matrix(fread("825_filter_maf001_geno01_mind01_chr6_82513_82665.ld"))
dim(LD)
bim <- fread("825_filter_maf001_geno01_mind01_chr6_82513_82665.bim")
setnames(bim, c("CHR","SNP","CM","BP","A1","A2"))
gwas <- fread("chr6_82513_82665_mlma.txt")
setnames(gwas, c("CHR","SNP","BP","A1","A2","Freq","b","se","p"))
# merge 使得gwas2里包含bim顺序的SNP
gwas2 <- merge(
    bim[, .(CHR, SNP, BP, A1, A2)],
    gwas,
    by = c("CHR","SNP","BP","A1","A2")
)

# 根据bim SNP顺序排序
gwas2 <- gwas2[match(bim$SNP, gwas2$SNP), ]
stopifnot(
    nrow(gwas2) == nrow(LD),
    nrow(LD) == ncol(LD)
)
Z <- gwas2$b / gwas2$se
n <- 852
fit <- susie_rss(
    z = Z,
    R = LD,
    n = n,
    L = 10    # 最多10个独立因果信号，可改大可改小
)
pip <- fit$pip

result_susie <- data.table(
    Chrom   = gwas2$CHR,
    Position= gwas2$BP,
    Name    = gwas2$SNP,
    PIP     = pip
)

fwrite(result_susie, "chr6_82513_82665_susie_PIP.txt", sep = "\t")
fit$sets$cs
susie_plot(fit, y="PIP")

#提取某个位点的基因分型
plink --bfile 825_filter_maf001_geno01_mind01 --snp 6:82610088 --recode A --chr-set 31 --out 6_82610088
       
#LDblockshow
/angr/wangsn/jyz/szd_data/ld/LDBlockShow/bin/LDBlockShow -InVCF 825_filter_maf001_geno01_mind01_chr6_82513_82665.vcf.gz -OutPut out -InGWAS chr6_82513_82665_mlma.txt -InGFF genes_mRNA.gff -SpeSNPName hmga2.txt -Region 6:82513010:82664922  -OutPng -SeleVar 3 -TopSite
