##prepare results in shell
##awk '{print $2,$1,$3,$9}' 546_gcta_sex_pc1-3_result.mlma > 546_gcta_sex_pc1-3_result.txt #准备为cmplot绘图的格式
##输入文件格式为SNPID CHR POS P
##在shell中去除-nan的数值，并按照P值从小到大排列
awk '$4 != "nan" && $4 != "-nan" && $4 != "NA"' 546_gcta_sex_pc1-3_result.txt | sort -k4,4g > 546_gcta_sex_pc1-3_result_sorted.txt


setwd("F:/caas/毕业课题/第四章_GWAS/horse/852/")
library(CMplot)

# 1. 读取数据
df <- read.delim("825_gcta_sex_pc1-5_loco_result_sorted.txt",
                 header = FALSE, sep = "", stringsAsFactors = FALSE)

# 2. 命名列（注意顺序要和你的文件一致）
colnames(df) <- c("SNP", "Chromosome", "Position", "P.value")

# 3. 选择是否使用默认Bonferroni阈值或自定义阈值
use_bonferroni <- TRUE  # 设为TRUE表示使用Bonferroni阈值，设为FALSE表示自定义阈值

# 4. 如果选择使用Bonferroni，计算Bonferroni阈值；否则，使用自定义阈值
if (use_bonferroni) {
  threshold_value <- 0.05 / nrow(df)
  cat("Using Bonferroni threshold. Bonferroni threshold =", threshold_value, "\n")
} else {
  # 设定自定义阈值（请根据需要修改此处）
  threshold_value <- 0.01  # 自定义的P值阈值，例如：0.01
  cat("Using custom threshold. Custom threshold =", threshold_value, "\n")
}

# 5. Manhattan plot 带阈值线
CMplot(df,
       plot.type = "m",
       threshold = threshold_value,
       threshold.col = "red",     # 阈值线颜色
       threshold.lty = 2,         # 虚线
       threshold.lwd = 1.5,       # 线宽
       LOG10 = TRUE,              # -log10(P)
       amplify = TRUE,            # 高亮显著点
       file.output = TRUE,        # 输出图像文件
       file = "jpg",              # 输出格式
       dpi = 300,
       main = "GCTA GWAS Manhattan Plot")

# 6. QQ plot
CMplot(df,
       plot.type = "q",
       LOG10 = TRUE,
       main = "GCTA GWAS QQ Plot",
       file.output = TRUE,
       file = "jpg",
       dpi = 300)

# === 计算 λ（genomic inflation factor）===
p <- as.numeric(df$P.value)
p <- p[is.finite(p) & !is.na(p)]
# 处理 p=0（避免 qchisq(1-p,1) 变成 Inf）
p[p == 0] <- 1/(length(p) + 1)

chisq <- qchisq(1 - p, df = 1)
lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
cat("Genomic inflation factor (lambda) =", round(lambda, 3), "\n")
