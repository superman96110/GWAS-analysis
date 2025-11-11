##prepare results in shell
##awk '{print $2,$1,$3,$9}' 546_gcta_sex_pc1-3_result.mlma > 546_gcta_sex_pc1-3_result.txt #准备为cmplot绘图的格式
##输入文件格式为SNPID CHR POS P

library(CMplot)

# 1. 读取数据
df <- read.delim("F:/caas/毕业课题/第四章_GWAS/horse/546_data/546_gcta_sex_pc1-3_result.txt",
                 header = FALSE, sep = "", stringsAsFactors = FALSE)

# 2. 命名列（注意顺序要和你的文件一致）
colnames(df) <- c("SNP", "Chromosome", "Position", "P.value")

# 3. 计算Bonferroni阈值
threshold_value <- 0.05 / nrow(df)
cat("Bonferroni threshold =", threshold_value, "\n")

# 4. Manhattan plot 带阈值线
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

# 5. QQ plot
CMplot(df,
       plot.type = "q",
       LOG10 = TRUE,
       main = "GCTA GWAS QQ Plot",
       file.output = TRUE,
       file = "jpg",
       dpi = 300)
