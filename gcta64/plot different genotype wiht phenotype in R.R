# 加载必要的包
library(readxl)
library(ggplot2)
library(dplyr)
library(ggsignif)

# 读取Excel文件中的指定sheet
geno_data <- read_excel("your_file.xlsx", sheet = "6_82664891")

# 选择需要的列:第一列(个体ID)和第8列(基因型)
geno_data <- geno_data[, c(1, 8)]
colnames(geno_data) <- c("ID", "Genotype")

# 将基因型转换为因子,确保顺序正确
geno_data$Genotype <- factor(geno_data$Genotype, levels = c(0, 1, 2))

# 读取表型文件
tryCatch({
  pheno_data <- read.table("825_pheno.txt", header = TRUE, comment.char = "")
}, error = function(e) {
  message("使用readLines方法读取文件...")
  lines <- readLines("825_pheno.txt")
  lines <- lines[!grepl("^#", lines)]
  temp_file <- tempfile()
  writeLines(lines, temp_file)
  pheno_data <- read.table(temp_file, header = TRUE)
  unlink(temp_file)
})

# 设置列名
if(ncol(pheno_data) >= 3) {
  colnames(pheno_data) <- c("ID", "V2", "Phenotype")
} else {
  stop("表型文件列数不足,请检查文件格式")
}

# 检查数据结构
print("基因型数据前几行:")
print(head(geno_data))
print("表型数据前几行:")
print(head(pheno_data))

# 合并基因型和表型数据
merged_data <- merge(geno_data, pheno_data, by = "ID")

# 移除缺失值
merged_data <- merged_data[complete.cases(merged_data), ]

# 检查合并后的数据
print("合并后数据前几行:")
print(head(merged_data))
print("各基因型样本数:")
print(table(merged_data$Genotype))

# ========== 统计分析部分 ==========

# 1. 计算每种基因型的平均值和标准差
mean_heights <- merged_data %>%
  group_by(Genotype) %>%
  summarise(
    Mean_Height = mean(Phenotype, na.rm = TRUE),
    SD_Height = sd(Phenotype, na.rm = TRUE),
    N = n()
  )

print("各基因型统计描述:")
print(mean_heights)

# 2. 进行ANOVA分析
anova_result <- aov(Phenotype ~ Genotype, data = merged_data)
anova_summary <- summary(anova_result)
print("ANOVA结果:")
print(anova_summary)

# 提取F值和p值
f_value <- anova_summary[[1]]$`F value`[1]
p_value <- anova_summary[[1]]$`Pr(>F)`[1]

# 3. 进行事后两两比较(Tukey HSD检验)
tukey_result <- TukeyHSD(anova_result)
print("Tukey HSD 事后检验结果:")
print(tukey_result)

# 提取两两比较的p值
p_01 <- tukey_result$Genotype["1-0", "p adj"]  # 基因型0 vs 1
p_02 <- tukey_result$Genotype["2-0", "p adj"]  # 基因型0 vs 2
p_12 <- tukey_result$Genotype["2-1", "p adj"]  # 基因型1 vs 2

# 4. 定义函数将p值转换为显著性标记
get_significance_label <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# 转换p值为显著性标记
sig_01 <- get_significance_label(p_01)
sig_02 <- get_significance_label(p_02)
sig_12 <- get_significance_label(p_12)

print(sprintf("基因型0 vs 1: p = %.4f (%s)", p_01, sig_01))
print(sprintf("基因型0 vs 2: p = %.4f (%s)", p_02, sig_02))
print(sprintf("基因型1 vs 2: p = %.4f (%s)", p_12, sig_12))

# ========== 绘图部分 ==========

# 计算显著性标记的位置
y_max <- max(merged_data$Phenotype, na.rm = TRUE)
y_min <- min(merged_data$Phenotype, na.rm = TRUE)
y_range <- y_max - y_min
y_step <- y_range * 0.05

# 绘制箱线图散点图,并添加显著性标记
p <- ggplot(merged_data, aes(x = Genotype, y = Phenotype)) +
  geom_boxplot(aes(fill = Genotype), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = Genotype), width = 0.2, alpha = 0.6, size = 1.5) +
  geom_point(data = mean_heights, 
             aes(x = Genotype, y = Mean_Height), 
             color = "red", shape = 18, size = 4) +
  
  # 添加显著性标记 - 使用计算得到的显著性
  geom_signif(
    comparisons = list(c("0", "1")),
    annotations = sig_01,
    y_position = y_max + y_step * 2,
    tip_length = 0.01,
    vjust = 0.5
  ) +
  geom_signif(
    comparisons = list(c("0", "2")),
    annotations = sig_02,
    y_position = y_max + y_step * 3,
    tip_length = 0.01,
    vjust = 0.5
  ) +
  geom_signif(
    comparisons = list(c("1", "2")),
    annotations = sig_12,
    y_position = y_max + y_step * 1,
    tip_length = 0.01,
    vjust = 0.5
  ) +
  
  labs(title = "三种基因型的体高分布",
       x = "基因型",
       y = "体高",
       subtitle = "红点表示每种基因型的平均值；***: p < 0.001, **: p < 0.01, *: p < 0.05, ns: 不显著") +
  scale_fill_manual(values = c("0" = "#E69F00", "1" = "#56B4E9", "2" = "#009E73")) +
  scale_color_manual(values = c("0" = "#E69F00", "1" = "#56B4E9", "2" = "#009E73")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

print(p)

# 在图中添加统计结果文本 - 使用计算得到的统计值
p_with_text <- p + 
  annotate("text", x = 1.5, y = y_max + y_step * 4.5, 
           label = sprintf("ANOVA: F = %.2f, p %s", 
                          f_value, 
                          ifelse(p_value < 0.001, "< 0.001", 
                                sprintf("= %.3f", p_value))), 
           size = 3, fontface = "italic")

print(p_with_text)

# 可选:保存图片
ggsave("genotype_phenotype_plot.png", plot = p_with_text, 
       width = 10, height = 8, dpi = 300)

# 输出完整的统计报告
cat("\n========== 统计分析摘要 ==========\n")
cat(sprintf("样本总数: %d\n", nrow(merged_data)))
cat(sprintf("基因型0样本数: %d, 平均值: %.2f ± %.2f\n", 
            mean_heights$N[1], mean_heights$Mean_Height[1], mean_heights$SD_Height[1]))
cat(sprintf("基因型1样本数: %d, 平均值: %.2f ± %.2f\n", 
            mean_heights$N[2], mean_heights$Mean_Height[2], mean_heights$SD_Height[2]))
cat(sprintf("基因型2样本数: %d, 平均值: %.2f ± %.2f\n", 
            mean_heights$N[3], mean_heights$Mean_Height[3], mean_heights$SD_Height[3]))
cat(sprintf("\nANOVA: F(2, %d) = %.2f, p %s\n", 
            nrow(merged_data) - 3, f_value,
            ifelse(p_value < 0.001, "< 0.001", sprintf("= %.3f", p_value))))
cat("\n两两比较(Tukey HSD):\n")
cat(sprintf("  基因型0 vs 1: p = %.4f (%s)\n", p_01, sig_01))
cat(sprintf("  基因型0 vs 2: p = %.4f (%s)\n", p_02, sig_02))
cat(sprintf("  基因型1 vs 2: p = %.4f (%s)\n", p_12, sig_12))
cat("==================================\n")
