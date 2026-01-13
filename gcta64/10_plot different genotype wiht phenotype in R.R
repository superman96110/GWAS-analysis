setwd("F:/caas/毕业课题/第四章_GWAS/sheep/height/1_178008696/")

# ===================== 加载包 =====================
library(readxl)
library(ggplot2)
library(dplyr)
library(ggsignif)

# ===================== 可配置区域 =====================
# 0/1/2 对应的基因型名称（你可以改这里）
geno_map <- c("0" = "AA", "1" = "AG", "2" = "GG")
geno_levels <- unname(geno_map)  # c("AA","AG","GG")

geno_file  <- "snp.xlsx"
geno_sheet <- "Sheet1"
pheno_file1 <- "664_pheno_sd.txt"
pheno_file2 <- "825_pheno.txt"

out_png <- "genotype_phenotype_plot.png"

# ===================== 读取基因型 =====================
geno_data <- read_excel(geno_file, sheet = geno_sheet)
geno_data <- geno_data[, c(1, 8)]
colnames(geno_data) <- c("ID", "Genotype")

# ---- 更稳健：兼容 0/1/2 或已经是 AA/AG/GG ----
geno_data$Genotype <- as.character(geno_data$Genotype)

# 如果是0/1/2，就映射；如果已经是AA/AG/GG，就直接保留
geno_data$Genotype <- ifelse(
  geno_data$Genotype %in% names(geno_map),
  geno_map[geno_data$Genotype],
  geno_data$Genotype
)

# 只保留我们关心的水平，其余设为NA
geno_data$Genotype <- ifelse(geno_data$Genotype %in% geno_levels, geno_data$Genotype, NA)

# 设置因子水平顺序
geno_data$Genotype <- factor(geno_data$Genotype, levels = geno_levels)

# ===================== 读取表型 =====================
pheno_data <- NULL
tryCatch({
  pheno_data <- read.table(pheno_file1, header = TRUE, comment.char = "")
}, error = function(e) {
  message("读取 ", pheno_file1, " 失败，尝试 ", pheno_file2, " ...")
  lines <- readLines(pheno_file2)
  lines <- lines[!grepl("^#", lines)]
  temp_file <- tempfile()
  writeLines(lines, temp_file)
  pheno_data <<- read.table(temp_file, header = TRUE)
  unlink(temp_file)
})

if (is.null(pheno_data)) stop("表型文件读取失败，请检查文件名和路径。")

if (ncol(pheno_data) >= 3) {
  colnames(pheno_data)[1:3] <- c("ID", "V2", "Phenotype")
} else {
  stop("表型文件列数不足(至少3列)，请检查文件格式。")
}

# Phenotype 强制转数值（防止读成字符）
pheno_data$Phenotype <- suppressWarnings(as.numeric(pheno_data$Phenotype))

# ===================== 合并 & 清洗 =====================
merged_data <- merge(geno_data, pheno_data[, c("ID","Phenotype")], by = "ID")

merged_data <- merged_data %>%
  filter(!is.na(Genotype), !is.na(Phenotype), is.finite(Phenotype))

# ===================== 统计分析 =====================
mean_heights <- merged_data %>%
  group_by(Genotype) %>%
  summarise(
    Mean_Height = mean(Phenotype, na.rm = TRUE),
    SD_Height   = sd(Phenotype, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

anova_result <- aov(Phenotype ~ Genotype, data = merged_data)
tukey_result <- TukeyHSD(anova_result)
tuk_p <- tukey_result$Genotype[, "p adj"]

# ---- 更稳健：自动识别三组比较对应的名字 ----
get_pair_p <- function(a, b, tuk_p_vec) {
  nm1 <- paste0(b, "-", a)
  nm2 <- paste0(a, "-", b)
  if (nm1 %in% names(tuk_p_vec)) return(tuk_p_vec[nm1])
  if (nm2 %in% names(tuk_p_vec)) return(tuk_p_vec[nm2])
  return(NA)
}

p_AA_AG <- get_pair_p("AA", "AG", tuk_p)
p_AA_GG <- get_pair_p("AA", "GG", tuk_p)
p_AG_GG <- get_pair_p("AG", "GG", tuk_p)

get_significance_label <- function(pv) {
  if (is.na(pv)) return("NA")
  if (pv < 0.001) return("***")
  if (pv < 0.01)  return("**")
  if (pv < 0.05)  return("*")
  return("ns")
}

sig_AA_AG <- get_significance_label(p_AA_AG)
sig_AA_GG <- get_significance_label(p_AA_GG)
sig_AG_GG <- get_significance_label(p_AG_GG)

# ===================== 绘图参数 =====================
y_max <- max(merged_data$Phenotype, na.rm = TRUE)
y_min <- min(merged_data$Phenotype, na.rm = TRUE)
y_range <- y_max - y_min
y_step <- ifelse(y_range == 0, 1, y_range * 0.08)

# 显著性线位置（从低到高）
y1 <- y_max + y_step * 1.0
y2 <- y_max + y_step * 2.0
y3 <- y_max + y_step * 3.0

# 为了避免显著性线被裁掉：给y轴上方留空间
y_top_lim <- y_max + y_step * 3.8

# ===================== 作图（去掉subtitle备注、去掉ANOVA文字）=====================
p <- ggplot(merged_data, aes(x = Genotype, y = Phenotype)) +
  geom_boxplot(aes(fill = Genotype),
               alpha = 0.65, outlier.shape = NA, width = 0.55) +
  geom_jitter(aes(color = Genotype),
              width = 0.18, alpha = 0.35, size = 1.2) +
  geom_point(data = mean_heights,
             aes(x = Genotype, y = Mean_Height),
             color = "red", shape = 18, size = 3.8) +

  # 显著性标记（Tukey HSD）
  geom_signif(comparisons = list(c("AA", "AG")),
              annotations = sig_AA_AG,
              y_position = y1, tip_length = 0.01, textsize = 4) +
  geom_signif(comparisons = list(c("AG", "GG")),
              annotations = sig_AG_GG,
              y_position = y2, tip_length = 0.01, textsize = 4) +
  geom_signif(comparisons = list(c("AA", "GG")),
              annotations = sig_AA_GG,
              y_position = y3, tip_length = 0.01, textsize = 4) +

  labs(
    title = "三种基因型的体高分布",
    x = "基因型",
    y = "体高"
  ) +
  scale_fill_manual(values = c("AA" = "#E69F00", "AG" = "#56B4E9", "GG" = "#009E73")) +
  scale_color_manual(values = c("AA" = "#E69F00", "AG" = "#56B4E9", "GG" = "#009E73")) +
  coord_cartesian(ylim = c(y_min, y_top_lim)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank()
  )

print(p)

ggsave(out_png, plot = p, width = 10, height = 8, dpi = 300)

# ===================== 输出统计报告（保留在控制台，不画在图里）=====================
cat("\n========== 统计分析摘要 ==========\n")
cat(sprintf("样本总数: %d\n", nrow(merged_data)))
cat("各基因型样本数:\n")
print(table(merged_data$Genotype))

cat("\n两两比较(Tukey HSD):\n")
cat(sprintf("  AA vs AG: p = %s (%s)\n", format(p_AA_AG, digits = 4), sig_AA_AG))
cat(sprintf("  AA vs GG: p = %s (%s)\n", format(p_AA_GG, digits = 4), sig_AA_GG))
cat(sprintf("  AG vs GG: p = %s (%s)\n", format(p_AG_GG, digits = 4), sig_AG_GG))
cat("==================================\n")
