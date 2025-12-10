#读取的表型文件为三列数据
#前两列分别为FID和IID，最后一列为表型的txt文件


# ===============================
# 读取数据（重要：禁用注释符号）
# ===============================
df <- read.table("825_pheno.txt",
                 header = FALSE,
                 comment.char = "",
                 sep = "",
                 stringsAsFactors = FALSE)

# 第三列为表型
phenotype <- df$V3

# ===============================
# 1. 绘制表型分布直方图
# ===============================
hist_res <- hist(phenotype,
                 breaks = 30,
                 main = "Phenotype Distribution",
                 xlab = "Phenotype Value",
                 ylab = "Frequency",
                 col = "skyblue",
                 border = "white")

# 添加拟合正态分布曲线
x <- seq(min(phenotype), max(phenotype), length = 200)
y <- dnorm(x, mean = mean(phenotype), sd = sd(phenotype))

# 将密度转换到直方图的频数尺度
y_scaled <- y * length(phenotype) * diff(hist_res$breaks)[1]

lines(x, y_scaled, col = "red", lwd = 2)

# ===============================
# 2. 正态性检验
# ===============================
cat("\n===== Shapiro-Wilk Normality Test =====\n")
print(shapiro.test(phenotype))

# ===============================
# 3. Q-Q Plot（正态性检查辅助图）
# ===============================
qqnorm(phenotype,
       main = "Q-Q Plot of Phenotype")
qqline(phenotype, col = "red", lwd = 2)
