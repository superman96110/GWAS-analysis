## R code 
## 默认是基于品种的平均数据的正态分布的表型数值，随机挑选100次，只要挑选的某一次在平均值上下5%内，及赋值给个体
##ind.txt
##FID	ID	Father	Mother	Sex.x	Variety	Sex.y
##A543	A543	0	0	1	Arab	female
##A544	A544	0	0	0	Arab	female
##A548	A548	0	0	0	Arab	male
##A555	A555	0	0	1	Arab	female
##A556	A556	0	0	0	Arab	female

##breed_avg_height.txt
##Breed	Height(male)	sd	Height(female)	sd
##Akhal-Teke	160		145	
##American Saddlebred	163		152	
##Arab	146.2		141.1	
##Chakouyi	134.5	3.72	130.57	6.52
##Dali	121.18	3.25	118.31	3.94
##Debao	97.42	3.76	98.35	4.55
##Falabella	80		70	

## --------------------------------------------
## 路径设置
## --------------------------------------------
setwd("F:/caas/毕业课题/第四章_GWAS/horse/pheno/")

## --------------------------------------------
## 读入数据（修复 # 注释问题）
## --------------------------------------------
ind <- read.table("ind.txt",
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                  comment.char = "")

bh  <- read.table("breed_avg_height.txt",
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                  check.names = FALSE, fill = TRUE, comment.char = "")

## --------------------------------------------
## 参数：当 sd 缺失或不合理时，用 sd = mean * cv_fill
## --------------------------------------------
cv_fill <- 0.05  # 你要求的 5%

num <- function(x) suppressWarnings(as.numeric(x))

## 均值列（必须存在这两列）
male_mean   <- num(bh[["Height(male)"]])
female_mean <- num(bh[["Height(female)"]])

## 若存在 sd 列，则先取其数值；不存在就用 NA 占位
male_sd_raw   <- if ("sd"   %in% names(bh))   num(bh[["sd"]])   else rep(NA_real_, nrow(bh))
female_sd_raw <- if ("sd.1" %in% names(bh))   num(bh[["sd.1"]]) else rep(NA_real_, nrow(bh))

## 对 NA 或 ≤0 的 sd 逐行用 mean * cv_fill 填充
male_sd   <- ifelse(is.na(male_sd_raw)   | male_sd_raw   <= 0, male_mean   * cv_fill, male_sd_raw)
female_sd <- ifelse(is.na(female_sd_raw) | female_sd_raw <= 0, female_mean * cv_fill, female_sd_raw)

## 构成长表
height_male <- data.frame(
  Breed     = trimws(bh$Breed),
  Sex       = "male",
  Height    = male_mean,
  Height_sd = male_sd,
  stringsAsFactors = FALSE
)
height_female <- data.frame(
  Breed     = trimws(bh$Breed),
  Sex       = "female",
  Height    = female_mean,
  Height_sd = female_sd,
  stringsAsFactors = FALSE
)
height_long <- rbind(height_male, height_female)

## --------------------------------------------
## 整理个体表：用 Variety(品种) + Sex.y(性别: male/female)
## --------------------------------------------
trimlow <- function(x) tolower(trimws(x))

if (!("Variety" %in% names(ind))) stop("ind.txt 中未找到列 'Variety'。")
if (!("Sex.y"   %in% names(ind))) stop("ind.txt 中未找到列 'Sex.y'。")

indlist <- data.frame(
  SampleID = ind$ID,
  Breed    = trimws(ind$Variety),
  Sex      = trimlow(ind$Sex.y),
  stringsAsFactors = FALSE
)

## --------------------------------------------
## 合并（左连接）
## --------------------------------------------
merged_df <- merge(indlist, height_long, by = c("Breed","Sex"), all.x = TRUE)

## --------------------------------------------
## 拒绝采样：在 mean ± sd 范围内抽样；最多试 100 次，否则回到 mean
## --------------------------------------------
simulate_value <- function(mean, sd) {
  if (is.na(mean)) return(NA_real_)
  if (is.na(sd) || sd <= 0) return(round(mean, 2))  # 双保险
  for (i in 1:100) {
    val <- rnorm(1, mean, sd)
    if (val >= (mean - sd) && val <= (mean + sd)) return(round(val, 2))
  }
  round(mean, 2)
}

set.seed(123)

## --------------------------------------------
## 生成模拟身高
## --------------------------------------------
merged_df$Sim_Height <- mapply(simulate_value, merged_df$Height, merged_df$Height_sd)

## --------------------------------------------
## 导出结果
## --------------------------------------------
out <- merged_df[, c("SampleID","Breed","Sex","Sim_Height")]
write.table(out, "ind_pheno_height_sim.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

## 简单自检：看看 Franches-Montagnes / male 的前几行是否不再全是 160
subset_check <- subset(out, Breed == "Franches-Montagnes" & Sex == "male")
print(head(subset_check, 10))
