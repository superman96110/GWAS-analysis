# 1. 读入原始表型
pheno <- read.table(
    "ind_pheno_weight_sim.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    comment.char = ""
)

# 2. 对 Sim_weight 做标准化（z-score）
pheno$PHENO_STD <- as.numeric(scale(pheno$Sim_weight))

# 3. 组织成：ID  ID  标准化表型  三列
out <- data.frame(
    pheno$SampleID,   # 第一列：ID
    pheno$SampleID,   # 第二列：ID（再重复一次）
    pheno$PHENO_STD   # 第三列：标准化后的表型
)

# 4. 输出为 txt：
#    - 无标题（col.names = FALSE）
#    - 制表符分隔
write.table(
    out,
    file = "ind_pheno_weight_sim_std.txt",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

#在linux中需要调整一下文本格式
dos2unix ind_pheno_height_sim_std.txt
