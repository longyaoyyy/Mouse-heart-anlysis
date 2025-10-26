library(ggplot2)


plot_gene_group_anova <- function(mat_row, group_labels, gene_name) {
  # mat_row: 矩阵的一行（某个基因在所有样本中的表达值）
  # group_labels: 向量，长度=ncol(mat)，表示每个样本属于哪个组
  # gene_name: 当前基因的名称（用于标题）
  cols <- c( "#444576", "#4682B4", "#AEDEEE","#FFA500", "#FFD790","#C65762",'#FBDFDE', "#F6EFCF","#BCB99F")
  pal <- colorRampPalette(cols)
  # 将表达值与对应的组别关联起来
  expressions <- as.numeric(mat_row)          # 当前基因在所有样本中的值
  groups <- group_labels                      # 每个样本的组别
  
  # 构造数据框
  df <- data.frame(
    Expression = expressions,
    Group = as.factor(groups)  # 转为因子，用于绘图和ANOVA
  )
  
  # 检查是否有至少2个组
  unique_groups <- unique(df$Group)
  num_groups <- length(unique_groups)
  
  if (num_groups < 2) {
    stop("组数少于2，无法进行方差分析。")
  }
  
  # 进行单因素方差分析（ANOVA）
  anova_result <- aov(Expression ~ Group, data = df)
  p_value <- summary(anova_result)[[1]]$"Pr(>F)"[1]
  
  # 计算每组的均值和标准误
  library(dplyr)
  stats <- df %>%
    group_by(Group) %>%
    summarise(
      Mean = mean(Expression),
      SD = sd(Expression),
      SE = SD / sqrt(n()),
      .groups = 'drop'
    )
  

  # 绘图
  p <- ggplot(stats, aes(x = Group, y = Mean, fill = Group)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = pal(8))+
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
    labs(
      title = paste("Gene:", gene_name),
      x = "",
      y = "FoldChange"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = 'none',
          axis.text.x = element_text(angle = 45,
                                     vjust = 0.5,hjust = 0.5))
  # 🟩 在图上标注 F 值 和 p 值
  p <- p + 
    annotate("text", x = 1, y = max(stats$Mean) + max(stats$SE) + 0.5,
             label = paste0("ANOVA p =", format.pval(p_value, digits = 2, eps = 0.001)),
             hjust = 0, size = 3.5)
  
  # 打印 p 值到控制台（可选）
  cat("—— Gene:", gene_name, "——\n")
  cat("ANOVA P-value:", p_value, "\n\n")
  
  return(p)
}
pdf("Gene_Expression_By_Group_ANOVA_Plots.pdf", width = 8, height = 6)
for (i in 1:nrow(hypertrophy)) {
  gene_name <- rownames(hypertrophy)[i]
  mat_row <- hypertrophy[i, ]
  p <- plot_gene_group_anova(mat_row, group, gene_name)
  print(p)  # 显示图形
}
dev.off()
