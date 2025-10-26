library(Seurat)
library(ggplot2)
library(ggunchained) 
library(reshape2) 
library(SCP)
set.seed(4180)
setwd("E:\\2.工作/YF/mHeart/")
#########color
cols <- c( "#444576", "#4682B4","#87CEEB","#FFD700", "#FFA500", "#C65762","#4D4D4D", "#8C4F00" ,"#BCA88F")
pal <- colorRampPalette(cols)
Mac <-  readRDS("E:/2.工作/YF/mHeart/Mac.Rds")
obj <- subset(Mac,subtype %in% c('Macrophage-1'))
#########miloR######
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
DefaultAssay(Mac) <- 'RNA'
Mac@assays <- Mac@assays['RNA']
traj_milo_sce <- as.SingleCellExperiment(Mac)
traj_milo <- Milo(traj_milo_sce)
traj_milo <- buildGraph(traj_milo)
traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = 10, d=30, refined = TRUE)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="orig.ident")
traj_milo@colData$preg <- ifelse(traj_milo@colData$group %in% c('NP','PP'),'Non-Preg.','Preg.')
traj_design <- data.frame(colData(traj_milo))[,c("orig.ident", "preg")]
traj_design <- distinct(traj_design)
rownames(traj_design)<- c(1:12)
traj_milo <- calcNhoodDistance(traj_milo, d=30)
row.names(traj_design) <- traj_design$orig.ident
da_results <- testNhoods(traj_milo, design = ~preg, design.df = traj_design)
da_results %>%
  arrange(- SpatialFDR) %>%
  head()
traj_milo <- buildNhoodGraph(traj_milo)
p <- plotNhoodGraphDA(traj_milo, da_results, alpha=0.05) +scale_size(
  range = c(0, 10))+
  plot_layout(guides="collect")
p

############scdist######
library(dplyr)
library(scDist)
library(ggplot2)
library(ggpubr)
##run Dist
sim <- list(Y=Mac@assays$RNA$scale.data %>% as.data.frame(),
            meta.data=Mac@meta.data %>% as.data.frame())
out <- scDist(normalized_counts = sim$Y, # 标准化的数据矩阵
              meta.data = sim$meta.data, # metadata表格
              d = 20, # 指定用于PCA分析的维度数量为前20
              fixed.effects = "group", # 你感兴趣的分组条件，对应metadata表格中的列名称
              random.effects = c('orig.ident'),  # 需要去除的潜在的影响因素，比如不同的样本、年龄、批次等等，对应metadata表格中的列名称
              clusters="subtype" # 待分析的细胞类型，对应metadata表格中的列名称
)
DistPlot(out)+ theme_bw()
#plot
distGenes(out, cluster = "Mac-2") # 指定4这一群细胞
df <- data.frame(value = out$vals[["Mac-2"]]$beta.hat, 
                 label = out$gene.names) %>% top_n(100, abs(value)) # 可视化前20个
df$color <- ifelse(df$value>0, "Positive", "Negative")
ggbarplot(df,
          x="label", 
          y="value", 
          fill = "color",
          color = "white",
          palette = "npg",
          sort.val = "asc",
          sort.by.groups = FALSE,
          xlab = "",
          legend.title = "")+ theme_bw() + ylab("Condition difference") + coord_flip()
################通路富集########
library(dplyr)
library(BiocParallel)
library(clusterProfiler)
library(org.Mm.eg.db)
register(BPPARAM = SnowParam(1))
samplegroup <-c('NP','MP','LP','PP')
obj <- subset(Mac,subtype %in% c('Macrophage-1'))
obj <- RunDEtest(srt = obj , group_by = "group", assay = 'aUCell',
                 fc.threshold = 1, only.pos = FALSE,min.pct = 0.25,
                 group1 = 'NP')
DEGs_group <- obj@tools$DEtest_custom$AllMarkers_wilcox
for (j in 2:4){
  obj<- RunDEtest(srt = obj, group_by = "group", assay = 'aUCell',
                  fc.threshold = 1, only.pos = FALSE,min.pct = 0.25,
                  group1 = samplegroup[j],group2 = 'NP')
  DEGs_group <- rbind(DEGs_group,obj@tools$DEtest_custom$AllMarkers_wilcox)
}
obj@tools$DEtest_group$AllMarkers_wilcox <- DEGs_group
top_genes <- DEGs_group %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  group_by(group1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 60)
df_sig <- top_genes
group <- data.frame(gene=df_sig$gene,
                    group=df_sig$group1)
gene_ID <- bitr(group$gene, fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Mm.eg.db")
data  <- merge(gene_ID ,group,by.x='SYMBOL',by.y='gene')
levels(data$group) <- c('NP','MP','LP','PP')
data_GO <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Mm.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)
res <- data_GO@compareClusterResult
library(dplyr)
for (i in 1:dim(res)[1]) {
  arr = unlist(strsplit(as.character(res[i,"geneID"]), split="/"))
  gene_names = paste(unique(data$SYMBOL[data$ENTREZID %in% arr]), collapse="/")
  res[i,"geneID"] = gene_names
}
writexl::write_xlsx(res,'FigS9i_Mac1_四组GO通路富集结果.xlsx')
dotplot(data_GO, showCategory= 5,font.size = 10)

###########promoter ATAC######
library(ggrastr)
setwd("E:/2.工作/YF/mHeart/8.终稿调整/2050314/figS9_Mac")
avdata <- AggregateExpression(obj,group.by = 'group',assays = 'activities',slot = 'counts')
avdata <- as.data.frame(avdata$activities)
for(j in c('MP','LP','PP')){
  avdata_df <- data.frame(NP=avdata$NP,test=avdata[,j])
  p <- ggplot(avdata_df , aes(x = log10(NP), y = log10(test))) +
    geom_point_rast(size = 0.001)  +
    # 添加 x = y 斜线
    geom_abline(intercept = 0, slope = 1, color = "#E83945") +
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    coord_fixed()+
    labs(x = "log10(gene activity) of NP",
         y = paste("log10(gene activity) of ",j,sep = ''),
         title = i)
  ggsave(paste(i,'_NPvs',j,'所有细胞染色质开放片段数散点图.pdf',sep = ''),p,
         width = 3.84,height = 3.68)
}
ratiodt <- data.frame(
  Column = c( "MP",'LP','PP'),  # 列名称
  Number =c(sum(avdata[, 2] > avdata[, 1]),
            sum(avdata[, 3] > avdata[, 1]),
            sum(avdata[, 4] > avdata[, 1])) # 
)
ratiodt <- data.frame(
  Column = c('NP',"MP",'LP','PP'),  # 列名称
  Number = table(Mac$group)# 
)
ratiodt$Column <- factor(ratiodt$Column,levels = c( 'NP',"MP",'LP','PP'))
ggplot(ratiodt, aes(x = Column, y = Number.Freq, fill = Column)) +
  geom_bar(stat = "identity", width = 0.6) +  # 绘制条形图
  scale_fill_manual(values = pal(4)) +  # 自定义颜色
  geom_text(aes(label = Number.Freq),  # 标记数字，保留两位小数
            vjust = -0.5, color = "black", size = 3)+
  labs(
    title = "",  # 图表标题
    x = "",                     # x轴标签
    y = "Number of Macropages",
  ) +
  theme_minimal() +  # 使用简洁主题
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5),  # 标题居中
        panel.grid.major = element_blank(),  # 去除主网格线
        panel.grid.minor = element_blank(),  # 去除次网格线
        axis.line = element_line(color = "black"),  # 添加 X 轴和 Y 轴线
        axis.ticks = element_line(color = "black")  # 添加刻度线
  ) +
  scale_x_discrete(expand = c(0, 0)) +  # x轴范围从0开始
  scale_y_continuous(expand = c(0, 0),limits = c(0, 3000))  # 标题居中

###############figR靶基因#########
library(ggraph)
library(tidygraph)
library(igraph)
library(FigR)
library(dplyr)
library(patchwork)
figR.d <-  readRDS("E:/2.工作/YF/mHeart/Mac1_figR.d.Rds")
rank <- rankDrivers(figR.d,rankBy = "meanScore")$data
rank <- rank %>%
  arrange(Score) %>% # 按 score 升序排序
  slice(c(1:4, (n() - 3):n())) # 筛选前 10 和后 10 行
figR.d<- figR.d[figR.d$Score > 0.5 | figR.d$Score < -0.5, ]
Deg.dt <- FindMarkers(obj ,group.by = 'group',ident.1 = 'PP',ident.2 = 'NP',logfc.threshold = 0,min.pct = 0)
Deg.dt$gene <- rownames(Deg.dt)
# 创建网络
tf_target_matrix <- data.frame(
  TF = figR.d$Motif,
  Target = figR.d$DORC,
  log2FC = Deg.dt[figR.d$DORC,2], # 表达量
  Ratio = Deg.dt[figR.d$DORC,3] # 百分率
)
tf_target_matrix <- na.omit(tf_target_matrix)
# 将矩阵转换为图对象
graph <- as_tbl_graph(tf_target_matrix, dirMacted = TRUE)
# 添加节点属性：表达量和百分率
specified_TFs <- rank$Motif
graph <- graph %>%
  activate(nodes) %>%
  mutate(
    Type = ifelse(name %in% tf_target_matrix$TF, "TF", "Target"),
    log2FC = ifelse(Type == "TF", 
                    tf_target_matrix$log2FC[match(name, tf_target_matrix$TF)], 
                    tf_target_matrix$log2FC[match(name, tf_target_matrix$Target)]),
    Ratio = ifelse(Type == "TF", 
                   tf_target_matrix$Ratio[match(name, tf_target_matrix$TF)], 
                   tf_target_matrix$Ratio[match(name, tf_target_matrix$Target)]),
    Size = ifelse(Type == "TF", Ratio, 0),
    Label = ifelse(name %in% specified_TFs, name, NA) # 只标记 TF
  )
# 绘制网络图
ggraph(graph, layout = "circle") + 
  geom_edge_link(
    arrow = arrow(length = unit(2, 'mm')), 
    color = "gray50", 
    alpha = 0.6
  ) +
  geom_node_point(
    aes(size = Size, color = log2FC), # 使用调整后的 Size
    alpha = 0.8
  ) +
  scale_color_gradient(low = "#FFFFFE", high = "#444576") + 
  scale_size(
    range = c(0, 10)
  ) +
  geom_text_repel(
    aes(x = x, y = y, label = Label), # 仅显示 TF 标签
    size = 4, 
    color = "black",
    box.padding = 0,  # 调整标签间距
    max.overlaps = Inf, # 允许无限避让
    min.segment.length = 0.1 # 调整标签连接线长度
  ) +
  theme_void() +
  labs(
    title = "Transcription Factor - Target Network",
    color = "log2FC"
  )

############相关性散点图#########
library(Hmisc)
library(Signac)
library(FigR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BiocParallel)
library(dplyr)
aUCell <- readRDS("E:/2.工作/YF/mHeart/aUCell.Rds")
mHeart_ATAC <- readRDS("E:/2.工作/YF/mHeart/mHeart_ATAC.Rds")
mHeart_ATAC[["aUCell"]] <- aUCell
obj[["aUCell"]] <-  mHeart_ATAC[,colnames(obj)][["aUCell"]]
obj[["motif"]] <- mHeart_ATAC[,colnames(obj)][["motif"]]
rm(aUCell);rm(mHeart_ATAC)
figR.d <-  readRDS("E:/2.工作/YF/mHeart/Mac1_figR.d.Rds")
DefaultAssay(obj) <- 'motif'
ranklist <- rankDrivers(figR.d,rankBy = "meanScore")$data
rownames(ranklist) <- ranklist$Motif
#'fatty acid beta oxidation' 'cardiac ventricle development'
dt <- data.frame(t(GetAssayData(obj[ranklist$Motif,],assay = 'motif',slot = 'data')),
                 cardiac = GetAssayData(obj,assay = 'aUCell',slot = 'data')['phagocytosis',])
res <- rcorr(as.matrix(dt),type	= 'pearson')
res_p <- res$P
plot.dt <- data.frame(score = ranklist[colnames(dt),]$Score,cor_r= res$r[,'cardiac'],genes =  colnames(dt) )
plot.dt$color <- ifelse(plot.dt$genes %in% 'Mafg','Mafg','Others')
ranklist <- as.data.frame(ranklist)
esrra_data <- subset(plot.dt, genes %in% c('Mafg'))
ggplot(data = plot.dt, mapping = aes(x = score, y = cor_r, color = color)) +
  geom_point() +
  xlim(-0.2, 0.2) +
  ylim(-0.6, 0.6) +
  theme_light() +
  geom_text_repel(
    data = esrra_data,  # 指定子集数据
    aes(label = genes),  # 标签内容
    color = "black",  # 标签颜色（与散点区分）
    nudge_x = 0.00,  # 水平偏移量
    nudge_y = 0.02,  # 垂直偏移量
    size = 4 )+
  scale_color_manual(values = c('#C65762','#444576')) +
  theme(
    panel.grid.major = element_blank(), # 移除主网格线
    panel.grid.minor = element_blank(), # 移除次网格线
    axis.line = element_line(color = "black"),
    legend.position = 'None') +
  xlab('Regulation score')+ylab('Pearson correlation coefficient \n (Motif activities vs Endothelial Cell Proliferation)')

############UMAP########
FeatureDimPlot(
  srt = subset(Mac),
  , features = c('innate immune response'), assay="aUCell",slot='data',bg_color = 'gray95',split.by = 'group',
  reduction = "umap", theme_use = "theme_blank",raster = F,ncol = 2)
FeatureStatPlot(subset(Mac),group.by='group',assay = 'aUCell',
                stat.by = c('phagocytosis','response to interleukin 12','immature t cell proliferation'),
                xlab = '',palcolor = pal(4),
                box_color = 'grey40',box_width = 0.05,
                legend.position = 'none',
                add_box = TRUE, stack = TRUE,ylab = '')
########GESAA#######
obj <- RunDEtest(srt = obj , group_by = "group", assay = 'activities',
                 fc.threshold = 1, only.pos = FALSE,min.pct = 0.25,
                 group1 = 'NP')
DEGs_group <- obj@tools$DEtest_custom$AllMarkers_wilcox
for (j in 2:4){
  obj<- RunDEtest(srt = obj, group_by = "group", assay = 'activities',
                  fc.threshold = 1, only.pos = FALSE,min.pct = 0.25,
                  group1 = samplegroup[j],group2 = 'NP')
  DEGs_group <- rbind(DEGs_group,obj@tools$DEtest_custom$AllMarkers_wilcox)
}
obj@tools$DEtest_group$AllMarkers_wilcox <- DEGs_group
obj <- RunGSEA(obj,
               group_by = "group", DE_threshold = "p_val_adj < 0.05",
               scoreType = "std", db = "GO_BP", species = "Mus_musculus"
)
gseadt <- obj@tools$GSEA_group_wilcox$enrichment
GSEAPlot(
  srt = obj, group_by = "group", group_use = 'PP', plot_type = "bar",
  direction = "both", topTerm = 20,pvalueCutoff = 0.5,padjustCutoff = NULL
)
GSEAPlot(srt = obj, group_by = "group", group_use = "LP", id_use = "GO:0080090")

##############motif差异分析#############
library(dplyr)
library(BiocParallel)
library(clusterProfiler)
library(org.Mm.eg.db)
setwd('E:/2.工作/YF/mHeart/8.终稿调整/20250315/差异motif')
register(BPPARAM = SnowParam(1))
samplegroup <-c('NP','MP','LP','PP')
obj <- subset(Mac,subtype %in% c('Macrophage-1'))
obj$sample <- obj$orig.ident
obj[['motif']] <- CreateAssayObject(mHeart_ATAC[,colnames(obj)]@assays$motif$data)
obj <- RunDEtest(srt = obj , group_by = "group", assay = 'motif',
                 fc.threshold = 1.5, only.pos = FALSE,min.pct = 0.25,
                 group1 = 'NP')
DEGs_group <- obj@tools$DEtest_custom$AllMarkers_wilcox
for (j in 2:4){
  obj<- RunDEtest(srt = obj, group_by = "group",  assay = 'motif',
                  fc.threshold = 1.5, only.pos = FALSE,min.pct = 0.25,
                  group1 = samplegroup[j],group2 = 'NP')
  DEGs_group <- rbind(DEGs_group,obj@tools$DEtest_custom$AllMarkers_wilcox)
}
obj@tools$DEtest_group$AllMarkers_wilcox <- DEGs_group
p_dt <- DEGs_group %>%
  filter(p_val_adj < 0.05) %>%
  filter(pct.1 > 0.25) %>%
  group_by(group1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 5)
GroupHeatmap(obj,
             features = unique( p_dt$gene),nlabel = 0,show_row_names = T,row_names_side = 'left',
             group.by = c("group"),assay = 'motif',group_palette = 'jama',height = 4,width = 1
)
writexl::write_xlsx(DEGs_group,'Mac_DEMotif.xlsx')
#gene
obj <- RunDEtest(srt = obj , group_by = "group", assay = 'RNA',
                 fc.threshold = 1.5, only.pos = FALSE,min.pct = 0.25,
                 group1 = 'NP')
DEGs_group <- obj@tools$DEtest_custom$AllMarkers_wilcox
for (j in 2:4){
  obj<- RunDEtest(srt = obj, group_by = "group",  assay = 'RNA',
                  fc.threshold = 1.5, only.pos = FALSE,min.pct = 0.25,
                  group1 = samplegroup[j],group2 = 'NP')
  DEGs_group <- rbind(DEGs_group,obj@tools$DEtest_custom$AllMarkers_wilcox)
}
obj@tools$DEtest_group$AllMarkers_wilcox <- DEGs_group
p_dt <- DEGs_group %>%
  filter(p_val_adj < 0.05) %>%
  filter(pct.1 > 0.25) %>%
  group_by(group1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 5)
GroupHeatmap(obj,
             features = unique( p_dt$gene),nlabel = 0,show_row_names = T,row_names_side = 'left',
             group.by = c("group"),assay = 'RNA',group_palette = 'jama',height = 4,width = 1
)

writexl::write_xlsx(DEGs_group,'Mac1_DERNA.xlsx')

#########rougue###
library(ROGUE)
library(tibble)
expr <- data.frame(GetAssayData(Mac,assay = "RNA",slot = "data"))
rogue.res <- rogue(expr, labels = Mac$subtype, samples = Mac$orig.ident, platform = "UMI",
                   span = 0.6)
p <- rogue.boxplot(rogue.res) 
df <- p[["data"]]
df$clusters <- gsub('Macrophage','Mac',df$clusters )
# 绘制箱型图
ggplot(na.omit(df), aes(x=clusters, y=ROGUE)) +
  geom_boxplot(aes(color=clusters)) +  
  geom_jitter(aes(color=clusters), width=0, size=1.5) +  
  scale_color_manual(values=pal(4)) + 
  theme_classic() + 
  theme(axis.text = element_text(color = 'black')) + 
  labs(title="", x="", y="ROGUE index") + 
  # scale_y_continuous(limits = c(0.75, 0.95), breaks = seq(0.75, 0.95, by = 0.1)) +  # 设置 y 轴范围和间隔
  theme(legend.position = "none") 
