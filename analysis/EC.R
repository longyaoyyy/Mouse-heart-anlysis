library(Seurat)
library(ggplot2)
library(ggunchained) 
library(reshape2) 
library(SCP)
set.seed(4180)
setwd("E:\\2.工作/YF/mHeart/")
#########color
cols <- c( "#444576", "#4682B4","#87CEEB","#FFD700", "#FFA500", "#C65762", "#4D4D4D", "#A9A9A9", "#FFFFEA")
pal <- colorRampPalette(cols)
EC <-  readRDS("E:/2.工作/YF/mHeart/EC.Rds")
obj <- subset(EC,subtype %in% c('EC-1'))
#########miloR######
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
DefaultAssay(EC) <- 'RNA'
EC@assays <- EC@assays[1]
traj_milo_sce <- as.SingleCellExperiment(EC)
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
sim <- list(Y=EC@assays$RNA$scale.data %>% as.data.frame(),
            meta.data=EC@meta.data %>% as.data.frame())
out <- scDist(normalized_counts = sim$Y, # 标准化的数据矩阵
              meta.data = sim$meta.data, # metadata表格
              d = 20, # 指定用于PCA分析的维度数量为前20
              fixed.effects = "group", # 你感兴趣的分组条件，对应metadata表格中的列名称
              random.effects = c('orig.ident'),  # 需要去除的潜在的影响因素，比如不同的样本、年龄、批次等等，对应metadata表格中的列名称
              clusters="subtype" # 待分析的细胞类型，对应metadata表格中的列名称
)
DistPlot(out)+ theme_bw()
p
#plot
distGenes(out, cluster = "EC-2") # 指定4这一群细胞
df <- data.frame(value = out$vals[["EC-2"]]$beta.hat, 
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
obj <- subset(EC,subtype %in% c('EC-1','EC-2'))
obj <- RunDEtest(srt = obj , group_by = "group", 
                   fc.threshold = 1, only.pos = FALSE,min.pct = 0.25,
                   group1 = 'NP')
DEGs_group <- obj@tools$DEtest_custom$AllMarkers_wilcox
 for (j in 2:4){
    obj<- RunDEtest(srt = obj, group_by = "group", 
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
    slice_head(n = 50)
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
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
res <- data_GO@compareClusterResult
library(dplyr)
for (i in 1:dim(res)[1]) {
  arr = unlist(strsplit(as.character(res[i,"geneID"]), split="/"))
  gene_names = paste(unique(data$SYMBOL[data$ENTREZID %in% arr]), collapse="/")
  res[i,"geneID"] = gene_names
}
writexl::write_xlsx(res,'EC1+2_四组GO通路富集结果.xlsx')
dotplot(data_GO, showCategory= 3,font.size = 10)

##############差异通路bar图#####4.71*4.57#####
library(BiocParallel)
library(dplyr)
library(patchwork)
library(clusterProfiler)
register(BPPARAM = SnowParam(1))
samplegroup <- c('NP','MP','LP','PP')
#EC
EC <- RunDEtest(EC, group_by = "subtype", only.pos = FALSE, fc.threshold = 1,min.pct = 0.1)
DEGs_group <- EC@tools$DEtest_subtype$AllMarkers_wilcox
#obj
obj <- RunDEtest(obj, group_by = "group", only.pos = FALSE, fc.threshold = 1)
obj <- RunDEtest(srt = obj , group_by = "group",
                 fc.threshold = 1, only.pos = FALSE,min.pct = 0.1,
                 group1 = 'NP')
DEGs_group <- obj@tools$DEtest_custom$AllMarkers_wilcox
for (j in 2:4){
  obj <- RunDEtest(srt = obj , group_by = "group", 
                   fc.threshold = 1, only.pos = FALSE,min.pct = 0.1,
                   group1 = samplegroup[j],group2 = 'NP')
  DEGs_group <- rbind(DEGs_group,obj@tools$DEtest_custom$AllMarkers_wilcox)
}
obj@tools$DEtest_group$AllMarkers_wilcox <- DEGs_group
#plot
allmarker <- DEGs_group %>%
  filter(p_val_adj < 0.05) %>%
  filter(pct.1 > 0.1) %>%
  filter(avg_log2FC > log2(1.5)) 
df_sig <- subset(allmarker)
group <- data.frame(gene=df_sig$gene,
                    group=df_sig$group1)
Gene_ID <- bitr(group$gene, fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Mm.eg.db")
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
levels(data$group) <-  c('vCM-1','vCM-2','vCM-3','vCM-4')
#levels(data$group) <- c('NP','MP','LP','PP')
data_GO <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Mm.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
res <- data_GO@compareClusterResult
for (i in 1:dim(res)[1]) {
  arr = unlist(strsplit(as.character(res[i,"geneID"]), split="/"))
  gene_names = paste(unique(data$SYMBOL[data$ENTREZID %in% arr]), collapse="/")
  res[i,"geneID"] = gene_names
}
enrich <- subset(res,Description %in% c('calcium ion homeostasis','vascular endothelial growth factor signaling pathway','regulation of endothelial cell chemotaxis',
                                        'wound healing','sprouting angiogenesis',
                                        'regulation of protein stability','ribosomal large subunit assembly','cellular response to oxidative stress',
                                        'actin filament organization','small GTPase-mediated signal transduction','extracellular matrix organization',
                                        'epithelial cell migration','regulation of chemotaxis'))
enrich <- subset(res,Description %in% c('regulation of glutamate receptor clustering','excitatory chemical synaptic transmission','cell-cell adhesion mediated by cadherin',
                                        'fatty acid metabolic process','carnitine metabolic process','skeletal muscle contraction',
                                        'regulation of epithelial cell migration','regulation of endothelial cell migration','regulation of vascular permeability',
                                        'myofibril assembly','actin filament organization','sarcomere organization'))
colnames(enrich ) <- gsub('Groups','Cluster',colnames(enrich))
dt <- enrich
dt$order <- gsub('NP','1',dt$Cluster);dt$order <- gsub('MP','2',dt$order);dt$order <- gsub('LP','3',dt$order);dt$order <- gsub('PP','4',dt$order)
dt$order <- gsub('EC-1','1',dt$Cluster);dt$order <- gsub('EC-2','2',dt$order);dt$order <- gsub('EC-3','3',dt$order);dt$order <- gsub('EC-4','4',dt$order);dt$order <- gsub('EC-5','5',dt$order)

dt <- dt[order(dt$pvalue), ]
dt <-  dt %>%
  group_by(group) %>%
  arrange(desc(pvalue)) %>%
  slice_tail(n = 3)
dt <- dt[order(dt$order,decreasing =F), ]
table(dt$Cluster)
cols <- RColorBrewer::brewer.pal(4,'Set1')
cols <- c('#4682B4','#C65762', '#7C9895','#DAA87C')
pal <-  colorRampPalette(cols)
#plot
dt$color <- factor(c(rep(pal(5)[3],each=3),rep(pal(5)[2],each=3),rep(pal(5)[1],each=3),rep(pal(5)[4],each=3),rep(pal(5)[5],each=3)),levels = pal(5))
dt$Description<- factor(dt$Description, levels = dt$Description)
dt$geneID <-  paste(substr(dt$geneID, start = 1, stop = 40),'...',sep = '')
dt$place <- factor(rownames(dt),levels = rownames(dt))
mytheme <- theme(
  axis.title = element_text(size = 15),
  axis.text = element_text(size = 15),
  axis.text.y = element_blank(), # 在自定义主题中去掉 y 轴通路标签:
  axis.ticks.length.y = unit(0,"cm"),
  plot.title = element_text(size = 15, hjust = 0.5),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
)
ggplot(data = dt, aes(x = -log10(pvalue), y = rev(place), fill = Cluster)) +
  scale_fill_manual(values =pal(5)) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.8,position = position_dodge()) +
  scale_x_continuous(expand = c(0,0)) + # 调整柱子底部与y轴紧贴 
  labs(x = "-Log10(pvalue)", y = "EC-5     EC-4      EC-3      EC-2     EC-1", 
       title = "EC GOBP pathway enrichment") + 
  # x = 0.61 用数值向量控制文本标签起始位置
  geom_text(size=4, aes(x = 0.05, label = Description), hjust = 0) + # hjust = 0,左对齐
  theme_classic() + 
  mytheme +
  NoLegend()
###########promoter ATAC######
library(ggrastr)
setwd("E:/2.工作/YF/mHeart/8.终稿调整/2050314/figS8_EC")
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
         title = 'EC')
  ggsave(paste('j_EC_NPvs',j,'所有细胞染色质开放片段数散点图.pdf',sep = ''),p,
         width = 3.84,height = 3.68)
}
ratiodt <- data.frame(
  Column = c("NP", "MP",'LP','PP'),  # 列名称
  Number = colSums(avdata < 10)  # 每列小于100的比例
)
ratiodt$Column <- factor(ratiodt$Column,levels = c("NP", "MP",'LP','PP'))
ggplot(ratiodt, aes(x = Column, y = Number, fill = Column)) +
  geom_bar(stat = "identity", width = 0.6) +  # 绘制条形图
  scale_fill_manual(values = pal(5)[1:4]) +  # 自定义颜色
  geom_text(aes(label = Number),  # 标记数字，保留两位小数
            vjust = -0.5, color = "black", size = 3)+
  labs(
    title = "",  # 图表标题
    x = "",                     # x轴标签
    y = "Number of inactive genes",
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
  scale_y_continuous(expand = c(0, 0),limits = c(0, 6000))  # 标题居中

###############figR靶基因#########
library(ggraph)
library(tidygraph)
library(igraph)
library(FigR)
figR.d <-  readRDS("E:/2.工作/YF/mHeart/EC1_figR.d.Rds")
rank <- rankDrivers(figR.d,rankBy = "meanScore")$data
rank <- rank %>%
  arrange(Score) %>% # 按 score 升序排序
  slice(c(1:5, (n() - 4):n())) # 筛选前 10 和后 10 行
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
graph <- as_tbl_graph(tf_target_matrix, directed = TRUE)
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
rm(rUCell);rm(mHeart_ATAC)
figR.d <-  readRDS("E:/2.工作/YF/mHeart/EC2_figR.d.Rds")
DefaultAssay(obj) <- 'motif'
ranklist <- rankDrivers(figR.d,rankBy = "meanScore")$data
rownames(ranklist) <- ranklist$Motif
#'fatty acid beta oxidation' 'sprouting angiogenesis''endothelial cell migration'
dt <- data.frame(t(GetAssayData(obj[ranklist$Motif,],assay = 'motif',slot = 'data')),
                 cardiac = GetAssayData(obj,assay = 'aUCell',slot = 'data')['sprouting angiogenesis',])
res <- rcorr(as.matrix(dt),type	= 'pearson')
res_p <- res$P
plot.dt <- data.frame(score = ranklist[colnames(dt),]$Score,cor_r= res$r[,'cardiac'],genes =  colnames(dt) )
plot.dt$color <- ifelse(plot.dt$genes %in% 'Gata3','Gata3','Others')
esrra_data <- subset(plot.dt, genes %in% c('Hif1a','Epas1','Gata2','Pknox1','Ppard'))

ggplot(data = plot.dt, mapping = aes(x = score, y = cor_r, color = color)) +
  geom_point(alpha=0.5) +
  xlim(-0.3, 0.3) +
  ylim(-0.45, 0.45) +
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
  xlab('Regulation score')+ylab('Pearson correlation coefficient \n (Motif activities vs sprouting angiogenesis)')

############UMAP########
mHeart[['rUCell']] <- rUCell
EC[['rUCell']] <- mHeart[,colnames(EC)][['rUCell']]
FeatureDimPlot(
  srt = subset(EC)
  , features = c('Vegfr2'), assay="RNA",slot = 'data',bg_color = 'gray95',split.by = 'group',
  reduction = "umap", theme_use = "theme_blank",raster = F,ncol = 2
)
FeatureStatPlot(subset(obj),group.by='group',assay = 'activities',
                stat.by = c('Nfyc'),
                xlab = '',palcolor = pal(4),
                box_color = 'grey40',box_width = 0.05,
                legend.position = 'none',
                add_box = TRUE, stack = TRUE,ylab = '')
library(BiocParallel)
library(dplyr)
register(BPPARAM = SnowParam(1))
allmarker <- FindAllMarkers(EC,min.pct = 0.25,logfc.threshold = 0.5,only.pos = T,group.by = 'subtype',assay = 'RNA',slot = 'data')
top_genes <-allmarker %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50)
p1=GroupHeatmap(
  srt = subset(EC), features = top_genes$gene, group.by = "subtype",show_row_names =F,,
  cluster_rows =FALSE, cluster_columns = FALSE,assay = 'RNA',slot = 'data',
  show_column_names = F,column_title = '',heatmap_palcolor = c('white','white','brown'),
  feature_split = top_genes$cluster,feature_split_palette = 'Safe',nlabel = 20,
  group_palette = 'Set2',height = 4,width = 1,anno_terms = T,species = "Mus_musculus",db='GO_BP'
)
########motif差异分析#######
library(dplyr)
library(BiocParallel)
library(clusterProfiler)
library(org.Mm.eg.db)
register(BPPARAM = SnowParam(1))
samplegroup <-c('NP','MP','LP','PP')
obj <- subset(EC,subtype %in% c('EC-1','EC-2'))
obj$sample <- obj$orig.ident
obj[['motif']] <- CreateAssayObject(obj@assays$motif$data)
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
ggsave(paste(i,samplegroup[j],'vs_NP_差异基因热图.pdf',sep = '_'),p1,
       width =3,height =5
)
writexl::write_xlsx(DEGs_group,'EC_DEMotif.xlsx')
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

writexl::write_xlsx(DEGs_group,'EC_DERNA.xlsx')

###########GSEA######
library(BiocParallel)
library(dplyr)
library(patchwork)
library(pheatmap)
library(ggrepel)
library(rtracklayer)
register(BPPARAM = SnowParam(1))
DefaultAssay(EC) <- "RNA"
EC <- RunDEtest(EC, group_by = "subtype", only.pos = FALSE, fc.threshold = 1)
EC <- RunGSEA(EC,
               group_by = "subtype", DE_threshold = "p_val_adj < 0.05",
               scoreType = "std", db = "GO_BP", species = "Mus_musculus"
)
gsea_data <- EC@tools$GSEA_subtype_wilcox$enrichment
writexl::write_xlsx(gsea_data,'EC亚群_gsea分析.xlsx')
GSEAPlot(
  srt = EC, group_by = "subtype", group_use = "EC-1", plot_type = "bar",
  direction = "both", topTerm = 20,pvalueCutoff = 0.5,padjustCutoff = NULL
)
#GO:0006635#GO:0006007#GO:0042773#GO:0009145#GO:0006099#GO:0042180#GO:0006089#GO:0006096
GSEAPlot(srt = vCM, group_by = "", group_use = "LP", id_use = "GO:0006635")
gsea_data <- vCM@tools$GSEA_subtype_wilcox$enrichment
gsea_data <- obj@tools$GSEA_group_wilcox$enrichment
writexl::write_xlsx(gsea_data,'vCM2_各时期_gsea分析.xlsx')
p.dt <- subset(gsea_data,Description %in% c( 'negative regulation of cell cycle G2/M phase transition','cell-cell adhesion mediated by cadherin','cellular extravasation',
                                             'mitochondrial ATP synthesis coupled electron transport','purine ribonucleotide biosynthetic process','fatty acid beta-oxidation',
                                             'regulation of cytokine-mediated signaling pathway','regulation of G protein-coupled receptor signaling pathway',
                                             'heart valve morphogenesis','myofibril assembly','extracellular matrix organization','actomyosin structure organization'))
p.dt$Description <- factor(p.dt$Description,levels =rev(c( 'negative regulation of cell cycle G2/M phase transition','cell-cell adhesion mediated by cadherin','cellular extravasation',
                                                           'mitochondrial ATP synthesis coupled electron transport','purine ribonucleotide biosynthetic process','fatty acid beta-oxidation',
                                                           'regulation of cytokine-mediated signaling pathway','regulation of G protein-coupled receptor signaling pathway',
                                                           'heart valve morphogenesis','myofibril assembly','extracellular matrix organization','actomyosin structure organization') ))
#气泡图
library(hrbrthemes)
library(viridis)
ggplot(p.dt,aes(x=Groups, y=Description, fill=NES)) +
  geom_tile(color = "black", linewidth = 0.3)+
  scale_size(range = c(.1, 8), name="-log10(pvalue)")+theme_light()+
  scale_fill_gradient2(low = '#444576',mid = 'white',high = '#FFA500',na.value = "white")+
  theme(legend.position="right",
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.title = element_text(size = 7)) +
  ylab("") +
  xlab("") 
library(ComplexHeatmap)

# 准备矩阵
heatmap_data <- reshape2::dcast(p.dt, Description ~ Groups, value.var = "NES")
rownames(heatmap_data) <- heatmap_data$Description
heatmap_data <- as.matrix(heatmap_data[, -1])

# 定义颜色映射
col_fun <- circlize::colorRamp2(
  breaks = c(min(heatmap_data, na.rm = TRUE), 0, max(heatmap_data, na.rm = TRUE)),
  colors = c('#444576', 'white', '#FFA500')
)

Heatmap(
  heatmap_data,
  col = col_fun,
  na_col = "white",          # 缺失值背景色设为白色
  rect_gp = gpar(col = "black"),  # 边框颜色
  row_names_side = "left",
  column_names_rot = 45,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  heatmap_legend_param = list(title = "NES"),
  width = unit(3, "cm"), 
  height = unit(9, "cm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (is.na(heatmap_data[i, j])) {
      grid.text("N/A", x, y, gp = gpar(fontsize = 8))  # 在NA单元格写"N/A"
    }
  }
)
#GO:0006635#GO:0009060#GO:0006099#GO:0006007#GO:0042773#GO:0009145#GO:0009260
GSEAPlot(srt = EC, group_by = "subtype", group_use = "EC-1", id_use = "GO:0001935")


#########rougue######
library(ROGUE)
expr <- data.frame(GetAssayData(EC,assay = "RNA",slot = "data"))
rogue.res <- rogue(expr, labels = EC$subtype, samples = EC$orig.ident, platform = "UMI",
                   span = 0.6)
p <- rogue.boxplot(rogue.res) 
df <- p[["data"]]
# 绘制箱型图
ggplot(df, aes(x=clusters, y=ROGUE)) +
  geom_boxplot(aes(color=clusters), fill=NA) +  
  geom_jitter(aes(color=clusters), width=0, size=1.5) +  
  scale_color_manual(values=pal(5)) + 
  theme_classic() + 
  theme(axis.text = element_text(color = 'black')) + 
  labs(title="", x="", y="ROGUE index") + 
  #scale_y_continuous(limits = c(0.75, 0.95), breaks = seq(0.75, 0.95, by = 0.1)) +  # 设置 y 轴范围和间隔
  theme(legend.position = "none") 

########################
library(ggrastr)
  avdata <- AggregateExpression(subset(obj)
                                ,group.by = 'group',assays = 'activities',slot = 'counts')
  avdata <- as.data.frame(avdata$activities)
  for(j in c('MP','LP','PP')){
    avdata_df <- data.frame(NP=avdata$NP,test=avdata[,j])
    p <- ggplot(avdata_df ) +
      geom_point_rast(aes(x = log10(NP), y = log10(test)),size = 0.001,alpha = 0.1,color="gray50")+
      # 添加 x = y 斜线
      geom_abline(intercept = 0, slope = 1, color = "#EF3945",alpha=0.7,lty=2) +
      theme_bw()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = 'none')+
      coord_fixed()+
      labs(x = "log10(gene activity) of NP",
           y = paste("log10(gene activity) of ",j,sep = ''),
           title = 'EC1&2')
    ggsave(paste('EC12_NPvs',j,'所有细胞染色质开放片段数散点图.pdf',sep = ''),p,
           width = 3.84,height = 3.68)
  }
  ratiodt <- data.frame( Column=c("MP",'LP','PP'),Number=colSums(avdata[,-1]  > avdata [, 1]))
  ratiodt$Column <- factor(ratiodt$Column,levels = c("MP",'LP','PP'))
  ggplot(ratiodt, aes(x = Column, y = Number, fill = Column)) +
    geom_bar(stat = "identity", width = 0.6) +  # 绘制条形图
    scale_fill_manual(values = pal(4)) +  # 自定义颜色
    geom_text(aes(label = Number),  # 标记数字，保留两位小数
              vjust = -0.5, color = "black", size = 3)+
    labs(
      title = "",  # 图表标题
      x = "",                     # x轴标签
      y = "Number of active genes",
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
    scale_y_continuous(expand = c(0, 0),limits = c(0, 25000))  # 标题居中

###########差异基因数量柱形图##########
  library(dplyr)
  library(BiocParallel)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  register(BPPARAM = SnowParam(1))
  samplegroup <-c('NP','MP','LP','PP')
  celltypes  <- levels(mHeart$celltype)
  setwd("E:/2.工作/YF/mHeart/8.终稿调整/202504/各细胞类型DEG数")
  #RNA
  DEG_nb <- data.frame()
  obj <- subset(EC,subtype %in% 'EC-1')
  DEGs_group <- data.frame()
  for (j in 2:4){
    obj<- RunDEtest(srt = obj, group_by = "group", assay = 'RNA',
                    fc.threshold = 1.5, only.pos = F,min.pct = 0.1,
                    group1 = samplegroup[j],group2 = 'NP')
    DEGs_group <- rbind(DEGs_group,obj@tools$DEtest_custom$AllMarkers_wilcox)
  }
  DEGs_group$col <- ifelse(DEGs_group$avg_log2FC>0,'up','down')
  DEG_dt <- DEGs_group %>%
    filter(p_val_adj < 0.05)%>%
    group_by(group1,col) %>%
    summarise(count = n())
  DEG_dt$count_adjusted <- ifelse(DEG_dt$col == "down", -DEG_dt$count, DEG_dt$count)
  DEG_dt$col <- factor(DEG_dt$col, levels = c("up", "down"))
  p <- ggplot(DEG_dt, aes(x = group1, y = count_adjusted, fill = col)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_manual(values = rev(pal(4)[1:2])) +  # 颜色需与 up/down 对应
    geom_text(aes(label = count),  # 显示原始值（非负）
              vjust = ifelse(DEG_dt$col == "up", -0.5, 1.2),  # 调整文字位置
              color = "black", size = 3) +
    labs(
      title = "",
      x = "",
      y = "Number of DEGs (.vs NP)",
      fill=''
    ) +
    theme_minimal() +
    theme(
      legend.position = 'right',
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(min(DEG_dt$count_adjusted ) - 300, max(DEG_dt$count) + 300),  # 对称范围
      labels = abs  # y轴标签显示绝对值
    )
  
  ggsave(paste('EC-1_组间DEG数.pdf',sep = ''),p,
         width = 3.3,height = 3)
  
  
  #atac
  setwd("E:/2.工作/YF/mHeart/8.终稿调整/202504/各细胞类型差异peak数")
  DEG_nb <- data.frame()
  obj <- subset(EC,subtype %in% 'EC-1')
  DEGs_group <- data.frame()
  for (j in 2:4){
    obj<- RunDEtest(srt = obj, group_by = "group", assay = 'ATAC',
                    fc.threshold = 1.5, only.pos = F,min.pct = 0.1,
                    group1 = samplegroup[j],group2 = 'NP')
    DEGs_group <- rbind(DEGs_group,obj@tools$DEtest_custom$AllMarkers_wilcox)
  }
  DEGs_group$col <- ifelse(DEGs_group$avg_log2FC>0,'up','down')
  DEG_dt <- DEGs_group %>%
    filter(p_val_adj < 0.05)%>%
    group_by(group1,col) %>%
    summarise(count = n())
  DEG_dt$count_adjusted <- ifelse(DEG_dt$col == "down", -DEG_dt$count, DEG_dt$count)
  DEG_dt$col <- factor(DEG_dt$col, levels = c("up", "down"))
  p <- ggplot(DEG_dt, aes(x = group1, y = count_adjusted, fill = col)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_manual(values = rev(pal(4)[1:2])) +  # 颜色需与 up/down 对应
    geom_text(aes(label = count),  # 显示原始值（非负）
              vjust = ifelse(DEG_dt$col == "up", -0.5, 1.2),  # 调整文字位置
              color = "black", size = 3) +
    labs(
      title = "",
      x = "",
      y = "Number of differential peaks (.vs NP)",
      fill=''
    ) +
    theme_minimal() +
    theme(
      legend.position = 'right',
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(min(DEG_dt$count_adjusted ) - 300, max(DEG_dt$count) + 300),  # 对称范围
      labels = abs  # y轴标签显示绝对值
    )
  
  ggsave(paste('EC-1_ATAC差异peak数.pdf',sep = ''),p,
         width = 3.3,height = 3)
  
  
########fig1e deg+atac num#####3*3######
  library(BiocParallel)
  library(dplyr)
  library(patchwork)
  library(pheatmap)
  library(ggrepel)
  register(BPPARAM = SnowParam(1))
  setwd("E:/2.工作/YF/mHeart/8.终稿调整/202504/RNA和ATAC相关性图")
  rna.dt <- AggregateExpression(EC, assays = 'RNA',group.by = 'subtype')[[1]]
  act.dt <- AggregateExpression(EC, assays = 'activities',group.by = 'subtype')[[1]]
  genes.dt <- intersect(row.names(rna.dt),row.names(act.dt))
  p.dt <- data.frame(rna = rna.dt[genes.dt,1],activities = act.dt[genes.dt,1])
  cor_value <- cor.test(log10(p.dt$rna + 1), log10(p.dt$activities + 1))
  r <- round(cor_value$estimate, 3)
  p_value <- scales::pvalue(cor_value$p.value)
  ggplot(data = p.dt, mapping = aes(x = log10(rna+1), y = log10(activities+1))) +
    geom_point(color='gray80', size=0.1) +
    geom_smooth(method = "lm", color = "#4682B4", formula = y ~ x) +  # 添加线性回归线
    theme_light() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = 'none'
    ) +
    xlab('Normalized RNA counts') +  # 修正坐标轴标签
    ylab('Normalized gene activities') +
    annotate("text", 
             x = Inf, y = 0,  # 将文本放在右上角
             label = paste0("r = ", r, "\n", "p = ", p_value), 
             hjust = 1.1, vjust = 0.1,  # 微调文本位置
             size = 4, color = "black")
  