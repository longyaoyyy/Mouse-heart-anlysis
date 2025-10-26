library(Seurat)
library(ggplot2)
library(ggunchained) 
library(reshape2) 
library(SCP)
library(Signac)
set.seed(4180)
setwd("E:\\2.工作/YF/mHeart/")
#########color
cols <- c( "#444576", "#4682B4", "#AEDEEE","#FFA500", "#FFD790","#C65762",'#FBDFDE', "#F6EFCF","#BCB99F")
pal <- colorRampPalette(cols)
vCM <-  readRDS("E:/2.工作/YF/mHeart/vCM.Rds")
obj <- subset(vCM,subtype %in% 'vCM-2')
###########UMAP######
library(mascarade)
library(tidyverse)
library(data.table)
#cellplot
df <- FetchData(object=vCM, vars=c("umap_1","umap_2","subtype"))
maskTable <- generateMask( dims=df[,1:2], cluster=df$subtype, minDensity = 0.5,smoothSigma = 0.05 )
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(2, "cm")
)

ggplot(df, aes(x=umap_1, y=umap_2)) + 
  geom_point(aes(color=subtype ),size = 0.05) + 
  geom_path(data=maskTable[], aes(group=group),linewidth=0.8,linetype = 2) +
  coord_fixed() + 
  scale_color_manual(values = pal(4))+
  guides(x = axis, y = axis,color = guide_legend(override.aes = list(size=5))) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(arrow = arrow(type = "closed",length = unit(0.1, "inches"))),
        axis.title = element_text(hjust = 0.05,size = 9))+
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL)
#四组 feature plot
#data
library(homologene)
mHeart[['rUCell']] <- rUCell
vCM[['rUCell']] <- mHeart[,colnames(vCM)][['rUCell']]
genes <- read_excel("8.终稿调整/20250315/figS7/基因总结(1).xlsx",sheet = "vCMs")
genes_sub <- subset(genes,文章 %in% 'Defining cardiac functional recovery in  end-stage heart failure at single-cell  esolution')
genelist <- human2mouse(genes_sub$基因)
plotgene <- genelist$mouseGene
FeatureDimPlot(
  srt = subset(vCM), features = c('Ppargc1a','Ppargc1b','Ncoa1','Ncoa2','Perm1','Ep300','Crebbp'),
  legend.position = 'right', assay="RNA",slot='data',bg_color ='gray95',split.by = 'group',pt.size = 0.5,
  reduction = "umap", theme_use = "theme_blank",raster = F,ncol =4,cells = sample(colnames(vCM),10000)
)
FeatureStatPlot(subset(vCM),stat.by ='Carns1',assay="activities",
                slot='data',group.by = 'group',palcolor = pal(4))
CellDimPlot(
  srt = vCM,group.by = 'subtype',split.by = 'group',legend.position = 'none',label_repel = T,label_point_size = 0.01,
  reduction = "umap", label = T,label_insitu = T,pt.size = 0.1,theme_use = 'theme_blank',
  label.fg = 'black',label.bg = 'grey95',label.bg.r = 0.1,title = '',palcolor =pal(4) ,ncol = 4,
  bg_color = 'grey90',raster = F
)

GroupHeatmap(subset(vCM,subtype %in% 'vCM-2'),assay = 'RNA',slot = 'data',
                    features = c('Ppargc1a','Ppargc1b','Ncoa1','Ncoa2','Perm1','Ep300','Crebbp'
                    ),group_palette = 'npg',
                    group.by = c("group"), show_column_names = F
)
CellDensityPlot(vCM, features = c('Ppargc1a','Ppargc1b','Ncoa1','Ncoa2','Perm1','Ep300','Crebbp'), ncol = 7,
                group.by = "group",palcolor = rep(pal(4)),y.min=0.01,legend.position = 'none',combine = T)+ylab('')
########all GRN#########
library(ggraph)
library(tidygraph)
library(igraph)
library(FigR)
library(dplyr)
figR.d <-  readRDS("E:/2.工作/YF/mHeart/vCM2_PPfigR.d.Rds")
rank.p <- rankDrivers(figR.d,rankBy = "meanScore")
rank <- rank.p$data
plot.TF <- c('Esrra','Esrrg','Ppara','Ppard','Stat5b','Gata6','Nfyc','Nfya')
plot.TF <- c('Twist1','Prrx2','Foxo4','Foxp1','Foxj2','Klf2','Klf4','Klf6','Klf10','Klf15')
rank$TF <- ifelse(as.character(rank$Motif) %in% plot.TF ,as.character(rank$Motif),'') 
rank.p$data$TF<- ''
rank.p +ggrepel::geom_text_repel(aes(label=TF),rank,max.overlaps = 200,box.padding = 0.6,  min.segment.length = 0.5,ylim = c(0.001, NA),
                                 force = T)
#plot.TF <- 
plot.TF <- c('Esrrg','Esrra','Nfyb','Nfyc','Tead1','Gata6')
plot.TF <- c('Ppara','Pparg','Esrra','Esrrg','Trp53','Lef1','Stat5b')
plot.TF <- c('Esrrb','Ppard','Esrrg','Esrra','Nr2c2','Prdm4','Esrra')
rank <- subset(rank,Motif %in% plot.TF )
Deg.dt <- FindMarkers(obj ,group.by = 'group',ident.1 = 'PP',ident.2 = 'NP',logfc.threshold = 0,min.pct = 0)
Deg.dt$gene <- rownames(Deg.dt)
# 创建网络
figR.d<- figR.d[figR.d$Score > 0.5 | figR.d$Score < -0.5, ]
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
    color = "gray90", 
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
    box.padding = 0.1,  # 调整标签间距
    max.overlaps = Inf, # 允许无限避让
    min.segment.length = 0.1 # 调整标签连接线长度
  ) +
  theme_void() +
  labs(
    title = "Transcription Factor - Target Network",
    color = "log2FC",
    size = 'cell ratio'
  )

########correlation score##########
library(Hmisc)
library(Signac)
library(FigR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BiocParallel)
aUCell <- readRDS("E:/2.工作/YF/mHeart/aUCell.Rds")
mHeart_ATAC[["aUCell"]] <- aUCell
obj[["aUCell"]] <-  mHeart_ATAC[,colnames(obj)][["aUCell"]]
obj[["motif"]] <- mHeart_ATAC[,colnames(obj)][["motif"]]
rm(aUCell);rm(mHeart_ATAC)
figR.d <-  readRDS("E:/2.工作/YF/mHeart/vCM2_figR.d.Rds")
DefaultAssay(obj) <- 'motif'
ranklist <- rankDrivers(figR.d,rankBy = "meanScore")$data
rownames(ranklist) <- ranklist$Motif
#'fatty acid beta oxidation' 'cardiac ventricle development'，’deoxyribonucleoside triphosphate biosynthetic process
dt <- data.frame(t(GetAssayData(obj[ranklist$Motif,],assay = 'motif',slot = 'data')),
                 cardiac = GetAssayData(obj,assay = 'aUCell',slot = 'data')['tricarboxylic acid cycle',])
res <- rcorr(as.matrix(dt),type	= 'pearson')
res_p <- res$P
plot.dt <- data.frame(score = ranklist[colnames(dt),]$Score,cor_r= res$r[,'cardiac'],genes =  colnames(dt) )
plot.dt$color <- ifelse(plot.dt$genes %in% c('Esrra','Esrrg','Esrrb'),'ERR','Others')
esrra_data <- subset(plot.dt, genes %in% c('Esrra','Esrrg','Esrrb','Ppara','Ppard','Stat5b','Gata6','Nfyc','Nfya'))
ggplot(data = plot.dt, mapping = aes(x = score, y = cor_r, color = color)) +
  geom_point(aes(alpha = 0.5)) +
  xlim(-0.25, 0.25) +
  ylim(-0.2, 0.2) +
  theme_light() +
  geom_text_repel(
    data = esrra_data,  # 指定子集数据
    aes(label = genes),  # 标签内容
    color = "black",  # 标签颜色（与散点区分）
    nudge_x = 0.02,  # 水平偏移量
    nudge_y = 0.02,  # 垂直偏移量
    size = 4,)+
  scale_color_manual(values = c('#C65762','#444576')) +
  theme(
    panel.grid.major = element_blank(), # 移除主网格线
    panel.grid.minor = element_blank(), # 移除次网格线
    axis.line = element_line(color = "black"),
    legend.position = 'None') +
  xlab('Regulation score')+ylab('Pearson correlation coefficient \n (Motif activities vs tricarboxylic acid cycle)')
###########MiloR#####
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(Seurat)
setwd('E:\\2.工作/YF/mHeart')
DefaultAssay(vCM) <- 'RNA'
vCM@assays <- vCM@assays[1]
traj_milo_sce <- as.SingleCellExperiment(vCM)
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

###########scDist####
library(dplyr)
library(Seurat) # v4版本
library(scDist)
library(ggplot2)
library(ggpubr)
##run Dist
sim <- list(Y=vCM@assays$RNA$scale.data %>% as.data.frame(),
            meta.data=vCM@meta.data %>% as.data.frame())
out <- scDist(normalized_counts = sim$Y, # 标准化的数据矩阵
              meta.data = sim$meta.data, # metadata表格
              d = 20, # 指定用于PCA分析的维度数量为前20
              fixed.effects = "group", # 你感兴趣的分组条件，对应metadata表格中的列名称
              random.effects = c('orig.ident'),  # 需要去除的潜在的影响因素，比如不同的样本、年龄、批次等等，对应metadata表格中的列名称
              clusters="subtype" # 待分析的细胞类型，对应metadata表格中的列名称
)
p <- DistPlot(out, return.plot = TRUE)
p + theme_bw()
#plot
distGenes(out, cluster = "vCM-2") # 指定4这一群细胞
df <- data.frame(value = out$vals[["vCM-2"]]$beta.hat, 
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
##########ATAC scatter######
library(ggrastr)
setwd("E:/2.工作/YF/mHeart/8.终稿调整/20250331/ATAC")
avdata <- AggregateExpression(obj,group.by = 'group',assays = 'activities',slot = 'counts')
avdata <- as.data.frame(avdata$activities)
  for(j in c('MP','LP','PP')){
    avdata_df <- data.frame(NP=avdata$NP,test=avdata[,j])
    p <-  ggplot(avdata_df ) +
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
           title = 'vCM-2')
    ggsave(paste('vCM-2_NPvs',j,'所有细胞染色质开放片段数散点图.pdf',sep = ''),p,
           width = 3.84,height = 3.68)
  }
ratiodt <- data.frame(
  Column = c( "MP",'LP','PP'),  # 列名称
  Number =c(sum(avdata[, 2] > avdata[, 1]),
            sum(avdata[, 3] > avdata[, 1]),
            sum(avdata[, 4] > avdata[, 1])) # 
)
ratiodt$Column <- factor(ratiodt$Column,levels = c( "MP",'LP','PP'))
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

############vCM pathway##########
library(dplyr)
library(BiocParallel)
library(clusterProfiler)
register(BPPARAM = SnowParam(1))
samplegroup <-c('NP','MP','LP','PP')
for(i in c('vCM-1','vCM-3','vCM-4')){
  obj <- subset(vCM,subtype %in% 'vCM-2')
  obj<- RunDEtest(obj , group_by = "group", only.pos = FALSE, fc.threshold = 1)
  obj <- RunDEtest(srt = obj , group_by = "group",
                   fc.threshold = 1, only.pos = FALSE,min.pct = 0.25,
                   group1 = 'NP')
  DEGs_group <- obj@tools$DEtest_custom$AllMarkers_wilcox
  for (j in 2:4){
    obj <- RunDEtest(srt = obj , group_by = "group", 
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
    slice_head(n = 100)
  df_sig <- top_genes
  group <- data.frame(gene=df_sig$gene,
                      group=df_sig$group1)
  Gene_ID <- bitr(group$gene, fromType="SYMBOL", 
                  toType="ENTREZID", 
                  OrgDb="org.Mm.eg.db")
  data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
  levels(data$group) <- 
    levels(data$group) <- c('NP','MP','LP','PP')
  data_GO <- compareCluster(
    ENTREZID~group, 
    data=data, 
    fun="enrichGO", 
    OrgDb="org.Mm.eg.db",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    
  )
  res <- data_GO@compareClusterResult
  library(dplyr)
  for (i in 1:dim(res)[1]) {
    arr = unlist(strsplit(as.character(res[i,"geneID"]), split="/"))
    gene_names = paste(unique(data$SYMBOL[data$ENTREZID %in% arr]), collapse="/")
    res[i,"geneID"] = gene_names
  }
  writexl::write_xlsx(res,'vCM2四组GO通路富集结果.xlsx')
  data_GO@compareClusterResult <- subset(data_GO@compareClusterResult,Description %in% c(unique(pathway_dt$Description),"muscle hypertrophy"))
  dotplot(data_GO, showCategory= 5,font.size = 10)
}

##############DEG pathway bar plot#####4.71*4.57#####
library(BiocParallel)
library(dplyr)
library(patchwork)
library(clusterProfiler)
register(BPPARAM = SnowParam(1))
samplegroup <- c('NP','MP','LP','PP')
#vCM
vCM <- RunDEtest(vCM, group_by = "subtype", only.pos = FALSE, fc.threshold = 1,min.pct = 0.1)
DEGs_group <- vCM@tools$DEtest_subtype$AllMarkers_wilcox
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
enrich <- subset(res,Description %in% c('regulation of cellular component size','wound healing','regulation of blood circulation',
                                        'circadian rhythm','muscle adaptation','autophagy of mitochondrion',
                                        'cardiac muscle contraction','cardiac muscle tissue development','cardiac muscle cell differentiation',
                                        'fatty acid beta-oxidation','acyl-CoA metabolic process','lipid catabolic process'))
enrich <- subset(res,Description %in% c('generation of precursor metabolites and energy','autophagy of mitochondrion','muscle adaptation',
                                        'cardiac muscle tissue development','cardiac muscle cell differentiation','cardiac ventricle development',
                                        'lipid catabolic process','heart morphogenesis','muscle cell development',
                                        'regulation of blood circulation','regulation of cellular component size','actin filament organization'))
enrich <- subset(res,Description %in% c('generation of precursor metabolites and energy','autophagy of mitochondrion','muscle adaptation')& group %in% 'NP')
enrich <- rbind(enrich , subset(res,Description %in% c( 'cardiac muscle tissue development','cardiac muscle cell differentiation','cardiac ventricle development')& group %in% 'MP'))
enrich <- rbind(enrich , subset(res,Description %in% c( 'lipid catabolic process','heart morphogenesis','muscle cell development')& group %in% 'LP'))
enrich <- rbind(enrich , subset(res,Description %in% c( 'regulation of blood circulation','regulation of cellular component size','actin filament organization')& group %in% 'PP'))
colnames(enrich ) <- gsub('Groups','Cluster',colnames(enrich))
dt <- enrich
dt$order <- gsub('NP','1',dt$Cluster);dt$order <- gsub('MP','2',dt$order);dt$order <- gsub('LP','3',dt$order);dt$order <- gsub('PP','4',dt$order)
dt$order <- gsub('vCM-1','1',dt$Cluster);dt$order <- gsub('vCM-2','2',dt$order);dt$order <- gsub('vCM-3','3',dt$order);dt$order <- gsub('vCM-4','4',dt$order)

dt <- dt[order(dt$pvalue), ]

dt <- dt[order(dt$order,decreasing =F), ]
table(dt$Cluster)
cols <- RColorBrewer::brewer.pal(4,'Set1')
cols <- c('#4682B4','#C65762', '#7C9895','#DAA87C')
pal <-  colorRampPalette(cols)
#plot
dt$color <- factor(c(rep(pal(4)[3],each=3),rep(pal(4)[2],each=3),rep(pal(4)[1],each=3),rep(pal(4)[4],each=3)),levels = pal(4))
dt$Description<- factor(dt$Description, levels = dt$Description)
dt$geneID <-  paste(substr(dt$geneID, start = 1, stop = 40),'...',sep = '')
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
ggplot(data = dt, aes(x = -log10(pvalue), y = rev(Description), fill = Cluster)) +
  scale_fill_manual(values =pal(4)) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0)) + # 调整柱子底部与y轴紧贴 
  labs(x = "-Log10(pvalue)", y = "PP          LP           MP          NP", 
       title = "vCM2 GOBP pathway enrichment") + 
  # x = 0.61 用数值向量控制文本标签起始位置
  geom_text(size=5, aes(x = 0.05, label = Description), hjust = 0) + # hjust = 0,左对齐
  theme_classic() + 
  mytheme +
  NoLegend()
#############satt plot########
cofactor <- c('Perm1','Scarb1','Scarf2','Crebbp','Prkn','Kank2','Lancl1','Ppp2r2a')
mHeart[["aUCell"]] <- aUCell
vCM[['aUCell']] <- mHeart[,colnames(vCM)][["aUCell"]]
FeatureDimPlot(
  srt = subset(vCM),
  , features = c('Ncoa1','Ncoa2'), assay="RNA",slot='data',bg_color = 'gray95',split.by = 'group',
  reduction = "umap", theme_use = "theme_blank",raster = F,ncol =4,cells = sample(colnames(vCM),10000))
FeatureStatPlot(subset(vCM),group.by='group',assay = 'RNA',
                stat.by = c('Atp1b1'),slot = 'data',
                xlab = '',palcolor = pal(4),
                box_color = 'grey40',box_width = 0.05,
                legend.position = 'right',comparisons = list(c('vCM-1','vCM-2'),
                                                             c('vCM-3','vCM-2'),
                                                             c('vCM-4','vCM-2')),
                add_box = F, stack = TRUE,ylab = '')
FeatureDimPlot(vCM,add_density = TRUE,assay = 'activities',slot = 'data',
               features = "Ckap5", split.by = "group", reduction = "umap",
               cells.highlight = TRUE, theme_use = "theme_blank"
)
DotPlot(vCM,features = 'Ckap5',group.by = 'subtype',split.by = 'group',cols = pal(4))
allmarkers <- FindAllMarkers(vCM,group.by = 'subtype',assay = 'rUCell',logfc.threshold = 0.5,only.pos = T,min.pct = 0.2)
ave <- data.frame(AverageExpression(vCM,assays = 'rUCell',group.by = 'subtype')$rUCell)
avexp <- as.matrix(AverageExpression(vCM,assays = "RNA",slot = "data",group.by = 'subtype')$RNA[m.cc.genes$g2m.genes,])
pheatmap(subset(avexp), 
               color = colorRampPalette(c("#92a5d1",'#c5dff4',"white",'#FBDFE2', "#E83945"))(100),  # 设置颜色渐变
               scale = "row",  # 按行进行缩放
               show_rownames = TRUE,  # 显示行名
               show_colnames = TRUE,cluster_rows = T,cluster_cols = F,main = 'G2M pahse genes',cellwidth = 15,  # 设置每个方格的宽度
               cellheight = 15 )

vCM$allgroup <- factor(paste(vCM$subtype,vCM$group,sep = '_'),levels = c('vCM-1_NP','vCM-1_MP','vCM-1_LP','vCM-1_PP',
                                                                         'vCM-2_NP','vCM-2_MP','vCM-2_LP','vCM-2_PP',
                                                                         'vCM-3_NP','vCM-3_MP','vCM-3_LP','vCM-3_PP',
                                                                         'vCM-4_NP','vCM-4_MP','vCM-4_LP','vCM-4_PP'))
CellDensityPlot(vCM, features = "Ckap5", group.by = "allgroup",palcolor = rep(pal(4),each=4),y.min=0.01,legend.position = 'none')+ylab('')
###############monocle#########
library(monocle)
#创建CDS
vCM <- vCM[,sample(colnames(vCM),size = 10000)]
data <- vCM@assays$RNA$counts
cell_metadata <- subset(vCM,celltype %in% c('vCM'))@meta.data
gene_annotation <-data.frame(gene_short_name = rownames(vCM),row.names = row.names(vCM))
pd <- new('AnnotatedDataFrame', data = cell_metadata) 
fd <- new('AnnotatedDataFrame', data = gene_annotation)
cds <- newCellDataSet(data,phenoData =pd, featureData =fd,
                      lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())
rm(data)
rm(vCM)
#估计size factor和离散度
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
#选择基因
library(msigdbr)
m_df = msigdbr( species = "Mus musculus",  category = "C5",subcategory = 'GO:BP')
m_df_sub <- subset(m_df,gs_name %in% 'GOBP_MUSCLE_HYPERTROPHY')
ordergene <- m_df_sub$gene_symbol
cds <- setOrderingFilter(cds,ordergene)
#降维
cds <- reduceDimension(cds,max_components = 2,method = 'DDRTree')
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 7)#把State7设成拟时间轴的起始点
plot_cell_trajectory(cds,color_by="Pseudotime",size=3,show_branch_points = F)+ 
  facet_wrap("~group", nrow = 1)+
  scale_color_gradientn(colours = pal(2))
plot_cell_trajectory(cds,color_by="State",size=1,show_backbone=TRUE,show_branch_points = F)
plot_cell_trajectory(cds,color_by="subtype",size=1,show_backbone=TRUE,show_branch_points = F)+
  scale_color_manual(values = pal(4))+facet_wrap("~group", nrow = 1)+labs(color= 'vCM subtypes')
#寻找与时间相关基因
Time_diff <-differentialGeneTest(cds[ordergene ,], cores = 4,
                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_genes <- Time_diff %>%pull(gene_short_name) %>% as.character()
plotdata <- plot_pseudotime_heatmap(cds[Time_genes,],num_clusters=6, show_rownames=F, return_heatmap=T,cores=6)
plotdf=pData(cds)
library(ggridges)
library(tidyverse)
ggplot(plotdf, aes(x=Pseudotime,y=cellname,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  scale_fill_gradientn(name="Pseudotime",colors =colorRampPalette(pal(3))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )
clusters <- cutree(plotdata$tree_row, k = 6)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
#折线图
plot_genes_in_pseudotime(cds[c('EPAS1','HIF1A'),], color_by = "cellname")
pData(cds)$EPAS1 =log2(exprs(cds)['EPAS1',]+1)
pData(cds)$VWF =log2(exprs(cds)['VWF',]+1)
plot_cell_trajectory(cds, color_by =c("VWF"),
                     show_branch_points = F)+scale_color_gradientn(colours = pal(2))

########vCM place#############
zone_marker <- data.frame(read_excel("E:/2.工作/YF/mHeart/8.终稿调整/20250328/参考文献/zone_marker.xlsx"))
rownames(zone_marker) <- zone_marker[,1]
genelist <- subset(zone_marker,zone %in% 'RV')[,1]
FeatureDimPlot(
  srt = subset(vCM), features = c("Nppb"),legend.position = 'right', assay="RNA",slot='data',bg_color ='gray95',split.by = 'group',
  reduction = "umap", theme_use = "theme_blank",raster = F,ncol =4#,cells = sample(colnames(vCM),10000)
)

###########GSEA######
library(BiocParallel)
library(dplyr)
library(patchwork)
library(pheatmap)
library(ggrepel)
library(rtracklayer)
register(BPPARAM = SnowParam(1))
DefaultAssay(vCM) <- "RNA"
vCM <- RunDEtest(vCM, group_by = "subtype", only.pos = FALSE, fc.threshold = 1)
vCM <- RunGSEA(vCM,
group_by = "subtype", DE_threshold = "p_val_adj < 0.05",
scoreType = "std", db = "GO_BP", species = "Mus_musculus"
)
GSEAPlot(
srt = vCM, group_by = "subtype", group_use = "vCM-2", plot_type = "bar",
direction = "both", topTerm = 20,pvalueCutoff = 0.5,padjustCutoff = NULL
)
#obj
samplegroup <-c('NP','MP','LP','PP')
obj <- RunDEtest(obj, group_by = "group", only.pos = FALSE, fc.threshold = 1)
obj <- RunDEtest(srt = obj , group_by = "group",
                 fc.threshold = 1, only.pos = FALSE,
                 group1 = 'NP')
DEGs_group <- obj@tools$DEtest_custom$AllMarkers_wilcox
for (j in 2:4){
  obj <- RunDEtest(srt = obj , group_by = "group", 
                   fc.threshold = 1, only.pos = FALSE,
                   group1 = samplegroup[j],group2 = 'NP')
  DEGs_group <- rbind(DEGs_group,obj@tools$DEtest_custom$AllMarkers_wilcox)
}
obj@tools$DEtest_group$AllMarkers_wilcox <- DEGs_group
obj <- RunGSEA(obj,
               group_by = "group", DE_threshold = "p_val_adj < 0.05",
               scoreType = "std", db = "GO_BP", species = "Mus_musculus"
)
GSEAPlot(
  srt = obj, group_by = "group", group_use = "PP", plot_type = "bar",
  direction = "both", topTerm = 20,pvalueCutoff = 0.5,padjustCutoff = NULL
)
#GO:0006635#GO:0006007#GO:0042773#GO:0009145#GO:0006099#GO:0042180#GO:0006089#GO:0006096
GSEAPlot(srt = vCM, group_by = "subtype", group_use = "vCM-2", id_use = "GO:0006119",pvalueCutoff=0.05,padjustCutoff = NULL,lab)
gsea_data <- vCM@tools$GSEA_subtype_wilcox$enrichment
gsea_data <- obj@tools$GSEA_group_wilcox$enrichment
writexl::write_xlsx(gsea_data,'vCM2_各时期_gsea分析.xlsx')
p.dt <- subset(gsea_data,Description %in% c( 'negative regulation of cell cycle G2/M phase transition','cell-cell adhesion mediated by cadherin','cellular extravasation',
                                             'mitochondrial ATP synthesis coupled electron transport','positive regulation of chromosome segregation','regulation of oxidative phosphorylation',
                                           'regulation of cytokine-mediated signaling pathway','regulation of cytokine-mediated signaling pathway',
                                        'cardiac muscle contraction','myofibril assembly','extracellular matrix organization','actomyosin structure organization'))
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
  xlab("") +coord_flip()
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
  t(heatmap_data),
  col = col_fun,
  na_col = "white",          # 缺失值背景色设为白色
  rect_gp = gpar(col = "black"),  # 边框颜色
  row_names_side = "left",
  column_names_rot = 45,
  cluster_rows = FALSE,
  cluster_columns = T,column_dend_side = 'bottom',column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = "NES"),column_names_side = 'top',
  width = unit(9, "cm"), 
  height = unit(3, "cm"),
)
#########rougue######
library(ROGUE)
expr <- data.frame(GetAssayData(vCM,assay = "RNA",slot = "data"))
rogue.res <- rogue(expr, labels = vCM$subtype, samples = vCM$orig.ident, platform = "UMI",
                   span = 0.6)
p <- rogue.boxplot(rogue.res) 
df <- p[["data"]]
# 绘制箱型图
ggplot(df, aes(x=clusters, y=ROGUE)) +
  geom_boxplot(aes(color=clusters), fill=NA) +  
  geom_jitter(aes(color=clusters), width=0, size=1.5) +  
  scale_color_manual(values=pal(4)) + 
  theme_classic() + 
  theme(axis.text = element_text(color = 'black')) + 
  labs(title="", x="", y="ROGUE index") + 
  #scale_y_continuous(limits = c(0.75, 0.95), breaks = seq(0.75, 0.95, by = 0.1)) +  # 设置 y 轴范围和间隔
  theme(legend.position = "none") 

#########footprint########
motifs <- c('Esrra')
vCM <- Footprint(
  object = vCM,
  motif.name = motifs ,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  assay = 'ATAC'
)
library(ggplot2)

PlotFootprint(
  object = vCM,
  features = 'Esrra',
  group.by = 'group',
  assay = 'ATAC',label.top = 0
)+ scale_fill_brewer(type = "qual", palette = 1)

############G2M###########
m.cc.genes <- readRDS("database/mouse_cell_cycle_genes.rds") 
vCM <- CellCycleScoring(vCM, s.features = m.cc.genes$s.genes, 
                           g2m.features = m.cc.genes$g2m.genes,
                           set.ident = TRUE)
CellDimPlot(
  srt = vCM,group.by = 'Phase',legend.position = 'right',label_repel = T,label_point_size = 0.01,
  reduction = "umap", label = F,label_insitu = T,pt.size = 0.1,theme_use = 'theme_blank',
  label.fg = 'black',label.bg = 'grey95',label.bg.r = 0.1,title = '',palcolor =pal(3) ,ncol = 4,
  bg_color = 'grey90',raster = F
)
FeatureDimPlot(
  srt = subset(vCM),
  , features = c('S.Score'), assay="RNA",slot='data',bg_color = 'gray95',split.by = 'group',
  reduction = "umap", theme_use = "theme_blank",raster = F,ncol = 2)
CellStatPlot(vCM, stat.by = "Phase", group.by = "subtype", label = F,palcolor = pal(4),flip = F,plot_type = "trend",xlab = '',ylab = 'Percentage of cell types')
FeatureStatPlot(subset(vCM),group.by='subtype',assay = 'RNA',split.by = 'group',
                stat.by = c('G2M.Score','S.Score'),slot = 'data',
                xlab = '',palcolor = pal(4),
                box_color = 'grey40',box_width = 0.05,
                legend.position = 'right',plot_type = 'box',
                add_box = TRUE, stack = TRUE,ylab = '')
########Epitrace###############
library(Seurat)
library(SeuratObject)
library(Signac)
library(EpiTrace)
initiated_peaks <- Init_Peakset(GetAssayData(vCM,assay = "ATAC",layer = "ranges"))
initiated_peaks_df <- as.data.frame(initiated_peaks,row.names = NULL)
rownames(vCM@meta.data) -> cellname_vec
paste0(initiated_peaks_df$seqnames,'_',initiated_peaks_df$start,'_',initiated_peaks_df$end) -> initiated_peaks_df$peakName
initiated_mm <- Init_Matrix(cellname = cellname_vec,peakname = initiated_peaks_df$peakName,matrix = vCM@assays$ATAC@counts)
epitrace_obj_age_estimated <- EpiTraceAge_Convergence(initiated_peaks,initiated_mm,celltype = NULL,qualnum = 10,Z_cutoff = 2.5,mean_error_limit = 0.01,iterative_time = 20,parallel = F,ncore_lim = 46,ref_genome = 'mm10',non_standard_clock = F)
mtx@colData[epitrace_obj_age_estimated@meta.data$cell,]$seurat_clusters -> epitrace_obj_age_estimated@meta.data$archR_cluster
epitrace_obj_age_estimated@meta.data$archR_cluster <- factor(epitrace_obj_age_estimated@meta.data$archR_cluster,levels=paste0('C',c(1:10)))
epitrace_obj_age_estimated@meta.data %>% as.data.frame() -> epitrace_obj_age_estimated_meta
vCM$EpiTrace_Age <- epitrace_obj_age_estimated_meta$EpiTraceAge_iterative

#########Cytotrace#############
library(CytoTRACE2)
library(patchwork)
vCM <- cytotrace2(vCM,
                     is_seurat = TRUE,
                     slot_type = "counts",
                     species = "mouse",
                     seed = 1234)
annotation <- data.frame(phenotype=vCM@meta.data$group) %>% set_rownames(., colnames(vCM))
plots <- plotData(cytotrace2_result=vCM, annotation=annotation, is_seurat=TRUE)
p1 <- plots$CytoTRACE2_UMAP
p2 <- plots$CytoTRACE2_Potency_UMAP
p3 <- plots$CytoTRACE2_Relative_UMAP
p4 <- plots$CytoTRACE2_Boxplot_byPheno
(p1+p2+p3+p4) + plot_layout(ncol = 2)

##########Slingshot#############
library(BiocParallel)
register(BPPARAM = SnowParam(1))
vCM <- RunSlingshot(srt = vCM, group.by = "subtype", reduction = "umap",start = 'vCM-1')
CellDimPlot(vCM, group.by = "subtype", reduction = "UMAP", lineages = paste0("Lineage", 1:3), lineages_span = 0.1,palcolor = pal(4),
            label.fg = 'black',label.bg = 'grey95',label.bg.r = 0.1,title = '',ncol = 4,label = T,label_insitu = T,theme_use = "theme_blank",
            bg_color = 'grey90',raster = F)
FeatureDimPlot(vCM, features = paste0("Lineage", 1:3), reduction = "umap", theme_use = "theme_blank",legend.position = 'none')
vCM <- RunDynamicFeatures(srt = vCM, lineages = c("Lineage1", "Lineage2"), n_candidates = 200)
ht <- DynamicHeatmap(
  srt = vCM, lineages = c("Lineage1", "Lineage2"),pvalueCutoff = 0.05,
  padjustCutoff = 0.1,
  use_fitted = TRUE, n_split = 6, reverse_ht = "Lineage1",
  species = "Mus_musculus", db = "GO_BP", 
  heatmap_palette = "viridis",
  separate_annotation = list("subtype", "Esrra"), separate_annotation_palcolor= list(pal(4), pal(3)[2]),
  pseudotime_label = 25, pseudotime_label_color = "red",
  height = 2, width = 2
)

###########ATAC CoveragePlot############
#"ATP1B1","FLNC",'PPARA','ACO2','Acacb','Add3','Ldb3','Mylk4','Myh6'
library(RColorBrewer)
library(ggplot2)
library(Signac)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
DefaultAssay(obj) <- 'ATAC'
p1 <- CoveragePlot(
  assay = 'ATAC',
  object =obj,
  region = 'Atp1b1',
  group.by='group',
  peaks = F,
  annotation = FALSE
) + scale_fill_manual(values = pal(4))
p2 <- AnnotationPlot(
  object = obj,
  region = 'Atp1b1',
  assay = 'ATAC'
)+ scale_x_continuous(
  breaks = scales::breaks_extended(n = 3),
  labels =  function(x) sprintf("%.1f Mb", x / 1e6)
) + scale_color_manual(values = c('black','black'))
CombineTracks(
  plotlist = list(p1, p2),
  heights = c(10, 1),
  widths = c(10,10)
)

###########DEG num##########
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
for(i in celltypes){
  obj <- subset(vCM,subtype %in% 'vCM-2')
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
  
  ggsave(paste('vCM-2_组间DEG数.pdf',sep = ''),p,
         width = 3.3,height = 3)
} 

#atac
setwd("E:/2.工作/YF/mHeart/8.终稿调整/202504/各细胞类型差异peak数")
DEG_nb <- data.frame()
  obj <- subset(vCM,subtype %in% 'vCM-2')
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
  
  ggsave(paste('vCM-2_ATAC差异peak数.pdf',sep = ''),p,
         width = 3.3,height = 3)


##########cuttag######
library(ggbio)
library(rtracklayer)
library(GenomicRanges)
library(biomaRt)
library(dplyr)
library(tidyr)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#"ATP1B1","FLNC",'PPARA','ACO2','Acacb','Add3','Ldb3','Mylk4','Myh6'
# 连接到Ensembl数据库
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# 1. 设置工作目录并获取所有.bdg文件
setwd("G:/cuttag_3/peaks/bdg/CM")
sample_files <- list.files(pattern = "\\.bdg$", full.names = TRUE)
sample_names <-  sub("\\_treat_pileup.bdg$", "", basename(sample_files))
grlist <- lapply(sample_files, function(file) {
  gr <- import(file, format = "bedGraph")
  return(gr)
})
names(grlist ) <- sub("\\_treat_pileup.bdg$", "", basename(sample_files))
# 获取ESRRA基因坐标
target_gene <- 'MYH7'
region_info <- getBM(attributes=c("chromosome_name", "start_position", "end_position", "strand","hgnc_symbol"),
                      filters="hgnc_symbol",
                      values=c(target_gene),
                      mart=ensembl)
 # 格式化坐标
gene_regions <- GRanges(seqnames = paste0("chr", region_info$chromosome_name),
                        ranges = IRanges(start = region_info$start_position-20000,
                                         end = region_info$end_position),
                        gene = region_info$hgnc_symbol)
bdg_data <- lapply(grlist, function(file) {
  subsetByOverlaps(file, gene_regions)
})
names(bdg_data) <- sub("\\_treat_pileup.bdg$", "", basename(sample_files))
plot_data <- do.call(rbind, lapply(names(bdg_data), function(sample){
  gr <- bdg_data[[sample]]
  data.frame(
    sample = sample,
    seqnames = as.character(seqnames(gr)),
    start = start(gr),
    end = end(gr),
    score = gr$score
  )
}))
# plot
library(ggplot2)
library(scales)
library(org.Hs.eg.db)
gene_symbols <- select(org.Hs.eg.db, 
                       keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
                       columns = "SYMBOL",
                       keytype = "ENTREZID")
plot_data$group <- ifelse(plot_data$sample %in% c('AC-neg','CM-neg'),plot_data$sample,
                          ifelse(plot_data$sample %in% c('AC-P1','AC-P2','AC-P3-1','AC-P3-2'),'AC16(n=3)',
                                 ifelse(plot_data$sample %in% c('CM-P1','CM-P2','CM-P3'),'hiPSC-CM-ESRRA\n(n=3)','Other'
                                        )))
plot_data$group <-gsub('CM-neg','hiPSC-CM-IGg',plot_data$group)
gene_coords <- GRanges(seqnames = paste0("chr", region_info$chromosome_name),
                       ranges = IRanges(start = region_info$start_position,
                                        end = region_info$end_position),
                       gene = region_info$hgnc_symbol)

# 转换为数据框用于绘图
gene_df <- data.frame(
  seqnames = as.character(seqnames(gene_coords)),
  start = start(gene_coords),
  end = end(gene_coords),
  strand = as.character(strand(gene_coords)),
  gene = gene_coords$gene
)
# 计算y轴最大值
y_max <- max(plot_data$score, na.rm = TRUE)
cols <- c(pal(4)[1],'gray',pal(4)[2],'gray')
p <- ggplot(plot_data, aes(x = start, y = score)) +
  geom_area(aes(fill = group),alpha = 0.8) +
  geom_line(aes(color = group) ,alpha = 0.8,size = 0.4) +
  facet_grid(group ~., switch = "y") +
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)+
  scale_y_continuous(
    breaks = c(0, y_max),  # 只显示0和最大值
    limits = c(0, y_max)   # 确保y轴从0开始
  ) +
  scale_x_continuous(
    breaks = NULL,
    labels = NULL  # 使用科学计数法格式化横坐标
    # 或者使用: labels = scientific_format()  # 另一种写法
  ) +
  theme(
    # 分面标签设置
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    strip.background = element_blank(),
    # 去除所有网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 去除图表边框
    panel.background = element_blank(),  # 去除绘图区域背景
    panel.border = element_blank(),      # 去除绘图区域边框
    # 去除坐标轴线
    axis.line = element_line(),         # 去除坐标轴线
    axis.ticks = element_blank(),         # 保留刻度线
    axis.ticks.length = unit(0, "cm"), # 设置刻度线长度
    # 调整边距
    plot.margin = margin(0, 0, 0, 0),   # 调整图形边距
    # 图例设置
    legend.position = "none"
  )+
  labs(x='',y='')
p_gene <- ggplot(gene_df) +
  # 基因范围线段
  geom_segment(
    aes(x = start, xend = end, y = 0, yend = 0),
    color = "black",
    size = 1
  ) +
  # 箭头（方向由 strand 决定）
  geom_segment(
    aes(
      x = ifelse(strand == "+", end, start),
      xend = ifelse(strand == "+", end - 0.1 * (end - start), start + 0.1 * (end - start)),
      y = 0,
      yend = 0
    ),
    color = "black",
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    size = 1
  ) +
  # 基因名称
  geom_text(
    aes(x = (start + end) / 2, y = 0.3, label = gene),
    vjust = 0,
    size = 3
  ) +
  annotate("text",
           x = min(plot_data$start),  # 最左侧坐标
           y = 0.3,                   # y轴位置
           label = paste0("chr", region_info$chromosome_name),  # 染色体编号
           hjust = 0,                 # 左对齐
           vjust = 0,                 # 下对齐
           size = 3) +
  # 调整坐标轴和主题
  scale_x_continuous(
    limits = c(min(plot_data$start), max(plot_data$end)),  # 与信号图对齐
    breaks = scales::breaks_extended(n = 3),
    labels =  function(x) sprintf("%.1f Mb", x / 1e6)
  ) +
  scale_y_continuous(limits = c(-0.5, 1)) +
  theme_void() +  # 空白背景
  theme(
    axis.line.x = element_line(),
    axis.ticks.x = element_line(),
    axis.text.x = element_text(size = 8),
    plot.margin = margin(0, 0, 0, 0), 
    axis.ticks.length = unit(0.1, "cm")
  ) +
  labs(x = "Genomic Position", y = "")
library(patchwork)
 p / p_gene + 
  plot_layout(
    heights = c(2, 0.3),       # 信号图占4份高度，基因图占0.1份
    guides = "collect"       # 统一图例（如果有）
  ) & 
  theme(
    plot.margin = margin(0, 10, 0, 0)  # 完全去除图形之间的空白
  )
         
 
########fig1e deg+atac num#####3*3######
 library(BiocParallel)
 library(dplyr)
 library(patchwork)
 library(pheatmap)
 library(ggrepel)
 register(BPPARAM = SnowParam(1))
 setwd("E:/2.工作/YF/mHeart/8.终稿调整/202504/RNA和ATAC相关性图")
 rna.dt <- AggregateExpression(vCM, assays = 'RNA',group.by = 'subtype')[[1]]
 act.dt <- AggregateExpression(vCM, assays = 'activities',group.by = 'subtype')[[1]]
 genes.dt <- intersect(row.names(rna.dt),row.names(act.dt))
 p.dt <- data.frame(rna = rna.dt[genes.dt,2],activities = act.dt[genes.dt,2])
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
 
############################
library(Hmisc)
library(Signac)
library(FigR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BiocParallel)
vCM <- AnnotateFeatures(vCM, species = "Mus_musculus", IDtype = "symbol",
                        db = c("Chromosome", "GeneType", "Enzyme", "TF", "CSPA", "VerSeDa"))
genes <- rownames(vCM)[ subset(vCM,subtype %in% 'vCM-2')@assays$RNA@meta.features$TF %in% 'TF cofactor']
dt <- t(as.matrix(GetAssayData(subset(vCM)[c('Esrra',genes),],assay = 'RNA',slot = 'data')))
cor_results <- cor(dt[,-1 ], dt[,1 ], use = "complete.obs") 
writexl::write_xlsx(data.frame(genes = colnames(dt)[-1],cor=cor_results),'Esrra_cor.xlsx')
hetmap <- GroupHeatmap(vCM ,features = genes ,group.by = 'group',cluster_rows = T)
