library(Seurat)
library(ggplot2)
library(ggunchained) 
library(reshape2) 
library(SCP)
library(Signac)
set.seed(4180)
setwd("E:\\2.工作/YF/mHeart/")
mHeart <- readRDS("E:/2.工作/YF/mHeart/mHeart.Rds")
#########color
cols <- c( "#444576", "#4682B4", "#AEDEEE", "#FFD790","#FFA500","#C65762",'#FBDFDE', "#F6EFCF","#BCB99F","#4D4D4D")
pal <- colorRampPalette(cols)
######fig1b uamp############
CellDimPlot(
  srt = mHeart,group.by = 'Phase',legend.position = 'none',label_repel = T,label_point_size = 0.01,
  reduction = "ATAC.umap", label = T,label_insitu = T,theme_use = 'theme_blank',
  label.fg = 'black',label.bg = 'grey95',label.bg.r = 0.1,title = '',palcolor =pal(4) ,ncol = 1,
  bg_color = 'grey90',raster = F
)
#######fig1de cellstat#############
mHeart$group <- factor(mHeart$group ,levels = c('PP','LP','MP','NP'))
mHeart$celltype <- factor(mHeart$celltype ,levels = rev(levels(mHeart$celltype )))
CellStatPlot(mHeart, stat.by = "celltype",plot_type = "pie",palcolor = pal(13),xlab = '',ylab = '')
CellStatPlot(mHeart, stat.by = "celltype", group.by = "group", label = F,palcolor = rev(pal(13)),
             flip = F,plot_type = "trend",xlab = '',ylab = 'Percentage of cell types')+coord_flip()
########fig1e deg+atac num#####3*3######
library(BiocParallel)
library(dplyr)
library(patchwork)
library(pheatmap)
library(ggrepel)
mHeart[['activities']] <- mHeart_ATAC[['activities']]
register(BPPARAM = SnowParam(1))
setwd("E:/2.工作/YF/mHeart/8.终稿调整/202504/RNA和ATAC相关性图")
rna.dt <- AggregateExpression(mHeart, assays = 'RNA',group.by = 'celltype')[[1]]
act.dt <- AggregateExpression(mHeart, assays = 'activities',group.by = 'celltype')[[1]]
genes.dt <- intersect(row.names(rna.dt),row.names(act.dt))
p.dt <- data.frame(rna = rowSums(rna.dt[genes.dt,]),activities = rowSums(act.dt[genes.dt,]))
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
############fig1f deg+atac num#####4.34*3.35##########
library(BiocParallel)
library(dplyr)
library(patchwork)
library(pheatmap)
library(ggrepel)
library(ggpubr)
mHeart[['activities']] <- mHeart_ATAC[['activities']]
register(BPPARAM = SnowParam(1))
mHeart$sample <- mHeart$orig.ident
celltypenames <-  unique(mHeart$celltype)
samplegroup <- c('NP','MP','LP','PP')
setwd("E:/2.工作/YF/mHeart/8.终稿调整/202504/RNA和ATAC相关性图")
rna.dt <- AggregateExpression(mHeart, assays = 'RNA',group.by = 'celltype')[[1]]
act.dt <- AggregateExpression(mHeart, assays = 'activities',group.by = 'celltype')[[1]]
genes.dt <- intersect(row.names(rna.dt),row.names(act.dt))
for (i in celltypenames){
  print(i)
  p.dt <- data.frame(rna = rna.dt[genes.dt,i],activities = act.dt[genes.dt,i])
  cor_value <- cor.test(log10(p.dt$rna + 1), log10(p.dt$activities + 1))
  r <- round(cor_value$estimate, 3)
  p_value <- scales::pvalue(cor_value$p.value)
  p <-  ggplot(data = p.dt, mapping = aes(x = log10(rna+1), y = log10(activities+1))) +
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
  ggsave(paste(i,'_RNA与ATAC相关性散点.pdf',sep = ''),p,
         width = 3,height = 3)
}
##############figS1p DEG bar plot#####4.71*4.57#####
library(BiocParallel)
library(dplyr)
library(patchwork)
library(clusterProfiler)
register(BPPARAM = SnowParam(1))
samplegroup <- c('NP','MP','LP','PP')
mHeart<- RunDEtest(mHeart, group_by = "group", only.pos = FALSE, fc.threshold = 1)
mHeart <- RunDEtest(srt = mHeart , group_by = "group",
                 fc.threshold = 1, only.pos = FALSE,min.pct = 0.1,
                 group1 = 'NP')
DEGs_group <- mHeart@tools$DEtest_custom$AllMarkers_wilcox
for (j in 2:4){
  mHeart <- RunDEtest(srt = mHeart , group_by = "group", 
                   fc.threshold = 1, only.pos = FALSE,min.pct = 0.1,
                   group1 = samplegroup[j],group2 = 'NP')
  DEGs_group <- rbind(DEGs_group,mHeart@tools$DEtest_custom$AllMarkers_wilcox)
}
mHeart@tools$DEtest_group$AllMarkers_wilcox <- DEGs_group
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
  qvalueCutoff = 0.05
)
res <- data_GO@compareClusterResult
enrich <- subset(res,Description %in% c('myeloid leukocyte migration','vascular endothelial growth factor signaling pathway','endothelial cell migration',
                                        'cellular response to starvation','cellular response to nutrient levels','autophagy of mitochondrion',
                                        'striated muscle cell differentiation','myofibril assembly','muscle cell development',
                                        'fatty acid beta-oxidation','regulation of heart contraction','sarcomere organization'))
colnames(enrich ) <- gsub('Groups','Cluster',colnames(enrich))
dt <- enrich
dt$order <- gsub('NP','1',dt$Cluster)
dt$order <- gsub('MP','2',dt$order)
dt$order <- gsub('LP','3',dt$order)
dt$order <- gsub('PP','4',dt$order)
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
  labs(x = "-Log10(pvalue)", y = "PP           LP          MP        NP ", 
       title = "GOBP pathway enrichment") + 
  # x = 0.61 用数值向量控制文本标签起始位置
  geom_text(size=4, aes(x = 0.05, label = Description), hjust = 0) + # hjust = 0,左对齐
  theme_classic() + 
  mytheme +
  NoLegend()
############DEG number########
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
  obj <- subset(mHeart,celltype %in% i)
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
  
  ggsave(paste(i,'_组间DEG数.pdf',sep = ''),p,
         width = 3.3,height = 3)
} 

#atac
setwd("E:/2.工作/YF/mHeart/8.终稿调整/202504/各细胞类型差异peak数")
DEG_nb <- data.frame()
for(i in celltypes){
  obj <- subset(mHeart,celltype %in% i)
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
 
  ggsave(paste(i,'_ATAC差异peak数.pdf',sep = ''),p,
         width = 3.3,height = 3)
}

############all celltypes workflow####
celltypenames <- 'aCM'
celltypenames <-  levels(mHeart$celltype)
celltypenames <- celltypenames[!celltypenames %in%  c('vCM','FB','EC','Macrophage','T','B','DC','SMC','Adipocyte')]
samplegroup <- c('NP','MP','LP','PP')
DefaultAssay(mHeart) <- 'RNA'
#组间差异
library(BiocParallel)
library(dplyr)
library(patchwork)
library(pheatmap)
library(ggrepel)
mHeart$sample <- mHeart$orig.ident
register(BPPARAM = SnowParam(1))
setwd("E:/2.工作/YF/mHeart/8.终稿调整/20250330/各细胞类型组间差异")
for (i in celltypenames[-c(1:8)]){
  print(i)
  obj <- subset(mHeart,celltype %in% i)
  obj <- RunDEtest(obj, group_by = "group", only.pos = FALSE, fc.threshold = 1)
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
  p_dt <- DEGs_group %>%
    filter(p_val_adj < 0.05) %>%
    filter(pct.1 > 0.25) %>%
    group_by(group1) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 10)
  obj$sample <- factor(obj$sample,levels = c('NP1','NP2','NP3',
                                             'MP1','MP2','MP3',
                                             'LP1','LP2','LP3',
                                             'PP1','PP2','PP3'))
  avexp <- as.matrix(AverageExpression(obj,assays = "RNA",slot = "data",group.by = 'sample')$RNA[p_dt$gene,])
  for (j in 2:4){
    #heat
    h_dt <- subset(p_dt,group1 %in% c('NP',samplegroup[j]))
    p1 <- pheatmap(subset(avexp[h_dt$gene,
                                c(grep("NP", colnames(avexp)),grep(samplegroup[j], colnames(avexp)))]), 
                   color = colorRampPalette(c("#92a5d1",'#c5dff4',"white",'#FBDFE2', "#E83945"))(100),  # 设置颜色渐变
                   scale = "row",  # 按行进行缩放
                   show_rownames = TRUE,  # 显示行名
                   show_colnames = TRUE,cluster_rows = F,cluster_cols = F,main = i,cellwidth = 15,  # 设置每个方格的宽度
                   cellheight = 15 )
    ggsave(paste(i,samplegroup[j],'vs_NP_差异基因热图.pdf',sep = '_'),p1,
           width =3,height =5
    )
    #vol
    vo_data <- subset(DEGs_group,group1 %in% samplegroup[j])
    vo_data$change <- as.factor(ifelse(vo_data$p_val_adj < 0.05 & abs(vo_data$avg_log2FC) > 0.5,
                                       ifelse(vo_data$avg_log2FC > 0,'UP','DOWN'),'NOT'))
    p2 <-  ggplot(data = vo_data, 
                  aes(x = avg_log2FC, y = -log10(p_val_adj),color=change))+
      geom_point(alpha=0.5, size = 2)+
      theme_bw()+theme(panel.border = element_rect(colour = "gray", fill = NA),  # 设置面板边框颜色变淡
                       panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank(),
                       axis.title.x.top = element_text(size = 2),
                       legend.position = 'none')+
      scale_color_manual(name = '', values = c("gray","#92a5d1" ,"#E83945"), limits = c('NOT','DOWN','UP'))+
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray")+
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray")+
      labs(title = i,
           x='log2FoldChange',
           subtitle = paste(samplegroup[j],'vs NP'))+
      geom_text_repel(data=subset(vo_data,gene %in% h_dt$gene,), aes(label=gene),show.legend = F)
    ggsave(paste(i,samplegroup[j],'vs_NP_差异基因火山图.pdf',sep = '_'),p2,
           width =4,height =4
    )
    #GO
    library(clusterProfiler)
    library(org.Mm.eg.db)
    group <- data.frame(gene=vo_data$gene,
                        group=vo_data$change)
    group <- subset(group,group != 'NOT')
    Gene_ID <- bitr(group$gene, fromType="SYMBOL", 
                    toType="ENTREZID", 
                    OrgDb="org.Mm.eg.db")
    data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
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
    p3 <- dotplot(data_GO, showCategory= 10,font.size = 10,label_format=50)+
      labs(title = i,
           subtitle = paste(samplegroup[j],'vs NP'))+
      scale_fill_continuous(high = "#92a5d1",
                            low = "#E83945")
    ggsave(paste(i,samplegroup[j],'vs_NP_差异基因GO通路.pdf',sep = '_'),p3,
           width =5,height =7
    )
  }  
}
#亚群差异
library(ComplexHeatmap)
setwd("E:/2.工作/YF/mHeart/8.终稿调整/202504/marker热图")
subtype.analysis <- function(seurat.obj,type.choose = NULL) { 
  data.list <- list()
  type.choose = type.choose
  for (i in 1:length(type.choose)) {
    print(type.choose[i])
    data <- subset(seurat.obj, celltype %in% type.choose[i])
    data[["RNA"]] <- CreateAssay5Object(GetAssayData(data,
                                                     assay = "RNA",
                                                     layer = "counts"
    ))
    data[["RNA"]] <- split(data[["RNA"]], f = data$orig.ident)
    data <- data %>%
      NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData(vars.to.regress = c("percent.mt", "percent.rb", "percent.hsp")) %>%
      RunPCA() %>%
      IntegrateLayers(.,
                      method = HarmonyIntegration,
                      orig.reduction = "pca",
                      new.reduction = "harmony",
                      verbose = FALSE
      ) %>%
      FindNeighbors(., reduction = "harmony", dims = 1:30) %>%
      FindClusters(., resolution = 0.3) %>%
      RunUMAP(.,
              dims = 1:30,
              reduction = "harmony",
              reduction.name = "umap"
      )
    data$subtype <- paste(type.choose[i], "-",
                          as.numeric(data$seurat_clusters),
                          sep = ""
    )
    data$subtype <- factor(data$subtype,levels = paste(type.choose[i],
                                                       sort(unique(as.numeric(data$seurat_clusters))),
                                                       sep = '-'))
    data <- JoinLayers(data)
    data[['RNA']] <- as(data[['RNA']],Class = 'Assay')
    data.list[[i]] <- data
  }
  names(data.list) <- type.choose
  return(data.list)
}
subsc.list <-  subtype.analysis(mHeart,
                                type.choose=celltypenames)
subsc.list <- c(subsc.list,EC=EC,vCM=vCM,FB=FB,Macrophsge=Mac)

for(i in 1:length(subsc.list)){
  palnumber <- length(unique(subsc.list[[i]]$subtype))
  p1 <- CellDimPlot(
    srt = subsc.list[[i]],group.by = 'subtype',legend.position = 'right',
    reduction = "umap", label = F,pt.size = 1,theme_use = 'theme_scp',
    palcolor =pal(palnumber))
  
  #ggsave(paste(names(subsc.list[i]),'_UMAP图.pdf',sep = ''), p1,width = 4.5,height = 3)  
  ###4.5*3
  p2 <- CellStatPlot(subsc.list[[i]], stat.by ='subtype' ,group.by ='group',palcolor = pal(palnumber),
                     label = T,position = "dodge",label.size = 3,label.bg.r = 0.1)###8*3
  #ggsave(paste(names(subsc.list[i]),'_各组中每个亚群所占比例bar图.pdf',sep = ''),p2,width = 8,height = 3)  
  subsc.list[[i]] <- RunDEtest(subsc.list[[i]], group_by = "subtype",fc.threshold = 1.3,only.pos = T,min.pct = 0.25)
  AllMarkers <- filter(subsc.list[[i]]@tools$DEtest_subtype$AllMarkers_wilcox, p_val_adj < 0.05,test_group_number ==1)
  AllMarkers <- AllMarkers %>%
    filter(p_val_adj < 0.05) %>%
    group_by(group1) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 50)
  AllMarkers$label <- rownames(AllMarkers)
  result <- AllMarkers %>%
    filter(p_val_adj < 0.05) %>%
    group_by(group1) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 5)
  ave.dt <- scale(t(as.matrix(AverageExpression(subsc.list[[i]][AllMarkers$gene,],group.by = 'subtype',assays = 'RNA',slot = 'data')$RNA)))
  ha = HeatmapAnnotation(foo = anno_mark(at = as.numeric(result$label) , 
                                   labels = result$gene,
                                   labels_gp = gpar(fontsize = 8,rot = 45)))
  pdf(paste(names(subsc.list[i]),'_各组中每个亚群maekr热图.pdf',sep = ''),width = 6,height = 3) #保存为pdf文件
  draw(Heatmap(ave.dt,cluster_columns = F,cluster_rows = F,col = c('grey80','white',"#444576"),show_column_names = F,
  heatmap_legend_param = list(
  title = "Expression", at = c(-2, 0, 2), 
  labels = c("Min", "", "Max")),
  top_annotation = ha, row_names_side = 'left',
  width = unit(10, "cm"), height = unit(4, "cm")
  ))
dev.off()
writexl::write_xlsx(AllMarkers,
                      paste(names(subsc.list[i]),'_各亚群所有的差异基因表(RNA).xlsx',sep = ''))
}

#################UMAP+vln########
FeatureDimPlot(
  srt = subset(mHeart_ATAC),
  , features = c('Ppargc1a'), assay="motif",bg_color = 'gray95',split.by = 'group',
  reduction = "umap", theme_use = "theme_blank",raster = F,ncol = 4,cells = sample(colnames(mHeart),10000))
FeatureStatPlot(subset(mHeart),group.by='celltype',assay = 'RNA',
                stat.by = c('S.Score'),
                xlab = '',palcolor = pal(13),
                box_color = 'grey40',box_width = 0.05,
                legend.position = 'none',plot_type = 'box',
                add_box = TRUE, stack = TRUE,ylab = '')

#######pesudo bulk##########
library(rgl)
library(mixOmics)
library(plotly)
avg <- AggregateExpression(mHeart,group.by = 'orig.ident')
pca_result <- prcomp(t(avg$RNA), scale. = F)
group <- factor(c(rep("LP", 3), rep("MP", 3), rep("NP", 3),rep("PP", 3)))
colors <- c("LP" = "#C65762", "MP" = "#FFD790", "NP" = "#4682B4",'PP'='#4D4D4D')
sample_colors <- colors[group]
pc_scores <- pca_result$x[, 1:3]
plot3d(pc_scores[, 1], pc_scores[, 2], pc_scores[, 3], 
       type = "s", size = 2, col = sample_colors, alpha = 0.8,  # 调整透明度和大小
       xlab = "PC1", ylab = "PC2", zlab = "PC3",
       box = TRUE, axes = TRUE, grid = TRUE)  # 添加网格和背景

# 添加图例
legend3d("topright", legend = levels(group), col = colors, pch = 16, cex = 1.5)

# 调整视角
view3d(theta = 30, phi = 20, zoom = 0.8)

###########DEG volcanol###########
library(BiocParallel)
library(dplyr)
library(patchwork)
library(pheatmap)
library(ggrepel)
library(rtracklayer)
register(BPPARAM = SnowParam(1))
mHeart$sample <- mHeart$orig.ident
celltypenames <-  levels(mHeart$celltype)
samplegroup <- c('NP','MP','LP','PP')
DefaultAssay(mHeart) <- 'RNA'
setwd("E:/2.工作/YF/mHeart/8.终稿调整/20250325/差异基因火山图")
# 读取GTF文件
gtf_file <- "E:/2.工作/YF/mHeart/database/genes.gtf.gz"
gtf_data <- import(gtf_file)  
gene_data <- as.data.frame(gtf_data@elementMetadata@listData)
lncgenes <- subset(gene_data,gene_type %in% "lncRNA")$gene_name
for (i in celltypenames[-c(6,8,11)]){
  print(i)
  obj <- subset(mHeart,celltype %in% i)
  obj <- RunDEtest(obj, group_by = "group", only.pos = FALSE, fc.threshold = 1)
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
  DEGs_group <- subset(DEGs_group,!gene %in% lncgenes)
    #vol
    vo_data <- DEGs_group
    vo_data$change <- as.factor(ifelse(vo_data$p_val_adj < 0.05 & abs(vo_data$avg_log2FC) > 0.5,
                                       ifelse(vo_data$avg_log2FC > 0,'UP','DOWN'),'NOT'))
    p_dt <- vo_data %>%
      filter(p_val_adj < 0.05) %>%  # 筛选显著差异基因
      filter(pct.1 > 0.25) %>%      # 筛选 pct.1 大于 0.25 的基因
      group_by(group1) %>%          # 按 group1 分组
      arrange(desc(avg_log2FC)) %>%  # 按 abs(avg_log2FC) 从大到小排序
      slice(c(1:5, (n() - 4):n()))  # 提取每组的前五个和后五个
    
    p <-  ggplot() +
      geom_point(data = vo_data,
                 aes(x = avg_log2FC, y = -log10(p_val_adj), color =change ,alpha = 0.5)) +
      coord_flip()+ #坐标轴翻转
      facet_grid(. ~ group1)+
      geom_vline(xintercept = c(-0.5, 0.5), 
                 size = 0.5, color = "grey50", lty = 'dashed')+ #添加阈值线
      scale_color_manual(values = c("#4682B4","gray80","#C65762"))+ #更改配色
      theme_bw()+
      theme(
        legend.position = 'none', #去掉图例
        panel.grid = element_blank(), #去掉背景网格
        axis.text = element_text(size = 10), #坐标轴标签大小
        axis.text.x = element_text(angle = 45, vjust = 0.8), #x轴标签旋转
        strip.text.x = element_text(size = 10, face = 'bold') #加粗分面标题
      )+
      geom_text_repel(
        data = p_dt,
        aes(x = avg_log2FC, y = -log10(p_val_adj), 
            label = gene, color = change),
        size = 3,
        fontface = 'bold.italic'
      )+labs(title = i)
    ggsave(paste(i,samplegroup[j],'vs_NP_差异基因火山图.pdf',sep = '_'),p,
           width =7,height =4.2
    )
  
}

########get markers##########
library(BiocParallel)
library(dplyr)
library(patchwork)
library(pheatmap)
library(ggrepel)
library(rtracklayer)
library(mascarade)
library(tidyverse)
library(data.table)
setwd("E:/2.工作/YF/mHeart/8.终稿调整/20250330")
register(BPPARAM = SnowParam(1))
subsc_list <- list(EC=EC,FB=FB,Mac=Mac,vCM=vCM)
celltypenames <-  names(subsc_list)
#读取GTF文件
gtf_file <- "E:/2.工作/YF/mHeart/database/genes.gtf.gz"
gtf_data <- import(gtf_file)  
gene_data <- as.data.frame(gtf_data@elementMetadata@listData)
lncgenes <- subset(gene_data,gene_type %in% "lncRNA")$gene_name
for (i in celltypenames){ 
  print(i)
  dir.create(i)
  obj <- subsc_list[[i]]
  obj <- RunDEtest(obj, group_by = "subtype", only.pos = T, fc.threshold = 1.1,min.pct = 0.3)
  DEGs_celltype <- obj@tools$DEtest_subtype$AllMarkers_wilcox
  DEGs_celltype <- subset(DEGs_celltype ,!gene %in% lncgenes)
  top_genes <- DEGs_celltype %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::filter(pct.1 > 0.3) %>%
    dplyr::filter(test_group_number == 1) %>%
    group_by(group1) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 6)
 # top_genes <- rbind(top_genes,subset(DEGs_celltype,gene %in% 'Esrra'))
 #top_genes <- top_genes[order(top_genes$group1),]
  for(j in unique(obj$subtype)){
    genelist <- subset(top_genes,group1 %in% j)$gene
    p<- FeatureDimPlot(
        srt = subset(obj),
        , features = c(genelist), assay="RNA",slot='data',bg_color = 'gray95',
        reduction = "umap", theme_use = "theme_blank",raster = F,ncol = 2,
        cells = sample(colnames(obj),5000))
    ggsave(paste(i,'/',j,'markers_RNA表达量UMAP图.pdf',sep = ''),p,
        width =5.3,height =6.3)
  p2 <- DotPlot(obj,dot.scale = 5,group.by = 'subtype',
          features = top_genes$gene,assay = 'RNA',
          cols = c('white',"#444576")
  ) + RotatedAxis() + # 来自Seurat
    theme(legend.position = 'right',
          legend.title = element_text(size = 8,angle = 90, vjust = 0.5, hjust = 0.5),
          legend.text = element_text(size = 8),
          legend.title.position = 'left',
          panel.border = element_rect(color = "black"),
          panel.spacing = unit(1, "mm"),
          axis.title = element_blank())+
    guides(colour = guide_colourbar(title.vjust = 0.1, title.hjust = 0,label.position = "right"),
           size = guide_legend(title.vjust = 0.1, title.hjust = 0,label.position = "right"))+
    labs(size = "Percent Expressed", color = "Average Expression")##7*3
  ggsave(paste(i,'/','markers_RNA表达量气泡图.pdf',sep = ''),p2,
         width =8.5,height =3)
  avexp <- AverageExpression(obj,assays = 'RNA',slot = 'data',group.by = 'subtype')[[1]][top_genes$gene,]
  top_genes <- cbind(top_genes,as.data.frame(avexp))
  writexl::write_xlsx(top_genes,
                      paste(i,'/',i,'亚群差异marker基因.xlsx',sep = ''))
  }
}
############gsea################
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

###########choose DEGs##########
library(BiocParallel)
library(dplyr)
library(patchwork)
library(pheatmap)
library(ggrepel)
library(rtracklayer)
library(mascarade)
library(tidyverse)
library(data.table)
register(BPPARAM = SnowParam(1))
mHeart[['activities']] <- mHeart_ATAC[['activities']]
mHeart[['motif']] <- mHeart_ATAC[['motif']]
mHeart[['ATAC']] <- mHeart_ATAC[['ATAC']]
mHeart <- RunDEtest(mHeart, group_by = "celltype", only.pos = T, fc.threshold = 1.5,min.pct = 0.3,assay = 'RNA')
DEG_RNA <- mHeart@tools$DEtest_celltype$AllMarkers_wilcox
top_genes_RNA <- DEG_RNA  %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.3) %>%
  dplyr::filter(gene %in% DEG_motif$gene)
mHeart <- RunDEtest(mHeart, group_by = "celltype", only.pos = T, fc.threshold = 1.5,min.pct = 0.3,assay = 'motif')
DEG_motif <- mHeart@tools$DEtest_celltype$AllMarkers_wilcox
top_genes_motif <- DEG_motif %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.3) %>%
  dplyr::filter(test_group_number == 1) %>%
  group_by(group1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 20)
gene <- c("Rxrg","Junb","Erg","Maf")
FeatureDimPlot(
  srt = subset(mHeart),
  , features = c('Ets1'), assay="motif",slot='data',bg_color = 'gray95',
  reduction = "umap", theme_use = "theme_blank",raster = F,ncol = 2,cells = sample(colnames(mHeart),20000))
MotifPlot(
  object = subset(mHeart),
  motifs =gene ,
  assay = 'ATAC'
)

###########G2m########
m.cc.genes <- readRDS("database/mouse_cell_cycle_genes.rds") 
mHeart <- CellCycleScoring(mHeart, s.features = m.cc.genes$s.genes, 
                       g2m.features = m.cc.genes$g2m.genes,
                       set.ident = TRUE)
##########cytotrace########
library(CytoTRACE2)
library(patchwork)
mHeart <- cytotrace2(mHeart,
                         is_seurat = TRUE,
                         slot_type = "counts",
                         species = "mouse",
                         seed = 1234)
annotation <- data.frame(phenotype=mHeart@meta.data$celltype) %>% set_rownames(., colnames(mHeart))
plots <- plotData(cytotrace2_result=mHeart, annotation=annotation, is_seurat=TRUE)
p1 <- plots$CytoTRACE2_UMAP
p2 <- plots$CytoTRACE2_Potency_UMAP
p3 <- plots$CytoTRACE2_Relative_UMAP
p4 <- plots$CytoTRACE2_Boxplot_byPheno
(p1+p2+p3+p4) + plot_layout(ncol = 2)

#########EpiTrace########
library(Seurat)
library(SeuratObject)
library(Signac)
library(EpiTrace)
initiated_peaks <- Init_Peakset(GetAssayData(vCM,assay = "ATAC",layer = "counts"))
initiated_peaks_df <- as.data.frame(initiated_peaks,row.names = NULL)

rownames(mtx@colData) -> cellname_vec
paste0(initiated_peaks_df$seqnames,'_',initiated_peaks_df$start,'_',initiated_peaks_df$end) -> initiated_peaks_df$peakName

as(assays(mtx)[['PeakMatrix']], "sparseMatrix") -> mtx2
initiated_mm <- Init_Matrix(cellname = cellname_vec,peakname = initiated_peaks_df$peakName,matrix = mtx2)

epitrace_obj_age_estimated <- EpiTraceAge_Convergence(initiated_peaks,initiated_mm,celltype = NULL,qualnum = 10,Z_cutoff = 2.5,mean_error_limit = 0.01,iterative_time = 20,parallel = T,ncore_lim = 46,ref_genome = 'hg19',non_standard_clock = F)

mtx@colData[epitrace_obj_age_estimated@meta.data$cell,]$seurat_clusters -> epitrace_obj_age_estimated@meta.data$archR_cluster
epitrace_obj_age_estimated@meta.data$archR_cluster <- factor(epitrace_obj_age_estimated@meta.data$archR_cluster,levels=paste0('C',c(1:10)))

epitrace_obj_age_estimated@meta.data %>% as.data.frame() -> epitrace_obj_age_estimated_meta
