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
###
raw.data <- read.csv("database/GEOdata/GSE121893/GSE121893_human_heart_sc_umi.csv.gz", header=T, row.names = 1)
metadata <- read.table("database/GEOdata/GSE121893/GSE121893_all_heart_cell_cluster_info.txt.gz", header=T, row.names = 1)
hHeart <- CreateSeuratObject(counts = raw.data,meta.data =metadata )
ct.dt <- unique(hHeart$ident)
hvCM <- subset(hHeart,ident %in% c("LV1","LV2","LV3","LV4","LV5",
                                   "LA1","LA2","LA3","LA4","LA5","AV"))
hvCM[["RNA"]] <- split(hvCM[["RNA"]], f = hvCM$sample)
hvCM <- hvCM %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
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
hvCM[['RNA']] <- as(hvCM[['RNA']] ,Class = 'Assay')
hvCM$group <- gsub('cHF_LV','cHF',hvCM$condition)
hvCM$group <- gsub('cHF_LA','cHF',hvCM$group)
hvCM$group <- gsub('dHF_LA','dHF',hvCM$group)
hvCM$group <- gsub('dHF_LV','dHF',hvCM$group)
hvCM$group <- gsub('N_LA','N',hvCM$group)
hvCM$group <- gsub('N_LV','N',hvCM$group)
CellDimPlot(
  srt = hvCM,group.by = 'ident',legend.position = 'None',bg_color = 'gray95',split.by = 'group',
  reduction = "umap", theme_use = "theme_blank",label =T,label_insitu = T,pt.size = 1,pt.alpha = 0.7,ncol = 3,
  label.fg = 'black',label.bg = 'grey90',label.bg.r = 0.1,title = '',palcolor = pal(4),raster = F
)
hvCM
FeatureDimPlot(
  srt = hvCM, features = c('ANKRD1',"NPPB","ACTA1","ESRRA"),legend.position = 'right', assay="RNA",slot='data',split.by = 'group',bg_color ='gray95',
  reduction = "umap", theme_use = "theme_blank",raster = FALSE,ncol =3
)
allmarker <- FindAllMarkers(hvCM,group.by = 'group')
hvCM <- CellCycleScoring(hvCM ,s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes,
                        set.ident = TRUE)
FeatureDimPlot(
  srt = hvCM, features = c('S.Score','G2M.Score'),legend.position = 'right', assay="RNA",slot='data',split.by = 'group',bg_color ='gray95',
  reduction = "umap", theme_use = "theme_blank",raster = FALSE,ncol =3
)
saveRDS(hvCM,'hvCM.Rds')
