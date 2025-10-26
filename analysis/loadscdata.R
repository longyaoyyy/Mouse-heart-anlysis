# 加载必要的库
library(Seurat)
library(parallel)
library(Signac)
library(ggplot2)
library(dplyr)
#########
library(EnsDb.Mmusculus.v79)
setwd("E:\\2.工作/YF/mHeart/")
plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2)
LoadScAtac <- function(sampleid = "control1",
                       data_dir = "path/to/sampleid",
                       group_prefixes = NULL,
                       group_names = NULL,
                       annotation = annotation,
                       combined.peaks = combined.peaks,
                       feature_name = c("percent.mt", "percent.rb", "percent.hsp"),
                       feature_patterns = c("^Mt-", "^Rp[sl]", "^Hsp")) {
  # 验证输入参数
  if (!is.character(sampleid) || length(sampleid) != 1) {
    stop("sampleid must be a single character string.")
  }
  if (!file.exists(file.path(data_dir, sampleid, "outs"))) {
    stop(paste0("Data directory does not exist: ", file.path(data_dir, sampleid, "outs")))
  }
  if (length(feature_name) != length(feature_patterns)) {
    stop("feature_name and feature_patterns must have the same length.")
  }

  # load the RNA and ATAC data
  counts <- tryCatch(Read10X_h5(file.path(data_dir, sampleid, "outs/filtered_feature_bc_matrix.h5")),
    error = function(e) stop("Failed to read 10X data.")
  )
  fragpath <- file.path(data_dir, sampleid, "outs/atac_fragments.tsv.gz")
  # create a Seurat object containing the RNA adata
  pbmc <- tryCatch(CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA",
    project = sampleid,
    min.cells = 3
  ), error = function(e) stop("Failed to create Seurat object."))

  # create ATAC assay and add it to the object
 frags.pbmc <- CreateFragmentObject(
  path = fragpath,
  cells = colnames(counts$`Gene Expression`)
)
 pbmc.ATAC <- FeatureMatrix(
  fragments = frags.pbmc,
  features = combined.peaks,
  cells = colnames(counts$`Gene Expression`)
)
  pbmc[["ATAC"]] <- tryCatch(CreateChromatinAssay(
    counts = pbmc.ATAC,
    fragments = frags.pbmc,
    annotation = annotation
  ), error = function(e) stop("Failed to create ATAC assay."))

  # Assign group based on prefixes, if provided
  if (!is.null(group_prefixes) && !is.null(group_names)) {
    safe_sampleid <- gsub("[^A-Za-z]", "", sampleid)
    pbmc$group <- group_names[which(group_prefixes %in% gsub("[0-9]", "", safe_sampleid))]
  }

  # Calculate feature percentages
  for (i in seq_along(feature_name)) {
    pbmc[[feature_name[i]]] <- tryCatch(PercentageFeatureSet(pbmc, pattern = feature_patterns[i]),
      error = function(e) stop(paste0("Failed to calculate feature percentage for: ", feature_name[i]))
    )
  }

  # Nucleosome signal and TSS enrichment calculations
  pbmc <- tryCatch(NucleosomeSignal(pbmc, assay = "ATAC"),
    error = function(e) stop("Failed to calculate nucleosome signal.")
  )
  pbmc <- tryCatch(TSSEnrichment(pbmc, assay = "ATAC"),
    error = function(e) stop("Failed to perform TSS enrichment analysis.")
  )
  pbmc$blacklist_ratio <- FractionCountsInRegion(
    object = pbmc,
    assay = "ATAC",
    regions = blacklist_mm10
  )
  return(pbmc)
}

sample_name <- c(
  "NP1", "NP2", "NP3",
  "MP1", "MP2", "MP3",
  "LP1", "LP2", "LP3",
  "PP1", "PP2", "PP3"
)
# get gene annotations
annotation <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation))
peak.list <- list()
for (i in 1:length(sample_name)) {
  # read in peak sets
  peaks.pbmc <- read.table(
    file = file.path("/mnt/i/mmCardiac", sample_name[i], "outs/atac_peaks.bed"),
    col.names = c("chr", "start", "end")
  )
  peak.list[i] <- makeGRangesFromDataFrame(peaks.pbmc)
}
# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(peak.list[[1]],peak.list[[2]],peak.list[[3]],
                               peak.list[[4]],peak.list[[5]],peak.list[[6]],
                               peak.list[[7]],peak.list[[8]],peak.list[[9]],
                               peak.list[[10]],peak.list[[11]],peak.list[[12]]))
# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]

sc_list <- lapply(sample_name, function(sample_id) {
  LoadScAtac(
    sampleid = sample_id,
    data_dir = "/mnt/i/mmCardiac",
    group_prefixes = c("NP", "MP", "LP", "PP"),
    group_names = c("NP", "MP", "LP", "PP"),
    annotation = annotation,
    combined.peaks = combined.peaks,
    feature_name = c("percent.mt", "percent.rb", "percent.hsp"),
    feature_patterns = c("^mt-", "^Rp[sl]", "^Hsp")
  )
})
names(sc_list) <- sample_name
#######################
library(decontX)
library(scDblFinder)
sc_list <- lapply(sc_list, FUN = function(x) {
  # decountX
  counts <- GetAssayData(x, assay = "RNA", layer = "counts")
  decontX_results <- decontX(counts)
  # scDblFinder
  RNA <- as.SingleCellExperiment(x, assay = "RNA")
  RNA <- scDblFinder(RNA)
  ATAC <- as.SingleCellExperiment(x, assay = "ATAC")
  ATAC <- scDblFinder(ATAC, aggregateFeatures=TRUE, nfeatures=25, processing="normFeatures")
  ## save results
  x$Contamination <- decontX_results$contamination
  x$RNA_Dbl_res <- RNA$scDblFinder.class
  x$RNA_Dbl_score <- RNA$scDblFinder.score
  x$ATAC_Dbl_res <- ATAC$scDblFinder.class
  x$ATAC_Dbl_score <- ATAC$scDblFinder.score
  
  return(x)
})

#######################
isOutlier <- function(qcidx,
                      type = c("both", "lower", "higher"),
                      threshold = 5) {
  library(outliers)
  type <- match.arg(type, c("both", "lower", "higher"))
  mad_value <- mad(qcidx, constant = 1.4826)
  # 计算每个数据点与中位数的绝对偏差
  deviations <- scale(qcidx, center = median(qcidx), scale = FALSE)
  if (type == "both") {
    out <- which(abs(deviations) > threshold * mad_value)
  } else if (type == "lower") {
    out <- which(deviations < -(threshold * mad_value))
  } else if (type == "higher") {
    out <- which(deviations > threshold * mad_value)
  }
  # 添加缺失值的索引并确保唯一性
  out <- unique(c(out, which(is.na(qcidx))))
  return(out)
}

for (i in 1:length(sample_name)) {
  dir.create(paste("QC/metrics/", sample_name[i], sep = ""))
  plotname <- c(
    "nCount_RNA",
    "nFeature_RNA",
    "nucleosome_signal",
    "Contamination"
  )
  for (j in 1:length(plotname)) {
    plotdata <- data.frame(ID = sc_list[[i]]@meta.data[, plotname[j]])
    outpt <- min(plotdata[, 1][isOutlier(plotdata[, 1], type = "higher")])
    plot <- plotdata |>
      ggplot(aes(x = ID)) +
      geom_density(alpha = 0.2) +
      ylab("Cell density") +
      xlab(plotname[j]) +
      geom_vline(xintercept = outpt) +
      geom_text(
        aes(
          x = outpt, y = Inf,
          label = paste("X =", outpt)
        ), # y=Inf使得文本尽可能往上显示
        vjust = 1, # 文本垂直对齐方式，1表示靠近上方
        hjust = 0, # 靠右
        inherit.aes = FALSE, # 不继承默认的aesthetics
        size = 3
      ) +
      theme_bw()
    ggsave(paste("QC/metrics/", sample_name[i], "/", plotname[j], ".pdf", sep = ""),
      plot,
      width = 5, height = 4
    )
  }
}
##
sc_list <- lapply(sc_list, function(x) {
  DefaultAssay(x) <- 'RNA'
  x$log10GenesPerUMI <- log10(x$nFeature_RNA)/log10(x$nCount_RNA)
  x$log10nCount_RNA <- log10(x$nCount_RNA)
  x$log10nFeature_RNA <- log10(x$nFeature_RNA)
  out.cell <- unique(c(
    isOutlier(x$log10nCount_RNA, type = "both"),
    isOutlier(x$log10nFeature_RNA, type = "both"),
    isOutlier(x$nucleosome_signal, type = "both"),
    isOutlier(x$TSS.enrichment, type = "both"),
    isOutlier(x$Contamination, type = "higher",threshold = 3),
    isOutlier(x$log10GenesPerUMI, type = "lower",threshold = 3)
  ))
  x <- x[, -out.cell]
  x <- subset(x, nFeature_RNA > 500 &
        RNA_Dbl_res %in% "singlet"&
        ATAC_Dbl_res %in% "singlet"  &
        Contamination  <  0.2  &
        percent.mt < 1 &
        percent.rb < 2 &
        percent.hsp < 1)

  return(x)
})
ATAC_list <- lapply(sc_list, function(x){
  DefaultAssay(x) <- "ATAC"
  x[['RNA']] <- NULL
  x <- x %>%
    FindTopFeatures( min.cutoff = 10)%>%
    RunTFIDF()%>%
   RunSVD()
  return(x)
}
)
for(i in 1:12){ 
  ATAC_list[[i]]@assays$ATAC@fragments[[1]]@path <- gsub('/mnt/i','G:', ATAC_list[[i]]@assays$ATAC@fragments[[1]]@path)
}
sc_list <- lapply(sc_list, function(x) {
  x[['ATAC']] <- NULL
  return(x)
})
####### integrate RNA###########
gc()
mHeart <- merge(sc_list[[1]], y = sc_list[-1], add.cell.ids = sample_name)
mHeart <- mHeart %>%
  NormalizeData(assay = 'RNA') %>%
  FindVariableFeatures(assay = 'RNA') %>%
  ScaleData(vars.to.regress = c("percent.mt", "percent.rb", "percent.hsp"),assay = 'RNA') %>%
  RunPCA(assay = 'RNA')
gc()
#####
mHeart <- IntegrateLayers(
  object = mHeart, method =  HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

############merge atac###################
mHeart_ATAC <- merge(ATAC_list[[1]], y = ATAC_list[-1], add.cell.ids = sample_name)
# process the combined dataset
mHeart_ATAC <- mHeart_ATAC %>%
  FindTopFeatures( min.cutoff = 10) %>%
  RunTFIDF() %>%
  RunSVD()

# find ATAC integration anchors
ATAC_list <-  lapply(ATAC_list, function(x) {
  colnames(x) <- paste(unique(x$orig.ident),colnames(x),sep = '_')
  return(x)
}
)
integration.anchors <- FindIntegrationAnchors(
  object.list = ATAC_list,
  anchor.features = rownames(ATAC_list[[1]]),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = mHeart_ATAC[["lsi"]],
  new.reduction.name = "ATAC_lsi",
  dims.to.integrate = 1:30
)
saveRDS(integrated,'ATAC.Rds')
####################
mHeart[["ATAC"]]  <- integrated[['ATAC']]
mHeart@reductions$rlsi <- integrated@reductions$ATAC_lsi
mHeart@reductions$lsi  <- mHeart_ATAC@reductions$lsi
mHeart[["RNA"]] <- JoinLayers(mHeart[["RNA"]])
#########
mHeart <-  FindMultiModalNeighbors(
  object = mHeart,
  reduction.list = list("harmony", "rlsi"), 
  dims.list = list(1:30, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
mHeart <- FindClusters(mHeart, resolution = 2,graph.name = 'wknn')
mHeart <- RunUMAP(mHeart, nn.name = "weighted.nn",assay = "RNA")
mHeart <- RunUMAP(mHeart,dims = 1:30, reduction = "harmony", reduction.name = "RNA.umap")
mHeart <- RunUMAP(mHeart,dims = 2:40, reduction = "rlsi", reduction.name = "ATAC.umap")
saveRDS(mHeart, "mHeart.rds")
######################
m.cc.genes <- readRDS("database/mouse_cell_cycle_genes.rds") 
mHeart <- CellCycleScoring(mHeart, s.features = m.cc.genes$s.genes, 
                             g2m.features = m.cc.genes$g2m.genes,
                             set.ident = TRUE)
DimPlot(mHeart,cols = pal(3),reduction = 'umap')+
tidydr::theme_dr(xlength = 0.2, 
 ylength = 0.2,
  arrow = arrow(length = unit(3, "mm"),type = "closed"))+
  theme(panel.grid = element_blank(), #移除背景网格线
        plot.title = element_text(hjust = 0.5,size = 16),
        axis.text.x = element_blank(), #x轴标签大小调整
        axis.text.y = element_blank(), #y轴标签大小调整
        axis.ticks = element_blank(),  #移除刻度
        axis.title.x = element_text(size = 8,hjust = 0.05), #x轴标题大小调整
        axis.title.y = element_text(size = 8,hjust = 0.05), #移除y轴标题
        legend.title = element_text(size = 12), #图例标题大小调整
        legend.text = element_text(size = 10),#图例标签大小调整
        legend.position="right")+
  guides(color = guide_legend(ncol =  1,override.aes = list(size = 4)))+ ##图例改成一行
  labs(title = "cell cycle",color= '',x='UMAP-1',y='UMAP-2')


