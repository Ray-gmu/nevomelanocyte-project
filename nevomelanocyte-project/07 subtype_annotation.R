
seurat <- qread("seurat_celltype.qs")


Idents(seurat) <- "celltype"
table(Idents(seurat))

CellDimPlot(seurat, 
            group.by = "celltype", 
            reduction = "UMAP")
ggsave("celltype.umap.pdf",width =9,height = 7)

CellDimPlot(seurat, 
            group.by = "celltype", 
            reduction = "TSNE")
ggsave("celltype.tsne.pdf",width =9,height = 7)


subset <- subset(seurat, subset = celltype == "Melanocytes")


table(subset$celltype)
DimPlot(subset, 
        group.by = "celltype", 
        reduction = "umap")


genes_to_check = c("MITF","MLANA","TYR","PMEL","TYRP1","DCT") #Melanocytes

folder_name <- "featureplot"
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

for (gene in genes_to_check) {
  
  if (gene %in% rownames(subset)) {
    p <- FeaturePlot(subset, features = gene) + ggtitle(gene)
    
    
    ggsave(
      filename = paste0(folder_name, "/", gene, ".png"),
      plot = p,
      width = 6, height = 5, units = "in", dpi = 300
    )
    message("已保存: ", gene)
  } else {
    message("跳过: ", gene, " (基因不在矩阵中)")
  }
}




subset <- NormalizeData(subset, normalization.method = "LogNormalize", scale.factor = 10000)
subset <- FindVariableFeatures(subset, selection.method = "vst", nfeatures = 2000)
subset <- ScaleData(subset, features = VariableFeatures(subset))
subset <- RunPCA(subset, npcs = 50, verbose = FALSE)


top10 <- head(VariableFeatures(subset), 10)

plot1 <- VariableFeaturePlot(subset)
plot2 <- LabelPoints(plot = plot1, points = top10, 
                     repel = TRUE);plot2
ggsave("Variablegenes.pdf",plot = plot2, width = 8,height = 5)


marker <- c("TYR","PMEL","DCT","MLANA") # Melanocytes
pdf(file="VariableFeaturePlotmarker.pdf",width = 8, height = 5)
LabelPoints(plot = VariableFeaturePlot(object = subset), points = marker, repel = TRUE)
dev.off()


DimPlot(subset, reduction = "pca", group.by = "orig.ident")
ggsave("pc_DimPlot.png",width = 8,height = 6)

subset@reductions$pca@cell.embeddings[1:4,1:4]


p <- VizDimLoadings(subset, dims = 1:4, reduction = "pca")
ggsave("pc_VizDimLoadings.png",width = 8,height = 16)


png("pc_DimHeatmap.png", width = 2100, height = 1800, res = 300)
p <- DimHeatmap(subset, dims = c(1:3,6:8,12:14), cells = 500, balanced = TRUE, combine = TRUE)
dev.off()


ElbowPlot(subset,ndims=50)
ggsave("ElbowPlot.png",width = 8,height = 6)


library(findPC)
stdev_sorted <- sort(subset@reductions$pca@stdev, decreasing = TRUE) # 主成分降序排列
stdev_sorted


findPC(sdev = stdev_sorted, 
       number = 50, 
       method = "all", 
       figure = T)
#harmony
subset <- RunHarmony(subset,group.by.vars ="orig.ident") 
#FindNeighbors、FindClusters
dims_N <- "12" 
subset <- FindNeighbors(subset, dims = 1:dims_N, 
                        reduction = "harmony", 
                        verbose = FALSE) 

subset <- FindClusters(subset, 
                       resolution = seq(0.1, 2, by = 0.1), 
                       verbose = FALSE)
head(Idents(subset), 5)


library(clustree)
library(patchwork)
p1 <- clustree(subset, prefix = 'RNA_snn_res.') + coord_flip();p1
ggsave("clustertree.png",plot = p1, width = 20,height = 12)


p2 <- DimPlot(subset, group.by = 'RNA_snn_res.0.2', label = T);p2 
Tree <- p1 + p2 + plot_layout(widths = c(3, 1));print(Tree)
ggsave("clustertree_p2.png",plot = Tree, width = 20,height = 12)


subset$seurat_clusters <- subset$RNA_snn_res.0.2

#Umap
subset <- RunUMAP(subset, dims = 1:dims_N,reduction = "harmony",
                  umap.method = "uwot", metric = "cosine")
umap <- DimPlot(subset,
                reduction = 'umap',
                label = T,
                group.by = "seurat_clusters");umap
ggsave("UMAP.pdf",plot = umap,width = 8,height = 6)

#Tsne
subset <- RunTSNE(subset,dims = 1:dims_N,reduction = "harmony")
tsne <- DimPlot(subset,
                label = T,
                reduction = 'tsne',
                pt.size =1,
                group.by = "seurat_clusters");tsne
ggsave("TSNE.pdf",plot = tsne,width = 8,height = 6)



qsave(subset,"subset_Melanocytes.qs")

#FindAllmarker

Idents(subset) <- "seurat_clusters"
table(Idents(subset))
sub.markers <- FindAllMarkers(subset, 
                              only.pos = TRUE,
                              logfc.threshold = 0.5,
                              min.pct=0.1)
write.csv(sub.markers,file="markers.RNA.csv")


top10_markers <- sub.markers %>% dplyr::group_by(cluster) %>%
  dplyr::slice_head(n = 10)
print(top10_markers,n = 100)


DotPlot(subset, features = top10_markers$gene %>% unique()) +
  RotatedAxis() + ggplot2::coord_flip() +
  ggplot2::theme(axis.text.y = element_text(size = 7))
ggsave(
  filename = "top10_Markers_DotPlot.pdf", 
  # plot = p, 
  width = 8,    
  height = 12   
)



cluster_ids <- sort(unique(subset@meta.data$seurat_clusters))


sub_celltype <- data.frame(
  ClusterID = cluster_ids,
  sub_celltype = rep("Unknown", length(cluster_ids))
)


sub_celltype[sub_celltype$ClusterID %in% c(0,1 ),2]='MX1+_MX2+_Melanocytes'
sub_celltype[sub_celltype$ClusterID %in% c(2 ),2]='Proliferative_Melanocytes'
sub_celltype[sub_celltype$ClusterID %in% c(3 ),2]='TGFβ_responsive_Melanocytes'
sub_celltype[sub_celltype$ClusterID %in% c(4 ),2]='ECM_Melanocytes'



table(sub_celltype$sub_celltype)

table(subset@meta.data$seurat_clusters)

subset@meta.data$sub_celltype = "NA"

for(i in 1:nrow(sub_celltype)){
  subset@meta.data[which(subset@meta.data$seurat_clusters == sub_celltype$ClusterID[i]),'sub_celltype'] <- sub_celltype$sub_celltype[i]}

table(subset@meta.data$sub_celltype)


CellDimPlot(subset, group.by = "sub_celltype", reduction = "UMAP")
ggsave("CellDimPlot.celltype.umap.pdf",width =9,height = 7)

CellDimPlot(subset, group.by = "sub_celltype", reduction = "TSNE")
ggsave("CellDimPlot.celltype.tsne.pdf",width =9,height = 7)


CellDimPlot(subset, group.by = "orig.ident", reduction = "UMAP")
ggsave("orig.ident.pdf",width =9,height = 7)


#基因可视化----
table(Idents(subset))
Idents(subset) <- "sub_celltype"
genes_to_check = c("MX1","MX2", #MX1+_MX2+_Melanocytes
                   "PCNA","MKI67","TOP2A", #Proliferative_Melanocytes
                   "HLA-G", "ABCB5","HSPA1A", #TGFβ_responsive_Melanocytes
                   "MMP1","MMP3","COL1A1","COL3A1","COL6A1" #ECM_Melanocytes
                   )

############ 横坐标为celltype纵坐标表达量的小提琴图
for (gene in genes_to_check) {
  p <- VlnPlot(
    subset,
    features = gene,
        pt.size = 0,
    cols = c("#BEB8DC","#82B0D2","#C1D37F","#F6C1A0", "#E098AE" )
  ) + 
    geom_boxplot(width = 0.12, outlier.shape = NA) + 
    ggtitle(gene)+
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),  # 调整主标题字体大小，并使其居中
      axis.title.x = element_blank(),             # 调整 X 轴标题字体大小 (去掉 "Identity")
      axis.title.y = element_text(size = 10),             # 调整 Y 轴标题字体大小 (如 "Expression Level")
      axis.text.x = element_text(size = 6),              # 调整 X 轴刻度标签字体大小 (如 "0, 1, 2, 3, 4")
      axis.text.y = element_text(size = 10),              # 调整 Y 轴刻度标签字体大小 (如具体的表达量数值)
      legend.text = element_text(size = 8)               # 调整右侧图例文字的字体大小 (如果未将其隐藏)
      )
  
  # 保存为 PNG（高分辨率）
  ggsave(
    filename = paste0(gene, ".png"),
    plot = p,
    width = 6, height = 5, units = "in", dpi = 300
  )
}


####### 气泡图 #横纵左边反转
DotPlot(subset, features = genes_to_check) + RotatedAxis()
ggsave("DotPlot.sub_celltype.pdf",width =9,height = 7)

#热图
DoHeatmap(subset, features = genes_to_check, group.by = "sub_celltype")
ggsave("DoHeatmap.sub_celltype.pdf",width =9,height = 7)

#等比例热图----
# 1. 确保当前的 Idents 是 celltype
Idents(subset) <- "sub_celltype"
table(Idents(subset))
# 2. 【关键步骤】每个组随机抽取 300 个细胞
# 如果某个组细胞少于 300，则会提取该组全部细胞
subset_sub <- subset(seurat, downsample = 300)

# 3. 绘制热图
# 使用抽样后的对象 seurat_sub
p <- DoHeatmap(subset_sub, 
               features = genes_to_check, 
               group.by = "sub_celltype",
               size = 4,           # 调大上方分组文字的大小
               angle = 0,          # 调整文字角度，0度更易读
               label = TRUE,       # 显示标签
               draw.lines = TRUE)  # 在组间画分割线

# 4. 保存图片
ggsave("DoHeatmap_balanced.pdf", plot = p, width = 9, height = 7)

# scop包热图
ht <- GroupHeatmap(
  srt = subset,
  features = genes_to_check,
  group.by = c("sub_celltype"),
  heatmap_palette = "YlOrRd",
  # cell_annotation = c("Phase", "G2M_score", "Cdh2"),
  # cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  # show_row_names = TRUE,
  row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE
)
pdf("GroupHeatmap_res.pdf", width = 10, height = 8)
print(ht$plot)
dev.off()

# averageHeatmap热图
# 按照上面的顺序对scRNA中的celltype进行重新排序
# 在热图中会呈现对应的顺序
subset$sub_celltype <- factor(subset$sub_celltype,
                          levels = c("MX1+_MX2+_Melanocytes", "Proliferative_Melanocytes",
                                     "TGFβ_responsive_Melanocytes", "ECM_Melanocytes"))
Idents(subset) <- "sub_celltype"
# 热图
pdf("average_heatmap.pdf", width = 8, height = 10)
averageHeatmap(object = subset,
               markerGene = genes_to_check,
               group.by = "sub_celltype",
               cluster_rows = FALSE)
dev.off()


#最终结果展示----(如果选择去除环境RNA污染的数据则使用seuratnew)----
CellDimPlot(subset, 
            group.by = "sub_celltype", 
            reduction = "UMAP")
ggsave("final_sub_celltype.umap.pdf",width =9,height = 7)

CellDimPlot(subset, 
            group.by = "sub_celltype", 
            reduction = "TSNE")
ggsave("final_sub_celltype.tsne.pdf",width =9,height = 7)

qsave(subset, "sub_celltype.qs")