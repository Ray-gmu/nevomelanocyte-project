
seurat$detailed_celltype <- as.character(seurat$celltype)

subset_barcodes <- rownames(subset@meta.data)
subset_annotations <- as.character(subset$sub_celltype)

seurat@meta.data[subset_barcodes, "detailed_celltype"] <- subset_annotations

seurat$detailed_celltype <- as.factor(seurat$detailed_celltype)

table(seurat$detailed_celltype)

qsave(seurat, "seurat_withsubmc.qc")

CellDimPlot(subset, 
group.by = "sub_celltype", 
reduction = "UMAP")
ggsave("celltype.umap.pdf",width =10,height = 7)

CellDimPlot(subset, 
            group.by = "sub_celltype", 
            reduction = "TSNE")
ggsave("celltype.tsne.pdf",width =10,height = 7)
# decontX ----
library(decontX)
set.seed(123)

counts <- GetAssayData(object = subset, layer = "counts")

decontX_results <- decontX(x = counts,
                       z = subset$sub_celltype,
                       batch = subset$orig.ident)
# decontX_results <- decontX(counts)

subset$contamination <- decontX_results$contamination

summary(subset$contamination)

library(ggplot2)
FeaturePlot(subset, 
            features = 'contamination', 
            raster=FALSE,
            reduction = "tsne"
) + 
  scale_color_gradient2(low = "#E5554D", high = "#80B1D3", midpoint = 0.5)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab('UMAP_1')+
  ylab('UMAP_2')
ggsave("contamination.tsne.pdf",width =8,height = 7)

clean_counts <- round(decontX_results$decontXcounts)
subset[["DecontX_RNA"]] <- CreateAssayObject(counts = clean_counts)

DefaultAssay(subset) <- "DecontX_RNA"

subset <- NormalizeData(subset)

library(patchwork)
Idents(subset) <- "sub_celltype"

gene <- "MLANA"

{

DefaultAssay(subset) <- "RNA"
p1 <- FeaturePlot(subset, features = gene, label = FALSE, reduction = "umap") + 
  ggtitle(paste0(gene, " Before DecontX (Raw RNA)")) 

DefaultAssay(subset) <- "DecontX_RNA"
p2 <- FeaturePlot(subset, features = gene, label = FALSE, reduction = "umap") + 
  ggtitle(paste0(gene, " After DecontX (Cleaned)")) 

p1 | p2

ggsave(filename = paste0(gene, "_decontX_clean.pdf"), width = 16, height = 7)
}
# qsave(seurat, "seurat_clean.qs")
# seurat <- qread("seurat_clean.qs")

# DotPlot----

library(ggplot2)


DefaultAssay(subset) <- "DecontX_RNA"
DefaultAssay(subset) <- "RNA"

Idents(subset) <- "sub_celltype"


melanocyte_markers <- c("MITF", "MLANA", "PMEL", "TYR", "TYRP1")
ecm_markers <- c("COL1A1", "COL3A1", "COL6A1")


features_to_plot <- c(melanocyte_markers, ecm_markers)
dot_plot <- DotPlot(subset, features = features_to_plot, scale = FALSE) + 
  
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic"),
        axis.text.y = element_text(face = "bold")) +
  
  
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  
  
  labs(title = "Expression of Melanocyte and ECM markers (After DecontX)",
       x = "Marker Genes", 
       y = "Cell Types")


print(dot_plot)

DefaultAssay(subset) <- "DecontX_RNA"  
DefaultAssay(subset) <- "RNA"   

feature1 = "MLANA"
feature2 = "COL3A1"


{
#   FeatureScatter(subset, 
#                feature1 = feature1, 
#                feature2 = feature2, 
#                pt.size = 1.5,
#                group.by = "sub_celltype") +
#   geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") + 
#   geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
#   ggtitle(paste("Co-expression of", feature1, "and", feature2, "in Melanocytes")) +
#   theme_bw() +
#   scale_color_manual(values = custom_colors)
# ggsave(filename = paste("Co-expression of", feature1, "and", feature2, "in Melanocytes.pdf"), width = 8, height = 7)

FeatureScatter(subset, 
               feature1 = feature1, 
               feature2 = feature2, 
               pt.size = 1.5,
               group.by = "sub_celltype") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") + 
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  ggtitle(paste("Co-expression of", feature1, "and", feature2, "in Melanocytes")) +
  theme_bw()+
  scale_color_manual(values = custom_colors)+
    xlim(0, 6.5) + 
    ylim(0, 6)
ggsave(filename = paste("sub_celltype_Co-expression of", feature1, "and", feature2, "in Melanocytes.pdf"), width = 9, height = 7)

# FeatureScatter(subset, 
#                feature1 = feature1, 
#                feature2 = feature2, 
#                pt.size = 1.5,
#                group.by = "orig.ident") +
#   geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") + 
#   geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
#   ggtitle(paste("Co-expression of", feature1, "and", feature2, "in Melanocytes")) +
#   theme_bw()+
#   scale_color_manual(values = custom_colors)
# ggsave(filename = paste("sample Co-expression of", feature1, "and", feature2, "in Melanocytes.pdf"), width = 8, height = 7)

FeatureScatter(subset,
               feature1 = feature1,
               feature2 = feature2,
               pt.size = 1.5,
               group.by = "sub_celltype",
               split.by = "sub_celltype") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  ggtitle(paste("Co-expression of", feature1, "and", feature2, "in Melanocytes")) +
  theme_bw()+
  scale_color_manual(values = custom_colors)+
  xlim(0, 6.5) + 
  ylim(0, 6)
ggsave(filename = paste("split Co-expression of", feature1, "and", feature2, "in Melanocytes.pdf"), width = 14, height = 7)
}



#  (Co-expression Analysis) for----
DefaultAssay(subset) <- "DecontX_RNA"  
DefaultAssay(subset) <- "RNA"   


feature1 <- "MLANA"


feature2_list <- c("MMP1", "TNC", "COL6A1", "MMP3", "COL1A1", "DCN", "THBS1", "SERPINE1", "COL3A1")


for (feature2 in feature2_list) {
  
  
  cat("正在绘制并保存:", feature1, "与", feature2, "的共表达图...\n")
  
  # ==========================================
 
  # ==========================================
  p1 <- FeatureScatter(subset, 
                       feature1 = feature1, 
                       feature2 = feature2, 
                       pt.size = 1.5,
                       group.by = "sub_celltype") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") + 
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
    ggtitle(paste("Co-expression of", feature1, "and", feature2, "in Melanocytes")) +
    theme_bw() +
    scale_color_manual(values = custom_colors) +
    xlim(0, 4.5) + 
    ylim(0, 6)
  
  # 使用 paste0 构建严谨的文件名并明确指定 plot = p1
  filename1 <- paste0("sub_celltype_Co-expression_of_", feature1, "_and_", feature2, ".pdf")
  ggsave(filename = filename1, plot = p1, width = 9, height = 7)
  
  
  # ==========================================
  
  # ==========================================
  p2 <- FeatureScatter(subset, 
                       feature1 = feature1, 
                       feature2 = feature2, 
                       pt.size = 1.5,
                       group.by = "sub_celltype",
                       split.by = "sub_celltype") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") + 
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
    ggtitle(paste("Co-expression of", feature1, "and", feature2, "in Melanocytes")) +
    theme_bw() +
    scale_color_manual(values = custom_colors) +
    xlim(0, 4.5) + 
    ylim(0, 6)
  
  
  filename2 <- paste0("split_Co-expression_of_", feature1, "_and_", feature2, ".pdf")
  ggsave(filename = filename2, plot = p2, width = 14, height = 7)
  
}

cat("所有基因批量绘制完成！\n")
# FeaturePlot blend----

genes_to_plot <- c("MLANA", "COL1A1")


FeaturePlot(
  object = subset, 
  features = genes_to_plot, 
  blend = TRUE,              
  blend.threshold = 0.5,     
  order = TRUE,              
  pt.size = 0.5,             
  cols = c("lightgrey", "#ff0000", "#00ff00") 
)



# nCount_RNA nFeature_RNA----
VlnPlot(subset,
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt",
                     "percent.rp",
                     "percent.hb"),
        ncol = 3,pt.size = 0, group.by = "sub_celltype")
w=length(unique(seurat_all$orig.ident))/3+5;w
ggsave(filename="Vlnplot_nCount_RNA_nFeature_RNA.png",plot=p1,width = w,height = 5)