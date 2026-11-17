
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(ROGUE)
library(clustree)
library(harmony)
library(SingleR)
library(dplyr)
library(monocle)
library(tidyverse)
library(CellChat)
library(SCopeLoomR)
library(SCENIC)
library(AUCell)
library(foreach)
library(KernSmooth)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(pheatmap)
library(ggheatmap)
library(reshape2)
library(SummarizedExperiment)
library(ggpubr)
library(ggrepel)
library(SeuratWrappers)
library(SeuratData)
library(copykat)
library(limma)
library(GSEABase)
library(GSVA)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(stringi)	
library(GOplot)
R.utils::setOption("clusterProfiler.download.method",'auto')
library(monocle3)




setwd("C:/Users/wangr/Desktop/痣组织单细胞2/8 subtype/Melanocytes/monocle3")

subset <- readRDS("C:/Users/wangr/Desktop/痣组织单细胞2/8 subtype/Melanocytes/sub/subtype.rds")

colnames(subset@meta.data)

p = DimPlot(subset, group.by = "sub_celltype", label = T, reduction = "umap")
ggsave("DimPlot.UMAP.pdf", plot = p, width = 8, height = 5)
p = DimPlot(subset, group.by = "sub_celltype", label = T, reduction = "tsne")
ggsave("DimPlot.TSNE.pdf", plot = p, width = 8, height = 5)

data <- GetAssayData(subset, assay = 'RNA', layer  = 'counts')

cell_metadata <- subset@meta.data

gene_annotation <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))

cds <- new_cell_data_set(data, 
                         cell_metadata = cell_metadata, 
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, method = "PCA")

cds <- reduce_dimension(cds, reduction_method = "UMAP",preprocess_method = 'PCA')

p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters", show_trajectory_graph = FALSE) + ggtitle('cds.umap')
ggsave("plot_cells_umap.pdf", plot = p1, width = 8, height = 5)

int.embed <- Embeddings(subset, reduction = "umap")

int.embed <- int.embed[rownames(cds@int_colData$reducedDims$UMAP),]

cds@int_colData$reducedDims$UMAP <- int.embed

p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters", show_trajectory_graph = FALSE) + ggtitle('seurat.umap')
p = p1|p2
ggsave("2.Reduction_Compare.pdf", plot = p, width = 10, height = 5)


cds <- cluster_cells(cds)

p1 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + ggtitle("partition")
ggsave("3.Cluster_Partition.pdf", plot = p1, width = 6, height = 5)


cds <- learn_graph(cds, learn_graph_control = list(euclidean_distance_ratio = 0.8))

p = plot_cells(cds, color_cells_by = "partition",label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
ggsave("4.Trajectory.pdf", plot = p, width = 6, height = 5)


cds <- order_cells(cds)

p = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
               label_leaves = FALSE,  label_branch_points = FALSE)
ggsave("5.Trajectory_Pseudotime.pdf", plot = p, width = 8, height = 6)

saveRDS(cds, file = "cds.rds")


Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=1)

write.csv(Track_genes, "Track_genes_all.csv", row.names = F)



saveRDS(Track_genes, "Track_genes_all.rds")


genes <- row.names(subset(Track_genes, morans_I > 0.25))


genes = rownames(top_n(Track_genes,n=60,morans_I))


genes <- c(
  "CDK1", "CDK2", "MKI67", "UBE2C", "TOP2A", "PCNA",
  "CDKN1A", "CDKN2A", "GADD45B", "SOD2", "JUN", "FOS",
  "MITF", "DCT", "TYR", "PMEL", "MYO5A", "EDNRB",
  "CCN1", "MMP1", "MMP3", "COL1A1", "COL3A1", "COL6A1",
  "IL6","CXCL8","COL6A3"
)
genes <- c(
  "COL1A1","COL3A1","COL6A1","COL1A2","COL6A2","COL6A3",
  "DCN","THBS1","TNC","FN1","CEMIP","MMP1","SERPINE1",
  "TIMP1","ID1","ID2","ID3","STC1","VIM"
)
genes <- c(
  "COL1A1","COL3A1","COL6A1",
  "DCN","FN1","SDC4","CD44","CD47","TNC","DAG1",
  "LAMA2","LAMA4","LAMB1","LAMC1","LAMC3",
  "THBS1","THBS2","PTN",
  "ITGB1","ITGA2","ITGA3","ITGA9-AS1","ITGA5",
  "PTPRZ1","MMP1","MMP3"
)

genes <- c("KIT","SLC45A2","MYC","TYR","TYRP1","PMEL") # pigmentation
genes <- c("COL1A1","COL3A1","COL6A1","MMP1","MMP3","PDGFRA","SFRP2","CRISPLD2") # ecm
genes <- c("MKI67","TOP2A","PCNA","UBE2C") # cell cycle
genes <- c("COX5B","ATP5F1E","CYCS") # cell cycle

genes <- c("KIT","SLC45A2","MYC","TYR","TYRP1","PMEL","COL1A1",
           "COL3A1","COL6A1","MMP1","MMP3","PDGFRA","SFRP2","CRISPLD2",
           "MKI67","TOP2A","PCNA","UBE2C","COX5B","CYCS") 

genes <- c("TYR","TYRP1","PMEL","COL1A1",
           "COL3A1","COL6A1","MMP1","MMP3",
           "MKI67","TOP2A","PCNA") 


num_clusters = 3

pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
mycol = rev(RColorBrewer::brewer.pal(11, "Spectral"))
p = ComplexHeatmap::Heatmap(
  pt.matrix, name = "z-score", show_row_names = T, show_column_names = F,
  col = circlize::colorRamp2(seq(from=-2,to=2,length=11), mycol),
  row_names_gp = gpar(fontsize = 6), row_title_rot= 0, km = num_clusters, 
  cluster_rows = TRUE, cluster_row_slices = FALSE, cluster_columns = FALSE,use_raster=F
)
pdf("PseudotimeGenes_heatmap.pdf", width = 8, height = 8)
plot(p)
dev.off()


for(g in genes){
  p <- plot_genes_in_pseudotime(cds[g, ], 
                                color_cells_by = "sub_celltype", 
                                min_expr = 0.5) +
    ggtitle(g)
  
  ggsave(paste0(g, ".pdf"), 
         plot = p, width = 6, height = 5)
}

set.seed(666)

genelist <- row.names(subset(Track_genes, morans_I > 0.1))
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 6)
table(gene_module$module)
write.csv(gene_module, "12.PseudotimeGenes_Module.csv", row.names = F)


cell_group <- tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$sub_celltype)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- str_c("Module ", row.names(agg_mat))

cell_order <- c("ECM_Melanocytes", "MX1+_MX2+_Melanocytes", "TGFβ_responsive_Melanocytes", "Proliferative_Melanocytes")

agg_mat <- agg_mat[, cell_order]

p <- pheatmap::pheatmap(
  agg_mat,
  scale = "column",
  cluster_cols = FALSE,            
  clustering_method = "ward.D2"    
)
ggsave("13.PseudotimeGenes_Module.pdf", plot = p, width = 8, height = 6.5)


ME.list <- lapply(1:length(unique(gene_module$module)), function(i){subset(gene_module,module==i)$id})
names(ME.list) <- paste0("module", 1:length(unique(gene_module$module)))
subset <- AddModuleScore(subset, features = ME.list, name = "module")

feat_cols <- grep("^module\\d+(_1)?$", colnames(subset@meta.data), value = TRUE)

for (f in sort(feat_cols)) {
  p <- FeaturePlot(subset, features = f, raster = TRUE) +   # raster=TRUE 更快更小（Seurat v5可用）
    scale_color_viridis_c(option = "C") +
    ggtitle(f)
  ggsave(paste0("ModuleScore_", f, ".png"),
         plot = p, width = 6, height = 5, dpi = 300)
}

saveRDS(subset,"sub_celltype.rds")
subset <-readRDS("sub_celltype.rds")



library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(Matrix)
library(grid)


module <- 5   
genes_mod <- gene_module$id[gene_module$module == module]


mat <- tryCatch(normalized_counts(cds), error = function(e) exprs(cds))


genes_mod <- intersect(genes_mod, rownames(mat))
ord <- order(pseudotime(cds), na.last = NA)
pt <- mat[genes_mod, ord, drop = FALSE]

if (length(genes_mod) < 2 || ncol(pt) < 10) {
  stop(paste0("Module", module, " 可画基因/细胞太少：genes=", 
              length(genes_mod), ", cells=", ncol(pt),
              "。请检查基因是否存在或放宽筛选。"))
}


smooth_row <- function(x, df = 5) {
  if (sum(is.finite(x)) <= df) return(x)
  y <- tryCatch(stats::smooth.spline(seq_along(x), x, df = df)$y, error = function(e) x)
  y
}
pt.s <- t(apply(pt, 1, smooth_row))


pt.z <- t(scale(t(pt.s)))
pt.z[!is.finite(pt.z)] <- 0
rownames(pt.z) <- genes_mod


bins <- 60
if (ncol(pt.z) > bins) {
  brks <- cut(seq_len(ncol(pt.z)), breaks = bins, labels = FALSE)
  pt.z <- sapply(split(seq_len(ncol(pt.z)), brks),
                 function(idx) Matrix::rowMeans(pt.z[, idx, drop = FALSE]))
  pt.z <- as.matrix(pt.z)
}


mycol <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
ht <- Heatmap(
  pt.z,
  name = paste0("z-score (Module", module, ")"),
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  show_column_names = FALSE,
  col = colorRamp2(seq(-2, 2, length.out = 11), mycol),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  use_raster = TRUE,
  raster_device = "png",
  raster_quality = 2
)


outfile <- paste0("PseudotimeGenes_heatmap_Module", module, ".pdf")
pdf(outfile, width = 7, height = 8)
draw(ht)
dev.off()


library(clusterProfiler)
library(org.Hs.eg.db) 
library(enrichplot)
library(ggplot2)
library(dplyr)

module_id <- 13   

genes_module <- subset(gene_module, gene_module$module == module_id)$id

gene_map <- bitr(genes_module, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
entrez_genes <- gene_map$ENTREZID

ego <- enrichGO(gene=entrez_genes, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
                ont="BP", pvalueCutoff=0.05, qvalueCutoff=0.05, readable=TRUE)

ekegg <- enrichKEGG(gene=entrez_genes, organism="hsa", pvalueCutoff=0.05)
ekegg <- setReadable(ekegg, OrgDb=org.Hs.eg.db, keyType="ENTREZID")


go_title   <- paste0("Module", module_id, " - GO BP")
kegg_title <- paste0("Module", module_id, " - KEGG")

p1 <- dotplot(ego, showCategory=20) + ggtitle(go_title)
p2 <- dotplot(ekegg, showCategory=20) + ggtitle(kegg_title)

ggsave(paste0("Module", module_id, "_GO_BP_dotplot.pdf"), p1, width=8, height=10)
ggsave(paste0("Module", module_id, "_KEGG_dotplot.pdf"),  p2, width=8, height=10)

write.csv(as.data.frame(ego),   paste0("Module", module_id, "_GO_BP.csv"),  row.names=FALSE)
write.csv(as.data.frame(ekegg), paste0("Module", module_id, "_KEGG.csv"),   row.names=FALSE)


library(ggplot2)
library(ggridges)
library(dplyr)

plot_data <- data.frame(
  cell_id = colnames(cds),
  sub_celltype = colData(cds)$sub_celltype,
  pseudotime = pseudotime(cds)
)


plot_data <- plot_data %>%
  filter(!is.na(pseudotime) & is.finite(pseudotime))


plot_data$sub_celltype <- factor(
  plot_data$sub_celltype,
  levels = rev(c("Proliferative_Melanocytes", 
                 "TGFβ_responsive_Melanocytes", 
                 "MX1+_MX2+_Melanocytes", 
                 "ECM_Melanocytes")) # 用 rev() 反转一下，让发育早期的在最上面
)


p_ridge <- ggplot(plot_data, aes(x = pseudotime, y = sub_celltype, fill = sub_celltype)) +
  geom_density_ridges(
    scale = 1.5,       
    alpha = 0.8,       
    color = "black",   
    linewidth = 0.5    
  ) +
  
  scale_fill_manual(values = c(
    "Proliferative_Melanocytes" = "#BF812D",
    "TGFβ_responsive_Melanocytes" = "#DFC27D",
    "MX1+_MX2+_Melanocytes" = "#80CDC1",
    "ECM_Melanocytes" = "#01665E"
  )) +
  theme_classic() +
  labs(x = "Pseudotime", y = "") +
  theme(
    legend.position = "none", # 去除图例，Y轴标签已经说明了身份
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )


print(p_ridge)
ggsave("Ridgeline_Plot_Pseudotime.pdf", plot = p_ridge, width = 8, height = 4)


my_colors <- c(
  "Proliferative_Melanocytes" = "#BF812D",
  "TGFβ_responsive_Melanocytes" = "#DFC27D",
  "MX1+_MX2+_Melanocytes" = "#80CDC1",
  "ECM_Melanocytes" = "#01665E"
)


for(g in genes){
  p <- plot_genes_in_pseudotime(
    cds[g, ], 
    color_cells_by = "sub_celltype", 
    min_expr = 0.2,        
    cell_size = 0.8        
  ) +
    
    scale_color_manual(values = my_colors) +
    scale_fill_manual(values = my_colors) +
    
    theme_classic() +
    
    labs(x = "Pseudotime", y = paste(g, "Expression")) +
    ggtitle(g) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, color = "black"),
      axis.line = element_line(linewidth = 0.8),
      legend.position = "right",
      legend.title = element_blank() 
    )
  
  
  ggsave(paste0("Pseudotime_Curve_", g, ".pdf"), plot = p, width = 6, height = 4)
}


library(ggplot2)
library(patchwork)


my_colors <- c(
  "Proliferative_Melanocytes" = "#BF812D",
  "TGFβ_responsive_Melanocytes" = "#DFC27D",
  "MX1+_MX2+_Melanocytes" = "#80CDC1",
  "ECM_Melanocytes" = "#01665E"
)


my_theme <- theme_void() + 
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )


p_pseudo <- plot_cells(
  cds, 
  color_cells_by = "pseudotime", 
  label_cell_groups = FALSE, 
  label_leaves = FALSE, 
  label_branch_points = FALSE,
  cell_size = 0.8,                     
  trajectory_graph_color = "black",    
  trajectory_graph_segment_size = 1.2  
) + 
  my_theme + 
  ggtitle("Pseudotime Trajectory")

ggsave("Trajectory_Pseudotime.pdf", plot = p_pseudo, width = 8, height = 6)



p_traj <- plot_cells(
  cds, 
  color_cells_by = "sub_celltype", 
  label_cell_groups = FALSE,           
  label_leaves = FALSE, 
  label_branch_points = FALSE,
  cell_size = 0.8,                     
  trajectory_graph_color = "black", 
  trajectory_graph_segment_size = 1.2 
) + 
  scale_color_manual(values = my_colors) + 
  my_theme + 
  ggtitle("Trajectory by Subtype") +
  theme(legend.position = "none") 

ggsave("Trajectory_by_Subtype.pdf", plot = p_traj, width = 8, height = 6)




p_umap <- plot_cells(
  cds, 
  color_cells_by = "sub_celltype", 
  show_trajectory_graph = FALSE,       
  label_cell_groups = FALSE,           
  label_leaves = FALSE,
  label_branch_points = FALSE,
  cell_size = 0.8,
  group_label_size = 4                 
) + 
  scale_color_manual(values = my_colors) +
  my_theme + 
  ggtitle("Cell Subtypes") +
  
  guides(color = guide_legend(override.aes = list(size = 4))) 

ggsave("Cell_Subtypes.pdf", plot = p_umap, width = 10, height = 6)
