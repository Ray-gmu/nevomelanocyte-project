
seurat <- RunPCA(seurat,features = VariableFeatures(seurat),verbose = FALSE)
DimPlot(seurat, reduction = "pca", group.by = "orig.ident")
ggsave("pc_DimPlot.png",width = 8,height = 6)


seurat@reductions$pca@cell.embeddings[1:4,1:4]


p <- VizDimLoadings(seurat, dims = 1:4, reduction = "pca")
ggsave("pc_VizDimLoadings.png",width = 8,height = 16)


png("pc_DimHeatmap.png", width = 2100, height = 1800, res = 300)
p <- DimHeatmap(seurat, dims = c(1:3,6:8,12:14), cells = 500, balanced = TRUE, combine = TRUE)
dev.off()


ElbowPlot(seurat,ndims=50)
ggsave("ElbowPlot.png",width = 8,height = 6)


library(findPC)
stdev_sorted <- sort(seurat@reductions$pca@stdev, decreasing = TRUE) # 主成分降序排列
stdev_sorted


findPC(sdev = stdev_sorted, 
       number = 50, 
       method = "all", 
       figure = T)


seurat <- RunHarmony(seurat,group.by.vars ="orig.ident")  




dims_N <- "14" 
seurat <- FindNeighbors(seurat, dims = 1:dims_N, 
                     reduction = "harmony", #单样本需要pca
                     verbose = FALSE) 


seurat <- FindClusters(seurat, 
                    resolution = seq(0.1, 2, by = 0.1), 
                    verbose = FALSE)
head(Idents(seurat), 5)


library(clustree)
library(patchwork)
p1 <- clustree(seurat, prefix = 'RNA_snn_res.') + coord_flip();p1
ggsave("clustertree.png",plot = p1, width = 20,height = 12)


p2 <- DimPlot(seurat, group.by = 'RNA_snn_res.0.8', label = T);p2 
Tree <- p1 + p2 + plot_layout(widths = c(3, 1));print(Tree)
ggsave("clustertree_p2.png",plot = Tree, width = 20,height = 12)




seurat$seurat_clusters <- seurat$RNA_snn_res.0.8


seurat <- RunUMAP(seurat, dims = 1:dims_N,reduction = "harmony",
               umap.method = "uwot", metric = "cosine")
umap <- DimPlot(seurat,
                reduction = 'umap',
                label = T,
                group.by = "seurat_clusters");umap
ggsave("UMAP.pdf",plot = umap,width = 7,height = 6)


seurat <- RunTSNE(seurat,dims = 1:dims_N,reduction = "harmony")
tsne <- DimPlot(seurat,
                label = T,
                reduction = 'tsne',
                pt.size =1,
                group.by = "seurat_clusters");tsne
ggsave("TSNE.pdf",plot = tsne,width = 7,height = 6)





p1 <- DimPlot(seurat,label = T,group.by = "orig.ident")
p2 <- DimPlot(seurat,label = T,group.by = "scDblFinder.class")
p1|p2
ggsave("double cell.png",plot = p1|p2,width = 12,height = 6,dpi = 300)

seurat_singlet <- subset(seurat, subset = scDblFinder.class == "singlet")
p3 <- DimPlot(seurat_singlet,label = T,group.by = "scDblFinder.class")
ggsave("singlet cell.png",plot = p3,width = 6,height = 6,dpi = 300)

qsave(seurat_singlet,"seurat_singlet.qs")
qsave(seurat,"seurat_singlet_double.qs")
