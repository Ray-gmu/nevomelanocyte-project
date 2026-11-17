
seurat <- NormalizeData(seurat, 
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000)


seurat <- FindVariableFeatures(seurat, 
                            selection.method = "vst", #mvp,disp
                            nfeatures = 2000) 

top10 <- head(VariableFeatures(seurat), 10)


plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, 
                     repel = TRUE);plot2
ggsave("Variablegenes.png",plot = plot2, width = 8,height = 5)


s.genes <- intersect(cc.genes$s.genes,rownames(seurat))
s.genes
g2m.genes <- intersect(cc.genes$g2m.genes,rownames(seurat))
g2m.genes
seurat <- CellCycleScoring(seurat,
                        s.features = s.genes,
                        g2m.features = g2m.genes,
                        set.ident = T)

seurat <- ScaleData(seurat,
                 # vars.to.regress = c("S.Score", "G2M.Score"),
                 features = rownames(seurat)) # features = NULL

# vars.to.regress去除某个因素的影响
# sce <- ScaleData(sce, vars.to.regress = "percent.mt")#去除线粒体影响
# sce <- ScaleData(sce, vars.to.regress = c("S.Score", "G2M.Score"))#去除周期影响

# counts <- as.matrix(GetAssayData(sce,slot = "counts"));counts
# data <- as.matrix(GetAssayData(sce,slot = "data"));data
# scale.data <- as.matrix(GetAssayData(sce,slot = "scale.data"));scale.data

qsave(seurat,file = "seurat_scaleData.qs")
seurat <- qread("sce_scaleData.qs")
