
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)
library(harmony)
library(qs)
library(scDblFinder)
library(BiocParallel) 



seurat_all <- qread("seurat_all.qs")


table(seurat_all@meta.data$orig.ident)



seurat_all[["percent.mt"]] <- PercentageFeatureSet(seurat_all, pattern = "^MT-")
seurat_all[["percent.rp"]] <- PercentageFeatureSet(seurat_all, pattern = "^RP[SL]") 
seurat_all[["percent.hb"]] <- PercentageFeatureSet(seurat_all, pattern = "^HB[^(P)]")

p1 <- VlnPlot(seurat_all,
              features = c("nFeature_RNA",
                           "nCount_RNA",
                           "percent.mt",
                           "percent.rp",
                           "percent.hb"),
              ncol = 3,pt.size = 0, group.by = "orig.ident")
w=length(unique(seurat_all$orig.ident))/3+5;w
ggsave(filename="Vlnplot_raw.png",plot=p1,width = w,height = 5)


seurat_qc <- subset(seurat_all,
              percent.mt < 15 &  
                nFeature_RNA > 1000 & nFeature_RNA < 9000 & 
                nCount_RNA < 60000 & nCount_RNA > 1000  &
                #percent.rp <15 &
                percent.hb <1
              )




table(seurat_qc$orig.ident)

p2 <- VlnPlot(seurat_qc,
              features = c("nFeature_RNA",
                           "nCount_RNA",
                           "percent.mt",
                           "percent.rp",
                           "percent.hb"),
              ncol = 3,pt.size = 0, group.by = "orig.ident")
w=length(unique(seurat_qc$orig.ident))/2+5;w
ggsave("VlnPlot_afterqc.png",plot=p2,width = w,height = 5)


plot1 <- FeatureScatter(seurat_qc, feature1 = "nCount_RNA",
                        feature2 = "percent.mt",
                        group.by = "orig.ident")


plot2 <- FeatureScatter(seurat_qc, feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA",
                        group.by = "orig.ident")

plot3 <- FeatureScatter(seurat_qc, feature1 = "nCount_RNA",
                        feature2 = "percent.rp",
                        group.by = "orig.ident")

plot4 <- FeatureScatter(seurat_qc, feature1 = "percent.rp",
                        feature2 = "percent.mt",
                        group.by = "orig.ident")

plot5 <- FeatureScatter(seurat_qc, feature1 = "nCount_RNA", 
                        feature2 = "percent.hb", 
                        group.by = "orig.ident")
p3 <- plot1 + plot2 + plot3 + plot4 + plot5
w=length(unique(seurat_qc$orig.ident))/2+8;w
ggsave("FeatureScatter_afterqc.png",plot = p3, width = w,height = 8)


qsave(seurat_qc,file = "seurat_qc.qs")

seurat_qc <- qread(file = "seurat_qc.qs")
