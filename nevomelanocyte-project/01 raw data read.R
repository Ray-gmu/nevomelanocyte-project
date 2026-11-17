
library(stringr)
library(Seurat)
library(data.table)
library(data.table)
library(qs)

samples <- list.files("C:/Users/wangr/Desktop/痣原代2/primary cell raw data")
samples


seurat_list <- list()


for (sample in samples) {
                        
  data.path <- paste0("C:/Users/wangr/Desktop/痣原代2/primary cell raw data/", sample)
  

  seurat_data <- Read10X(data.dir = data.path)
  
  
  seurat_obj <- CreateSeuratObject(counts = seurat_data,project = sample,min.features = 300,min.cells = 5)
  
  
  seurat_list <- append(seurat_list, seurat_obj)
}


seurat_combined <- merge(seurat_list[[1]], 
                         y = seurat_list[-1],
                         add.cell.ids = samples)


seurat_all = JoinLayers(seurat_combined)


qsave(seurat_all,file = "seurat_all.qs")

seurat_all <- qread(file = "seurat_all.qs")