#scDblFinder
library(scDblFinder)
library(SingleCellExperiment)

sc_data <- as.SingleCellExperiment(seurat_qc)
sc_data <- scDblFinder(sc_data, 
                       samples="orig.ident",
                    )
seurat_dc <- CreateSeuratObject(
  counts = assay(sc_data, "counts"),  
  assay = "RNA",                  
  meta.data = as.data.frame(colData(sc_data))  
)

table(seurat_dc$scDblFinder.class) 


qsave(seurat_dc,file = "seurat_dc.qs")

seurat_dc <- qread(file = "seurat_dc.qs")
