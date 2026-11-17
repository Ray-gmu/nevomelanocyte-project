#基因集打分
####单个基因集合----
#选择基因集
table(Idents(subset))
Idents(subset) <- "seurat_clusters"
Idents(subset) <- "sub_celltype"

#颜色
my_cols <- c("#E76254","#EF8A47","#EAC435","#8DD3C7","#32B897" )

#                      ↓修改参数
marker_list <- list(proliferative)
#            ↓修改参数
score <- "cola"

subset <- AddModuleScore(subset,features = marker_list, name = score)

#打分基因集VlnPlot
p <- VlnPlot(subset, 
             features = paste0(score, "1"),
             pt.size=0, flip=T,
             cols = my_cols) + 
  geom_boxplot(width = 0.15, outlier.shape = NA)+
  labs(title = paste("Score of", score, "in melanocytes"))  #标题
ggsave(
  filename = paste0("VlnPlot.", score, ".pdf"),
  plot = p,
  width = 6, height = 5, units = "in")

#打分基因集FeaturePlot
p <- FeaturePlot(subset,
                 features = paste0(score, "1"),
                 reduction = "umap",min.cutoff = "0.2")+
  labs(title = paste("Score of", score, "in Melanocytes"))  #标题
ggsave(
  filename = paste0("FeaturePlot.", score, ".pdf"),
  plot = p,
  width = 6, height = 5, units = "in")

#### 多基因集批量处理----
proliferative <- c("MKI67","TOP2A","PCNA","UBE2C")
cola <- c("COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL6A1","COL6A2","COL6A3")
melanogenesis <- c("TYR","DCT","PMEL","MITF")
stress <- c("GADD45B","HSPB1","FOS","JUN")
SASP <- c("IL6","CXCL2","CXCL8","MMP1","MMP3","SERPINE1","GADD45B")
SASP <- c("IL6","CXCL1","CXCL2","CXCL3","CXCL8","MMP1","MMP3","SERPINE1","GADD45B")
OXPHOS_stress <- c("NDUFA1","COX4I1", "ATP5F1A","TOMM6")
cell_cycle <- c(
  "MZT1","BOD1","KNSTRN","BCL7C","FBXO5","ZWINT","KMT5A","PSRC1",
  "MAD2L1","RACGAP1","TACC3","KIF23","CDCA5","UHRF1","TRIP13","SPC25","PRC1",
  "AURKB","CDT1","TPX2","KIFC1","BIRC5","CCNB1","TOP2A","NUSAP1","UBE2C","CDK1"
)
transcription <- c("POLR2A","GTF2B","TAF7","BRD2","EPC1","RLF","TSPYL2",
                   "HEXIM1","NR4A2","IFRD1","AKIRIN1","RBM39","CCNL1",
                   "TCERG1L-AS1","AL118516.1")
ecm <- c("COL1A1","COL3A1","COL4A1","COL6A1","FN1","DCN","LUM","BGN","SPARC",
         "MMP1","MMP3","SERPINE1","TIMP1","TIMP2","LOX","PLOD1","PLOD2",
         "SFRP2","SFRP4","CCN1","CCN2","GREM1","POSTN","TNC","THBS1"
)

# 1. 定义所有基因集 (放入一个命名的 list)
gene_sets <- list(
  proliferative = proliferative,
  cola = cola,
  melanogenesis = melanogenesis,
  stress = stress,
  SASP = SASP,
  OXPHOS_stress = OXPHOS_stress,
  cell_cycle = cell_cycle,
  transcription = transcription,
  ecm = ecm
)

# 2. 定义统一颜色
my_cols <- c("#BEB8DC", "#82B0D2", "#8491B0", "#C1D37F", "#E098AE")

# 3. 开始 for 循环处理
for (score_name in names(gene_sets)) {
  
  message(paste0("正在处理基因集: ", score_name))
  
  # --- 准备参数 ---
  marker_list <- list(gene_sets[[score_name]])
  # 统一列名格式，例如 SASP1, ecm1
  current_feature <- paste0(score_name, "1")
  
  # --- 安全计算分数 ---
  # 运行前先删除可能存在的同名列，防止生成 SASP2, SASP3 等
  subset@meta.data[[current_feature]] <- NULL
  subset <- AddModuleScore(subset, features = marker_list, name = score_name)
  
  # --- 1. 绘制并保存 VlnPlot ---
  p_vln <- VlnPlot(subset, 
                   features = paste0(score_name, "1"),
                   pt.size=0, flip=T,
                   cols = my_cols) + 
    geom_boxplot(width = 0.15, outlier.shape = NA)+
    labs(title = paste("Score of", score_name, "in melanocytes"))  #标题

  ggsave(
    filename = paste0("VlnPlot.", score_name, ".pdf"),
    plot = p_vln,
    width = 6, height = 5
  )
  
  # --- 2. 绘制并保存 FeaturePlot ---
  p_feat <- FeaturePlot(subset,
                        features = paste0(score_name, "1"),
                        reduction = "umap",min.cutoff = "0.2")+
    labs(title = paste("Score of", score_name, "in Melanocytes"))  #标题
  
  ggsave(
    filename = paste0("FeaturePlot.", score_name, ".pdf"),
    plot = p_feat,
    width = 6, height = 5
  )
  
  message(paste0(score_name, " 分析完成并已保存 PDF"))
}
