library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)


markers <- read.csv("markers.RNA.csv")
table(markers$cluster)

filter_markers <- markers %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 0.25)


table(filter_markers$cluster)



ids <- bitr(filter_markers$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

any(duplicated(ids$ENTREZID))
any(duplicated(ids$SYMBOL))

merge_markers <- inner_join(filter_markers, ids, by = c("gene" = "SYMBOL"))



target_cluster <- "4"
enrich_genes <- merge_markers %>% 
  filter(cluster == target_cluster) %>% 
  pull(ENTREZID) %>% 
  unique()

## GO----
ego <- enrichGO(
  gene = enrich_genes,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",  
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  # universe = bg_ids    
)


p <- dotplot(ego, 
            showCategory=10,
            orderBy="GeneRatio", 
            label_format=130, 
            split="ONTOLOGY",
            title = paste("Cluster", target_cluster ,"Enrichment of GO")
            ) + 
  facet_grid(ONTOLOGY~., scale='free')+
  scale_color_viridis_c(option = "D", direction = -1, name = "p.adjust")



print(p)










ggsave(paste("Clusters", target_cluster, "dotplot go.pdf"), 
       plot = p, width = 10, height = 8)



ego_df <- as.data.frame(ego) 


ego_cluster_readable <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType="ENTREZID")


ego_df$geneName <- as.data.frame(ego_cluster_readable)$geneID


write.csv(ego_df, 
          file = paste("Cluster_", target_cluster, "_GO_results.csv"), 
          row.names = FALSE)

## KEGG ----
ekegg <- enrichKEGG(
  gene = enrich_genes,  
  organism = "hsa",        
  pvalueCutoff = 0.05,
  # universe = bg_ids    
)


p <- dotplot(ekegg, 
        showCategory=20,
        orderBy="GeneRatio", 
        label_format=130, 
        title = paste("Cluster", target_cluster ,"Enrichment of KEGG"))
ggsave(paste("Clusters", target_cluster, "dotplot kegg.pdf"), 
       plot = p, width = 10, height = 8)


ekegg_df <- as.data.frame(ekegg)


ekegg_cluster_readable <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")


ekegg_df$geneName <- as.data.frame(ekegg_cluster_readable)$geneID


write.csv(ekegg_df, 
          file = paste("Cluster_", target_cluster, "_KEGG_results.csv"), 
          row.names = FALSE)

save(ego,ekegg,ego_df,ekegg_df, file = "enrich_genes.Rdata")
