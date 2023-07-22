if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.EcK12.eg.db")
library("org.EcK12.eg.db")
keytypes(org.EcK12.eg.db)


library("clusterProfiler")
gene <- scan(file = "ECOL_list", what = "list")

ggo <- groupGO(gene     = gene,
               OrgDb    = org.EcK12.eg.db,
               ont      = "BP",
               level    = 2,
               readable = TRUE,
               keyType = "SYMBOL")

head(ggo)

## GO over-representation analysis
ego <- enrichGO(gene          = gene,
                OrgDb         = org.EcK12.eg.db,
                ont           = "MF",
                pool = TRUE,
                pAdjustMethod = "fdr",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE,
                keyType = "SYMBOL")
head(ego)
ego_df <- as.data.frame(ego)
write.table(ego_df, file = "out.tab", sep = "\t", row.names = FALSE)

## Visualize enriched GO terms as a directed acyclic graph
enrichplot::goplot(ego)

## Bar Plot
library(enrichplot)
barplot(ego, showCategory=10) 

## Dot plot
dotplot(ego, showCategory=10) + ggtitle("dotplot for P_G MF")

## Gene-Concept Network
cnetplot(ego, layout = 'kk',
         node_label='category',
         cex_label_category = 0.5,
         categorySize='p.adjust',
         showCategory =10,
         circular= FALSE,
         colorEdge=TRUE)

## Heatmap-like functional classification
heatplot(ego, showCategory=10)

## UpSet Plot
## The upsetplot is an alternative to cnetplot for visualizing the complex association between genes and gene sets. 
library("enrichplot")
upsetplot(ego)

## pubmed trend of enriched terms
terms <- kk$Description[1:5]
p <- pmcplot(terms, 2010:2020)
p2 <- pmcplot(terms, 2010:2020, proportion=FALSE)
plot_grid(p, p2, ncol=2)
