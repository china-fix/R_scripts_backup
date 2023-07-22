## set the working directory
setwd("/home/fix/disk_2/R_working/GPE_KEGG_GO")

# KEGG enrichment analysis
library("clusterProfiler")

# KEGG pathway over-representation analysis
## give a vector of genes (gene)

gene <- scan(file = "input_list", what = "list")

kk <- enrichKEGG(gene         = gene,
                 organism = "set",
                 pvalueCutoff = 1,
                 keyType = "kegg",
                 pAdjustMethod = 'BH')

head(kk)
install.packages("./org.SEnteritidis.eg.db", repos=NULL)
library(org.SEnteritidis.eg.db)
kk_readable <- setReadable(kk, OrgDb = org.SEnteritidis.eg.db, keyType = "GID")
head(kk_readable)

kk_df <- as.data.frame(kk_readable)
kk_df <- as.data.frame(kk)
write.table(kk_df, file = "out.tab", sep = "\t", row.names = FALSE)

# KEGG module over-representation analysis

mkk <- enrichMKEGG(gene = gene,
                   organism = 'set',
                   keyType = "kegg",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   pAdjustMethod = 'fdr')

head(mkk)
mkk_readable <- setReadable(mkk, OrgDb = org.SEnteritidis.eg.db, keyType = "GID")
head(mkk_readable)

mkk_df <- as.data.frame(mkk_readable)
write.table(mkk_df, file = "out.tab", sep = "\t", row.names = FALSE)

# Visualize enriched KEGG pathways

browseKEGG(kk, 'set00470')

library("pathview")
hsa04110 <- pathview(gene.data  = gene,
                     pathway.id = "set02020",
                     species    = "set")


# Visualization of functional enrichment result

## Bar Plot
library(enrichplot)
barplot(kk, showCategory=10) + ggtitle("barplot for P_G core-SNP")

## Other variables that derived using mutate can also be used as bar height or color as demonstrated in Figure 15.1B.
mutate(kk, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")

## Dot plot
dotplot(kk, showCategory=5) + ggtitle("SGa&SPu vs. Sent (CSNP)")

## Gene-Concept Network
cnetplot(kk, layout = 'kk',
         node_label='category',
         cex_label_category = 0.5,
         categorySize='p.adjust',
         showCategory =30,
         circular= FALSE,
         colorEdge=TRUE)

## Heatmap-like functional classification
heatplot(kk, showCategory=30)

## UpSet Plot
## The upsetplot is an alternative to cnetplot for visualizing the complex association between genes and gene sets. 
library("enrichplot")
upsetplot(kk)

## pubmed trend of enriched terms
terms <- kk$Description[1:5]
p <- pmcplot(terms, 2010:2020)
p2 <- pmcplot(terms, 2010:2020, proportion=FALSE)
plot_grid(p, p2, ncol=2)




