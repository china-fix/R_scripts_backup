## make SE organism package
fSym <- read.csv('fSym.csv')
fChr <- read.csv('fChr.csv')
fGO <- read.csv('fGO.csv')

library('AnnotationForge')
makeOrgPackage(gene_info=fSym, chromosome=fChr, go=fGO,
               version="0.1",
               maintainer="Xiao Fei <xiao.fei@sund.ku.dk>",
               author="Xiao Fei <xiao.fei@sund.ku.dk>",
               outputDir = ".",
               tax_id="550537",
               genus="Salmonella",
               species="Enteritidis",
               goTable="go")

install.packages("./org.SEnteritidis.eg.db", repos=NULL)

library("org.SEnteritidis.eg.db")

# GO enrichment analysis
library("clusterProfiler")

## If an organism is not supported by AnnotationHub, user can use the AnnotationForge package to build OrgDb manually.

makeOrgPackageFromNCBI(version = "0.1",
                       author = "Xiao Fei <xiao.fei@sund.ku.dk>",
                       maintainer = "Xiao Fei <xiao.fei@sund.ku.dk>",
                       outputDir = ".",
                       tax_id = "550537",
                       genus = "Salmonella",
                       species = "Salmonella enterica")
## GO classification
gene <- scan(file = "input_list", what = "list")
library(org.SEnteritidis.eg.db)
ggo <- groupGO(gene     = gene,
               OrgDb    = org.SEnteritidis.eg.db,
               ont      = "BP",
               level    = 2,
               readable = TRUE,
               keyType = "GID")

head(ggo)

## GO over-representation analysis
ego <- enrichGO(gene          = gene,
                OrgDb         = org.SEnteritidis.eg.db,
                ont           = "ALL",
                pool = TRUE,
                pAdjustMethod = "fdr",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE,
                keyType = "GID")
head(ego)
ego_df <- as.data.frame(ego)
write.table(ego_df, file = "out.tab", sep = "\t", row.names = FALSE)

## filter the ego at a specific level


## Visualize enriched GO terms as a directed acyclic graph
enrichplot::goplot(ego)

## Bar Plot
library(enrichplot)
barplot(ego, showCategory=10) 
barplot(ego, split = "ONTOLOGY", label_format = 60)+facet_grid(ONTOLOGY~., scale = "free")+ggtitle("Barplot for GO")

## Dot plot
dotplot(ego, showCategory=15)
dotplot(ego, split = "ONTOLOGY", label_format = 60)+facet_grid(ONTOLOGY~., scale = "free")+ggtitle("Dotplot for GO")

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














