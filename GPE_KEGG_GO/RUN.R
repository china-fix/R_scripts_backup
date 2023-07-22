gene <- scan(file = "input_list", what = "list")

kk <- enrichKEGG(gene         = gene,
                 organism = "set",
                 pvalueCutoff = 1,
                 keyType = "kegg",
                 pAdjustMethod = 'fdr')

kk_readable <- setReadable(kk, OrgDb = org.SEnteritidis.eg.db, keyType = "GID")
kk_df <- as.data.frame(kk_readable)
write.table(kk_df, file = "out_KEGG.tab", sep = "\t", row.names = FALSE)

mkk <- enrichMKEGG(gene = gene,
                   organism = 'set',
                   keyType = "kegg",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   pAdjustMethod = 'fdr')

mkk_readable <- setReadable(mkk, OrgDb = org.SEnteritidis.eg.db, keyType = "GID")
mkk_df <- as.data.frame(mkk_readable)
write.table(mkk_df, file = "out_MKEGG.tab", sep = "\t", row.names = FALSE)

ego <- enrichGO(gene          = gene,
                OrgDb         = org.SEnteritidis.eg.db,
                ont           = "ALL",
                pool = TRUE,
                pAdjustMethod = "fdr",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE,
                keyType = "GID")
ego_df <- as.data.frame(ego)
write.table(ego_df, file = "out_GO.tab", sep = "\t", row.names = FALSE)
