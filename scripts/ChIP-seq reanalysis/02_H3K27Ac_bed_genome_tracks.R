library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
library(org.Hs.eg.db)
library(GenomicRanges)
library(dplyr)
bed_file <- "~/Downloads/Human_Orgs_NPCs_H3K27Ac/Replicated/ENCFF011BMV_organoids.bed"
#or "~/Downloads/Human_Orgs_NPCs_H3K27Ac/Replicated/ENCFF285PPN_NPC.bed"

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene



peakAnno <- annotatePeak(bed_file, 
                         tssRegion = c(-3000, 3000),
                         TxDb = txdb,
                         annoDb = "org.Hs.eg.db")
annotated_df <- as.data.frame(peakAnno)


colnames(annotated_df)[7] <- "score"

annotated_df$SYMBOL <- mapIds(org.Hs.eg.db,
                                   keys = annotated_df$geneId,
                                   column = "SYMBOL",
                                   keytype = "ENTREZID",
                                   multiVals = "first")

annotated_df_score <- annotated_df[, c("SYMBOL", "score")]

annotated_df_max <- annotated_df_score %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
  ungroup()

write.table(annotated_df, "~/Downloads/Human_Orgs_NPCs_H3K27Ac/orgs_annotated_df.tab",
            row.names = TRUE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")



