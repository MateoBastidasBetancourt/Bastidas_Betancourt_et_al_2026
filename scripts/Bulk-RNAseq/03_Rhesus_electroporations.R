library(stringr)
library(dplyr)
library(tidyr)
library(edgeR)
library(limma)
library(EnsDb.Hsapiens.v79)
library(openxlsx)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(edgeR)
library(RColorBrewer)
library(stringr)

pheno_group$gene <- factor(pheno_group$gene, levels = c("Empty", "OVOL2", "ZNF90"))
pheno_group$batch <- factor(pheno_group$batch)

design <- model.matrix(~ 0 + gene + batch, data = pheno_group)
colnames(design)


design <- model.matrix(~ 0 + gene + batch, data = pheno_group)
dge_obj <- DGEList(counts=dataset_rhesus[,1:9], group=pheno_group[,3], genes=rownames(dataset_rhesus))
identical(pheno_group$sample.id, colnames(dataset_rhesus)[1:9]) #Se if order in count table is identical than in the metadata
dge_TMM <- calcNormFactors(dge_obj, method="TMM")
keep <- filterByExpr(dge_TMM, design = design)
dge_TMM <- dge_TMM[keep, keep.lib.sizes=FALSE] ##Filter low count reads, recalculate lib sizes after filtering

plotMDS(dge_TMM)
disp <- estimateDisp(dge_TMM, design, robust=TRUE)
plotBCV(disp)

fit <- glmQLFit(disp, design, robust=TRUE)


CONTRASTS <- makeContrasts(ZNF90vsEmpty = geneZNF90 -geneEmpty, 
                           OVOL2vsEmpty = geneOVOL2 -geneEmpty, 
                           ZNF90vsOVOL2 = geneZNF90 - geneOVOL2 ,
                           EmptyvsZNF90_OVOL2 = geneEmpty - (geneZNF90+geneOVOL2)/2, 
                           Empty_ZNF90vsOVOL2 = (geneEmpty+geneZNF90)/2 - geneOVOL2, 
                           Empty_OVOL2vsZNF90 = (geneEmpty+geneOVOL2)/2 - geneZNF90, levels=design)

for (i in 1:ncol(CONTRASTS)){
  
  contrast_name <- colnames(CONTRASTS)[i]
  
  # Perform the test
  current.glmQLFTest <- glmQLFTest(fit, contrast = CONTRASTS[, i])
  
  # Assign the result to a variable with the contrast name
  assign(contrast_name, current.glmQLFTest)
}


#####PVALUE CORRECTION AND SORTING GENES ACCORDING TO THEIR PVALUE
result_ZNF90vsEmpty <- topTags(ZNF90vsEmpty,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_OVOL2vsEmpty <- topTags(OVOL2vsEmpty,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_ZNF90vsOVOL2 <- topTags(ZNF90vsOVOL2,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_EmptyvsZNF90_OVOL2 <- topTags(EmptyvsZNF90_OVOL2,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_Empty_ZNF90vsOVOL2 <- topTags(Empty_ZNF90vsOVOL2,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_Empty_OVOL2vsZNF90 <- topTags(Empty_OVOL2vsZNF90,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")



######################
filter_results <- function(result) {
  filtered_result <- result$table[
    (result$table$logFC > 0.5 | result$table$logFC < -0.5) & result$table$FDR < 0.1, ]##$FDR < 0.1 for ZNF90vsEmpty
  return(filtered_result)
}
########################3


# Apply the function to each topTags object
result_ZNF90vsEmpty <- filter_results(result_ZNF90vsEmpty)
result_OVOL2vsEmpty <- filter_results(result_OVOL2vsEmpty)
result_ZNF90vsOVOL2 <- filter_results(result_ZNF90vsOVOL2)
result_EmptyvsZNF90_OVOL2 <- filter_results(result_EmptyvsZNF90_OVOL2)
result_Empty_ZNF90vsOVOL2 <- filter_results(result_Empty_ZNF90vsOVOL2)
result_Empty_OVOL2vsZNF90 <- filter_results(result_Empty_OVOL2vsZNF90)


#####################################

#Transcripts to gene_symbols - Human, more reliable, do not use BioMart:


edb <- EnsDb.Hsapiens.v79

ids_to_symbols <- function(deg_df){
deg_df$ensembl_base <- sub("\\..*", "", deg_df$genes)
id_map <- AnnotationDbi::select(edb,
                                keys = deg_df$ensembl_base,
                                keytype = "GENEID",
                                columns = c("GENEID", "SYMBOL"))

# Merge
deg_annotated <- merge(deg_df, id_map, by.x = "ensembl_base", by.y = "GENEID", all.x = TRUE)
deg_annotated$SYMBOL[is.na(deg_annotated$SYMBOL)] <- deg_annotated$ensembl_base[is.na(deg_annotated$SYMBOL)]

return(deg_annotated)
}

result_ZNF90vsEmpty <- ids_to_symbols(result_ZNF90vsEmpty)
result_OVOL2vsEmpty <- ids_to_symbols(result_OVOL2vsEmpty)
result_ZNF90vsOVOL2 <- ids_to_symbols(result_ZNF90vsOVOL2)
result_EmptyvsZNF90_OVOL2 <- ids_to_symbols(result_EmptyvsZNF90_OVOL2)
result_Empty_ZNF90vsOVOL2 <- ids_to_symbols(result_Empty_ZNF90vsOVOL2)
result_Empty_OVOL2vsZNF90 <- ids_to_symbols(result_Empty_OVOL2vsZNF90)

###############3333 plotting candidates
CAM <- c("CLDN3", "TJP3", "CLDN4", "EPCAM", "CDH1", "F11R", "JUP")

# COL genes (Collagens)
COL <- c("COL16A1", "COL20A1", "COL22A1", "COL14A1", "COL11A2")

# Neuro genes (Neuronal transcription factors)
Neuro <- c("NEUROD1", "NEUROG2", "HEY1", "ROBO3", "GLUD1", "STMN2", "HES6")

# ECM genes (Extracellular Matrix)
ECM <- c("THBS2", "THBS3", "MMP1", "MMP14")

PAT <- c("PAX7",
         "DBX1",
         "ALX4",
         "FOXC1",
         "SALL2",
         "SOX12",
         "FRZB",
         "CHRD",
         "GAS1", "GAD1", "NKX2-2")

genes_of_interest <- c(
  CAM,
  COL,
  Neuro,
  ECM,
  PAT
)
#gene_list <- XXXXX all genes in my analysis

# Merge back in the same order as your count matrix
gene_length_rhesus_liftoff <- read.delim("~/Downloads/gene_length_rhesus_liftoff.tabular", header=FALSE)
gene_length_rhesus_liftoff$V1 <- sub("^gene-", "", gene_length_rhesus_liftoff$V1)
colnames(gene_length_rhesus_liftoff) <- c("external_gene_name","gene_length")

# Convert lengths from bp → kb
ids_to_symbols <- function(tpm_mat, edb) {
  
  df <- as.data.frame(tpm_mat)
  df$ensembl_id <- rownames(df)
  df$ensembl_base <- sub("\\..*", "", df$ensembl_id)
  
  id_map <- AnnotationDbi::select(
    edb,
    keys = unique(df$ensembl_base),
    keytype = "GENEID",
    columns = c("GENEID", "SYMBOL")
  )
  
  # make lookup vector (THIS avoids merge issues completely)
  symbol_map <- setNames(id_map$SYMBOL, id_map$GENEID)
  
  df$SYMBOL <- symbol_map[df$ensembl_base]
  
  # fallback for missing
  df$SYMBOL[is.na(df$SYMBOL)] <- df$ensembl_base[is.na(df$SYMBOL)]
  
  df
}
library(EnsDb.Hsapiens.v79)
edb <- EnsDb.Hsapiens.v79

dataset_rhesus_mapped <- ids_to_symbols(dataset_rhesus,edb)
rownames(dataset_rhesus_mapped) <- make.unique(dataset_rhesus_mapped$SYMBOL)

gene_order <- rownames(dataset_rhesus_mapped)

idx <- match(gene_order, gene_length_rhesus_liftoff$external_gene_name)

keep <- !is.na(idx)

dataset_rhesus_matched <- dataset_rhesus_mapped[keep, ]

gene_length_ordered <- gene_length_rhesus_liftoff[idx[keep], ]

rownames(gene_length_ordered) <- rownames(dataset_rhesus_matched)

########
meta_cols <- c("ensembl_id", "ensembl_base", "SYMBOL", "SYMBOL_unique")
expr_cols <- setdiff(colnames(dataset_rhesus_matched), meta_cols)

counts_mat <- as.matrix(dataset_rhesus_matched[, expr_cols])
mode(counts_mat) <- "numeric"
#####
gene_length_kb <- gene_length_ordered$gene_length / 1000

# Step 1: calculate reads per kilobase (RPK)
rpk <- counts_mat / gene_length_kb

# Step 2: calculate per-sample scaling factors
scaling_factors <- colSums(rpk)

# Step 3: calculate TPM
tpm <- sweep(rpk, 2, scaling_factors, FUN = "/") * 1e6

tpm_log <- log1p(tpm)   # log(TPM + 1)
tpm_log <- tpm_log[genes_of_interest, ]
mat_z <- t(scale(t(tpm_log)))

gene_group <- ifelse(rownames(mat_z) %in% CAM, "CAM",
                     ifelse(rownames(mat_z) %in% COL, "ECM", 
                            ifelse(rownames(mat_z) %in% ECM, "ECM",
                                    ifelse(rownames(mat_z) %in% Neuro, "Neuro",
                                          ifelse(rownames(mat_z) %in% PAT, "Patterning", "Other")))))

ha <- rowAnnotation(
  Group = gene_group,
  col = list(
    Group = c(
      CAM = "#1b9e77",
      Neuro = "#7570b3",
      ECM = "#e7298a",
      Patterning = "#66a61e",
      Other = "grey80"
    )
  ),
  annotation_legend_param = list(
    Group = list(
      title_gp = gpar(fontsize = 16, fontface = "bold"),
      labels_gp = gpar(fontsize = 14)
    )
  )
)

colnames(mat_z) <- c("Empty-1","OVOL2-1", "ZNF90-1", "Empty-2", "OVOL2-2", "ZNF90-2", "Empty-3" ,"OVOL2-3", "ZNF90-3")
# col_fun <- colorRamp2(
#   c(-2, 0, 2),
#   c("#2166ac", "white", "#b2182b")
# )
# lgd <- Legend(col_fun = col_fun, title = "Z-score")
# ComplexHeatmap::width(ldg)
Heatmap(
  mat_z[,c(3,6,9,2,5,8)],
  name = "Z-score",
  left_annotation = ha,
  show_row_names = TRUE,
  show_column_names = TRUE,
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), cluster_rows = F, cluster_columns = F,
  column_names_gp = gpar(fontsize = 16), heatmap_legend_param = list(
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 14)
  )
)

### group samples

# meta_cols <- c("ensembl_id", "ensembl_base", "SYMBOL", "SYMBOL_unique")
# expr_cols <- setdiff(colnames(tpm_log), meta_cols)
# 
#expr_mat <- as.matrix(tpm_log[, expr_cols])
mode(tpm_log) <- "numeric"

grp1 <- seq(1, ncol(tpm_log), by = 3)
grp2 <- seq(2, ncol(tpm_log), by = 3)
grp3 <- seq(3, ncol(tpm_log), by = 3)
##average and SD:
out <- data.frame(
  gene = rownames(tpm_log),
  
  mean_Empty = rowMeans(tpm_log[, grp1], na.rm = TRUE),
  sd_grp1   = apply(tpm_log[, grp1], 1, sd, na.rm = TRUE),
  
  mean_OVOL2 = rowMeans(tpm_log[, grp2], na.rm = TRUE),
  sd_grp2   = apply(tpm_log[, grp2], 1, sd, na.rm = TRUE),
  
  mean_ZNF90 = rowMeans(tpm_log[, grp3], na.rm = TRUE),
  sd_grp3   = apply(tpm_log[, grp3], 1, sd, na.rm = TRUE)
)

write.table(out,"candidates_ZNF90_OVOL2_electroporations_plots.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

mat <- out[, c("mean_OVOL2", "mean_ZNF90")]
rownames(mat) <- out$gene
mat <- as.matrix(mat)

mat <- log2(mat + 1)
mat_z <- t(scale(t(mat)))# 
gene_group <- ifelse(rownames(mat_z) %in% CAM, "CAM",
                     ifelse(rownames(mat_z) %in% COL, "COL",
                            ifelse(rownames(mat_z) %in% Neuro, "Neuro",
                                   ifelse(rownames(mat_z) %in% ECM, "ECM",
                                          ifelse(rownames(mat_z) %in% PAT, "PAT", "Other")))))

# 
# 
ha <- rowAnnotation(
  Group = gene_group,
  col = list(
    Group = c(
      CAM = "#1b9e77",
      COL = "#d95f02",
      Neuro = "#7570b3",
      ECM = "#e7298a",
      PAT = "#66a61e",
      Other = "grey80"
    )
  )
)



diff_mat <- mat[, "mean_ZNF90"] - mat[, "mean_OVOL2"]

diff_mat <- matrix(diff_mat, ncol = 1)
rownames(diff_mat) <- rownames(mat)
colnames(diff_mat) <- "ZNF90 - OVOL2"

col_fun_diff <- colorRamp2(
  c(min(diff_mat), 0, max(diff_mat)),
  c("#2166ac", "white", "#b2182b")
)

Heatmap(
  diff_mat,
  name = "Difference",
  col = col_fun_diff,
  left_annotation = ha,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE
)



gene_sets <- list(
  CAM = CAM,
  COL = COL,
  Neuro = Neuro,
  ECM = ECM,
  PAT = PAT
)

plot_gene_set <- function(geneset_name, genes) {
  
  df <- out %>%
    filter(gene %in% genes) %>%
    select(gene, mean_OVOL2, sd_grp2, mean_ZNF90, sd_grp3) %>%
    mutate(gene = factor(gene, levels = genes)) %>%
    pivot_longer(
      cols = c(mean_OVOL2, mean_ZNF90),
      names_to = "condition",
      values_to = "mean_expr"
    ) %>%
    mutate(
      sd_expr = ifelse(condition == "mean_OVOL2", sd_grp2, sd_grp3),
      condition = recode(condition,
                         mean_OVOL2 = "OVOL2",
                         mean_ZNF90 = "ZNF90")
    )
  
  ggplot(df, aes(x = gene, y = mean_expr, fill = condition)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(
      aes(ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
      position = position_dodge(width = 0.8),
      width = 0.25
    ) +
    theme_classic(base_size = 14) +
    labs(
      title = geneset_name,
      x = NULL,
      y = "Mean log1p(TPM)",
      fill = "Condition"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

plot_gene_set("CAM genes", CAM)
plot_gene_set("COL genes", COL)
plot_gene_set("Neuro genes", Neuro)
plot_gene_set("ECM genes", ECM)
plot_gene_set("PAT genes", PAT)


dfs <- list(result_ZNF90vsEmpty = result_ZNF90vsEmpty, result_OVOL2vsEmpty = result_OVOL2vsEmpty,
            result_ZNF90vsOVOL2 = result_ZNF90vsOVOL2, result_EmptyvsZNF90_OVOL2 = result_EmptyvsZNF90_OVOL2, result_Empty_ZNF90vsOVOL2=result_Empty_ZNF90vsOVOL2, 
            result_Empty_OVOL2vsZNF90=result_Empty_OVOL2vsZNF90)

# Write them to an Excel file
write.xlsx(dfs, file = "rhesus_electroporation_degs_260526.xlsx")
