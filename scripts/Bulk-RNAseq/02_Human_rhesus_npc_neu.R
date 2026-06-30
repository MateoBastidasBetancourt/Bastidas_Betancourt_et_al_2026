# Johnson: unstranded
# Fietz: unstranded
# Florio: unstranded
# ENCODE:stranded

# Luo: Unstranded
# Micali: reverse
# Zhao: unstranded


library(edgeR)
library(limma)
library(stringr)

design <- model.matrix(~ 0 + species_cells, data = pheno_group)
dge_obj <- DGEList(counts=dataset_merged, group=pheno_group[,4], genes=rownames(dataset_merged))
identical(pheno_group$sample.id, colnames(dataset_merged)) #See if order in count table is identical than in the metadata
dge_TMM <- calcNormFactors(dge_obj, method="TMM")
keep <- filterByExpr(dge_TMM, design = design)
dge_TMM <- dge_TMM[keep, keep.lib.sizes=FALSE] 


plotMDS(dge_TMM, method="bcv", col = c("blue", "red")[factor(pheno_group$species)], pch = 16, cex = 1.5,xlab = "MDS1", ylab = "MDS2")
legend(-0.5, 1.5, legend = c("Human","Rhesus"),
       col = c("blue", "red"), pch = 5, cex = 0.3, title = "Species")

library(ggplot2)
disp <- estimateDisp(dge_TMM, design, robust=TRUE)
#plotBCV(disp)

fit <- glmQLFit(disp, design, robust=TRUE)

CONTRASTS <- makeContrasts(NpcvsNeu = (species_cellshuman_npc-species_cellshuman_neu)-(species_cellsrhesus_npc-species_cellsrhesus_neu), 
                           NpcvsNpc = species_cellshuman_npc-species_cellsrhesus_npc, 
                           HumanvsRhesus = (species_cellshuman_npc+species_cellshuman_neu)/2-(species_cellsrhesus_npc+species_cellsrhesus_neu)/2,
                           NeuvsNeu = species_cellshuman_neu-species_cellsrhesus_neu , levels=design)

for (i in 1:ncol(CONTRASTS)){
  
  contrast_name <- colnames(CONTRASTS)[i]
  
  # Perform the test
  current.glmQLFTest <- glmQLFTest(fit, contrast = CONTRASTS[, i])
  
  # Assign the result to a variable with the contrast name
  assign(contrast_name, current.glmQLFTest)
}


#####PVALUE CORRECTION AND SORTING GENES ACCORDING TO THEIR PVALUE
result_NpcvsNeu <- topTags(NpcvsNeu,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_NpcvsNeu <- result_NpcvsNeu$table
result_NpcvsNeu$log_pvalue <- -log(result_NpcvsNeu$FDR) #volcano plot
write.table(result_NpcvsNeu[,c(1,2,7)], "up_human_npcs_vs_rhesus_volcano_plot.txt", row.names = F, quote = F, col.names = T, sep = "\t")
result_NpcvsNpc <- topTags(NpcvsNpc,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_HumanvsRhesus <- topTags(HumanvsRhesus,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result_NeuvsNeu <- topTags(NeuvsNeu,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")

######################
filter_results <- function(result) {
  filtered_result <- result$table[
    (result$table$logFC > 1 | result$table$logFC < -1) & result$table$FDR < 0.01, ]
  return(filtered_result)
}

result_NpcvsNeu <- filter_results(result_NpcvsNeu)
result_NpcvsNpc <- filter_results(result_NpcvsNpc)
result_HumanvsRhesus <- filter_results(result_HumanvsRhesus)
result_NeuvsNeu <- filter_results(result_NeuvsNeu)

###########333 WGCNA
library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr)      # provides the %>% operator
library(WGCNA)
library(edgeR)
library(igraph)

logCPM <- cpm(dge_TMM, log = TRUE, prior.count = 1)
logCPM <- t(logCPM)

#now to WGCNA, see if i just do it with the DEGs or not...
allowWGCNAThreads() 
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

sft = pickSoftThreshold(
  logCPM,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)


par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")
picked_power = 5
temp_cor <- cor       
cor <- WGCNA::cor
netwk <- blockwiseModules(logCPM,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor  

mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

MEs0 <- moduleEigengenes(logCPM, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
#module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name)
  )
# Start with the full list
all_treatments <- mME$treatment

# Group by your desired order
others   <- sort(all_treatments[!grepl("ENCFF|^CP|^Neu|NPC|NSCP|^RM", all_treatments)])
encff    <- sort(all_treatments[grepl("ENCFF", all_treatments) & !all_treatments %in% others])
cp       <- sort(all_treatments[grepl("^CP", all_treatments) & !all_treatments %in% c(others, encff)])
neu      <- sort(all_treatments[grepl("^Neu", all_treatments) & !all_treatments %in% c(others, encff, cp)])
npc_nscp <- sort(all_treatments[grepl("NPC|NSCP", all_treatments) & !all_treatments %in% c(others, encff, cp, neu)])
rm_samples <- sort(all_treatments[grepl("^RM", all_treatments) & !all_treatments %in% c(others, encff, cp, neu, npc_nscp)])

# Combine into final levels
# Make sure the custom levels vector has no duplicates
custom_levels <- unique(c(others, encff, cp, neu, npc_nscp, rm_samples))

# Apply to factor safely
mME <- mME %>%
  mutate(treatment = factor(treatment, levels = custom_levels))


# Plot
mME %>% 
  ggplot(aes(x = treatment, y = name, fill = value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)
  ) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill = "corr")

modules_of_interest <- c("royalblue","greenyellow")

genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = t(logCPM)[genes_of_interest$gene_id,]
# Only recalculate TOM for modules of interest (faster, although there's some online discussion if this will be slightly off)
TOM = TOMsimilarityFromExpr(t(expr_of_interest), power = picked_power)
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)


###### takes a bit of computing if too many genes...
edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

edge_list_filtered <- edge_list %>%
  filter(correlation > 0.20) %>% #0.2 if using both modules
  group_by(gene1) %>%
  slice_max(order_by = correlation, n = 10) %>%
  ungroup()
edge_list_filtered <- edge_list_filtered[edge_list_filtered$gene1 < edge_list_filtered$gene2, ]
edge_list_filtered <- edge_list_filtered[edge_list_filtered$gene1 != edge_list_filtered$gene2, ]
write.table(unique(c(edge_list_filtered$gene1, edge_list_filtered$gene2)), "WGCNA_genes.tsv", row.names = F, quote = F)


###quantitative estimate of hub genes:
gene_connectivity <- edge_list_filtered %>%
  # keep only non-NA correlations
  filter(!is.na(correlation)) %>%
  # gather connections for both genes
  group_by(gene1) %>%
  summarise(connectivity = sum(correlation, na.rm = TRUE)) %>%
  rename(gene = gene1)

gene_connectivity2 <- edge_list_filtered %>%
  filter(!is.na(correlation)) %>%
  group_by(gene2) %>%
  summarise(connectivity2 = sum(correlation, na.rm = TRUE)) %>%
  rename(gene = gene2)

gene_connectivity_total <- gene_connectivity %>%
  full_join(gene_connectivity2, by = "gene") %>%
  mutate(connectivity_total = rowSums(across(c(connectivity, connectivity2)), na.rm = TRUE)) %>%
  arrange(desc(connectivity_total))
top_hubs <- gene_connectivity_total %>% slice_max(connectivity_total, n = 50)

top_hubs
library(ggplot2)

ggplot(top_hubs, aes(x = reorder(gene, connectivity_total), y = connectivity_total)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = "Gene", y = "Connectivity", title = "Top Hub Genes")

##network plot
# Optional: remove pairs with correlation = 1 (likely duplicates)
edge_list_filtered <- edge_list_filtered[edge_list_filtered$correlation < 1, ]
g <- graph_from_data_frame(edge_list_filtered)

cell_cycle_genes <- c(
  "MCMBP", "CCNY", "PSRC1", "CDKN1B", "MCM4", "RMI2", "RAD23B", "PPP4R2", "DEK",
  "HNRNPU", "USP39", "SENP6", "SMAD5", "SUZ12"
)

# Color nodes
V(g)$color <- ifelse(V(g)$name %in% cell_cycle_genes, "red", "lightblue")

# Bigger nodes for highlighted genes
V(g)$size <- ifelse(V(g)$name %in% cell_cycle_genes, 10, 4)

# Show labels ONLY for selected genes
V(g)$label <- ifelse(V(g)$name %in% cell_cycle_genes, V(g)$name, NA)

# Bigger font for highlighted genes
V(g)$label.cex <- ifelse(V(g)$name %in% cell_cycle_genes, 1.2, 0)

# Bold labels for highlighted genes
V(g)$label.font <- ifelse(V(g)$name %in% cell_cycle_genes, 2, 1)

# Plot
png("igraph_plot.png", width = 2000, height = 2000, res = 150)
plot(
  g,
  vertex.label.color = "black",
  vertex.label.dist = 0,
  edge.width = 1,
  edge.color = "gray"
)
dev.off()

