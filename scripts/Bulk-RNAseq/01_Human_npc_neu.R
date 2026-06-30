library(edgeR)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v79)
library(ensembldb)
library(biomaRt)
library(GenomicRanges)
library(dplyr)
library(ggplot2)


design <- model.matrix(~ 0 + species_cells, data = pheno_group)
dge_obj <- DGEList(counts=dataset_human, group=pheno_group[,4], genes=rownames(dataset_human))
identical(pheno_group$sample.id, colnames(dataset_human)) #See if order in count table is identical than in the metadata
dge_TMM <- calcNormFactors(dge_obj, method="TMM")
keep <- filterByExpr(dge_TMM, design = design, min.prop = 0.7, min.count=3)
dge_TMM <- dge_TMM[keep, keep.lib.sizes=FALSE]

plotMDS(dge_TMM, method="bcv", col = c("blue", "red")[factor(pheno_group$species_cells)], pch = 16, cex = 1.5,xlab = "MDS1", ylab = "MDS2")
disp <- estimateDisp(dge_TMM, design, robust=TRUE)
plotBCV(disp)

fit <- glmQLFit(disp, design, robust=TRUE)
qlt.cells <- glmQLFTest(fit, contrast=makeContrasts(species_cellshuman_npc-species_cellshuman_neu, levels=design))

result <- topTags(qlt.cells,n=nrow(dge_TMM),adjust.method="BH",sort.by="PValue")
result <- result$table

library(org.Hs.eg.db)
library(ensembldb)
library(EnsDb.Hsapiens.v79)

ensembl_ids <- result$genes
gene_symbols <- mapIds(EnsDb.Hsapiens.v79, keys=ensembl_ids, 
                       column="GENENAME", keytype="GENEID", 
                       multiVals="first")
gene_symbols <- data.frame(gene_symbols)
###################################
# anno <- read.csv("biomart_anno_complete.txt", header = TRUE)
# anno <- unique(anno)
#####################################
result_2 <- merge(result, gene_symbols, by=0)

rownames(result_2)<- result_2[,2]
result_2 <- result_2[,-2]
colnames(result_2)[1] <- "Gene.name"
result_2[1] <- result_2[7]
result_2$log_pvalue <- -log(result_2$FDR) #volcano plot
write.table(result_2[,c(1,2,8)], "human_npcs_only_human_volcano_plot.txt", row.names = F, quote = F, col.names = T, sep = "\t")

result_final <- merge(result_2, dataset_human, by=0)

symbols <- AnnotationDbi::select(org.Hs.eg.db, keys=result_final[is.na(result_final$Gene.name),1], columns="SYMBOL", keytype="ENSEMBL")
symbols <- symbols[!duplicated(symbols$ENSEMBL),]
result_final[is.na(result_final$Gene.name),2] <- symbols$SYMBOL
result_final[is.na(result_final$Gene.name),2] <- result_final[is.na(result_final$Gene.name),1]

up_human <- subset(result_final, result_final$logFC>=1 & result_final$FDR<0.01)
down_human <- subset(result_final, result_final$logFC<=-1 & result_final$FDR<0.01)
non_de <- subset(result_final, 
                 FDR > 0.01 & logFC > -1 & logFC < 1)

anyNA(result_final$Gene.name)

#Coordinates - genes up, down and non-de
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

up_genes <- sub("\\..*", "", up_human$Row.names)
down_genes <- sub("\\..*", "", down_human$Row.names)
non_de_genes <- sub("\\..*", "", non_de$Row.names)

all_genes <- unique(c(up_genes, down_genes, non_de_genes))

gene_coords <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand"
  ),
  filters = "ensembl_gene_id",
  values = all_genes,
  mart = mart
)


gene_coords_clean <- gene_coords %>%
  rename(
    "Gene stable ID" = ensembl_gene_id,
    "Chromosome/scaffold name" = chromosome_name,
    "Gene start (bp)" = start_position,
    "Gene end (bp)" = end_position
  ) %>%
  filter(`Chromosome/scaffold name` %in% c(1:22, "X", "Y"))


coordinates_up <- gene_coords_clean %>%
  filter(`Gene stable ID` %in% up_genes)

coordinates_down <- gene_coords_clean %>%
  filter(`Gene stable ID` %in% down_genes)

coordinates_non_de <- gene_coords_clean %>%
  filter(`Gene stable ID` %in% non_de_genes)

write.table(coordinates_non_de, "coordinates_non_de.tsv", row.names = F, quote = F,  col.names = T, sep = "\t")
write.table(coordinates_up, "coordinates_up.tsv", row.names = F, quote = F,  col.names = T, sep = "\t")
write.table(coordinates_down, "coordinates_down.tsv", row.names = F, quote = F, col.names = T, sep = "\t")


##HAR analysis

setwd("/home/basti/Documents/PhD/Data/coordinates_human_vs_human/")
####INTRA-SPECIFIC APPROACH (HUMAN VS HUMAN)

#hars <- read.table("/home/basti/Documents/PhD/Data/zooHARs.txt", sep="\t", header=TRUE) [,]
# up_gene_loci <- read.table("coordinates_up.tsv", sep="\t", header=TRUE)
# down_gene_loci <- read.table("coordinates_down.tsv", sep="\t", header=TRUE)
# non_de_loci <- read.table("coordinates_non_de.tsv", sep="\t", header=TRUE)



hars$chrom <- sub("^chr", "", hars$chrom)

hars_gr <- GRanges(seqnames = hars$chrom, ranges = IRanges(start = hars$start, end = hars$end))
up_gene_loci_gr <- GRanges(seqnames = up_gene_loci$Chromosome.scaffold.name, ranges = IRanges(start = up_gene_loci$Gene.start..bp., end = up_gene_loci$Gene.end..bp.), strand = up_gene_loci$strand)
down_gene_loci_gr <- GRanges(seqnames = down_gene_loci$Chromosome.scaffold.name, ranges = IRanges(start = down_gene_loci$Gene.start..bp., end = down_gene_loci$Gene.end..bp.), strand = down_gene_loci$strand)

non_de_loci_gr <- GRanges(seqnames = non_de_loci$Chromosome.scaffold.name, ranges = IRanges(start = non_de_loci$Gene.start..bp., end = non_de_loci$Gene.end..bp.), strand = non_de_loci$strand)

distances_up <- mcols(distanceToNearest(up_gene_loci_gr, hars_gr))$distance
distances_down <- mcols(distanceToNearest(down_gene_loci_gr, hars_gr))$distance
distances_non_de<- mcols(distanceToNearest(non_de_loci_gr, hars_gr))$distance

distances_de <- c(distances_up, distances_down)


#### make a 3x2 contingency table for the three conditions, including down too:
# Set threshold
threshold <- 100000  # 100 kb

# Determine which genes have HARs within threshold
near_up <- distances_up < threshold
near_down <- distances_down < threshold
near_non_de <- distances_non_de < threshold

# Build 3x2 contingency table
table_data <- matrix(
  c(
    sum(near_up), sum(!near_up),
    sum(near_down), sum(!near_down),
    sum(near_non_de), sum(!near_non_de)
  ),
  nrow = 3,
  byrow = TRUE
)

# Add row and column names
rownames(table_data) <- c("up", "down", "non_DE")
colnames(table_data) <- c("HAR < 100kb", "HAR >= 100kb")

# Inspect table
table_data
chisq.test(table_data) #cool! u

prop_up <- mean(near_up) 
prop_down <- mean(near_down) 
prop_non_de <- mean(near_non_de) 
prop_up; prop_down; prop_non_de ####up_genes in npcs are driving the observed differences

##### see if hars are closer to upregulated genes than all of the rest

table_up_down <- matrix(c(
  sum(near_up), sum(!near_up),
  sum(near_down), sum(!near_down)
), nrow = 2)

fisher_up_down <- fisher.test(table_up_down)
pval_up_down <- fisher_up_down$p.value

### --------- Fisher test 2: UP vs NON-DE ---------
table_up_non_de <- matrix(c(
  sum(near_up), sum(!near_up),
  sum(near_non_de), sum(!near_non_de)
), nrow = 2)

fisher_up_non_de <- fisher.test(table_up_non_de)
pval_up_non_de <- fisher_up_non_de$p.value



### --------- Proportions ---------
plot_df <- data.frame(
  group = c("NPC-up", "NPC-down", "non-DE"),
  proportion = c(prop_up, prop_down, prop_non_de)
)

plot_df$se <- sqrt(plot_df$proportion * (1 - plot_df$proportion) / 
                     c(length(near_up), length(near_down), length(near_non_de)))

plot_df$percentage <- plot_df$proportion * 100
plot_df$se_perc <- plot_df$se * 100

### --------- P-value to stars ---------
p_to_stars <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")
}

stars_up_down <- p_to_stars(pval_up_down)
stars_up_non_de <- p_to_stars(pval_up_non_de)

### --------- Plot ---------
ggplot(plot_df, aes(x = group, y = percentage, fill = group)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(
    aes(ymin = percentage - se_perc, ymax = percentage + se_perc),
    width = 0.2
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(plot_df$percentage + plot_df$se_perc) * 1.15,
    label = stars_up_down,
    size = 8
  ) +
  annotate(
    "text",
    x = 2,
    y = max(plot_df$percentage + plot_df$se_perc) * 1.25,
    label = stars_up_non_de,
    size = 8
  ) +
  labs(
    y = "Genes with HAR < 100 kb (%)",
    x = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none")


library(VennDiagram)
library(grid)
load("~/Documents/PhD/Data/environment_rhesus_electroporations_260526.RData")
florio_genes <- unique(Florio_genes$V1)
our_genes <- unique(up_npc_human$Gene.name)
shared <- intersect(florio_genes, our_genes)
grid.newpage()
venn.plot <- draw.pairwise.venn(
  area1 = length(florio_genes),
  area2 = length(our_genes),
  cross.area = length(shared),
  category = c("", ""),
  fill = c("grey40", "grey75"),
  alpha = c(0.5, 0.5),
  lty = "blank",
  cex = 0,
  cat.cex = 0,
  scaled = TRUE
)



