library(Seurat)
library(dplyr)
library(ggplot2)
library(forcats)
library(tidyr)


seurat_obj_combined <- readRDS("/mnt/vast-standard/home/cbastid/u15756/seurat_obj_combined.Rds")
### pseudobulk analysis of apical progenitors
pseudo_ifnb <- AggregateExpression(seurat_obj_combined, return.seurat = T, group.by = c("species", "cell_line", "harmony_clusters"))
aRG1_pseudo <- FindMarkers(pseudo_ifnb, ident.1 = c("1_Marmoset","10_Marmoset"), ident.2 = c("1_Human","10_Human"), logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1)
aRG1_pseudo$gene <- rownames(aRG1_pseudo)
aRG1_pseudo$p_val.adj <- p.adjust(aRG1_pseudo$p_val, method = "BH")
#G2M_aRG_pseudo <- G2M_aRG_pseudo[G2M_aRG_pseudo$avg_log2FC>1 | G2M_aRG_pseudo$avg_log2FC < -1,]
aRG1_pseudo <- aRG1_pseudo[aRG1_pseudo$p_val.adj<0.1,]


aRG2_pseudo <- FindMarkers(pseudo_ifnb, ident.1 = c("8_Marmoset","12_Marmoset"), ident.2 = c("8_Human","12_Human"), logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1)
aRG2_pseudo$gene <- rownames(aRG2_pseudo)
aRG2_pseudo$p_val.adj <- p.adjust(aRG2_pseudo$p_val, method = "BH")

aRG3_pseudo <- FindMarkers(pseudo_ifnb, ident.1 = c("6_Marmoset","9_Marmoset"), ident.2 = c("6_Human","9_Human"), logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1)
aRG3_pseudo$gene <- rownames(aRG3_pseudo)
aRG3_pseudo$p_val.adj <- p.adjust(aRG3_pseudo$p_val, method = "BH")


aRG4_pseudo <- FindMarkers(pseudo_ifnb, ident.1 = c("2_Marmoset"), ident.2 = c("2_Human"), logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1)
aRG4_pseudo$gene <- rownames(aRG4_pseudo)
aRG4_pseudo$p_val.adj <- p.adjust(aRG4_pseudo$p_val, method = "BH")

#####################
#### plots
data <- FetchData(
  seurat_obj_combined,
  vars = c("OVOL2", "ZNF90", "Phase", "species", "cell_line")
)
data$ident <- Idents(seurat_obj_combined)

human_OVOL2 <- data %>%
  filter(OVOL2 > 0) %>%
  mutate(Gene = "OVOL2")

human_OVOL2 <- human_OVOL2 %>% filter(species == "Human")

human_ZNF90 <- data %>%
  filter(ZNF90 > 0) %>%
  mutate(Gene = "ZNF90")

human_ZNF90 <- human_ZNF90 %>% filter(species == "Human")

combined <- bind_rows(
  human_OVOL2,
  human_ZNF90
)

plot_data <- combined %>%
  group_by(Gene, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Gene) %>%
  mutate(prop = n / sum(n))

ggplot(plot_data, aes(x = Gene, y = prop, fill = Phase)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  ylab("Fraction of expressing cells") +
  xlab("") +
  scale_y_continuous(labels = scales::percent)
########
#Now also plot cell cycle per cluster for cells co-expressing the two genes
data <- FetchData(
  seurat_obj_combined,
  vars = c("OVOL2", "ZNF90", "Phase", "species", "ident")
)
coexpressing <- data %>%
  filter(
    species == "Human",
    OVOL2 > 0,
    ZNF90 > 0
  )

plot_data <- coexpressing %>%
  group_by(ident, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(ident) %>%
  mutate(prop = n / sum(n))

cluster_order <- coexpressing %>%
  count(ident) %>%
  arrange(desc(n)) %>%
  pull(ident)

plot_data$ident <- factor(plot_data$ident, levels = cluster_order)

ggplot(plot_data, aes(x = ident, y = prop, fill = Phase)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  ylab("Fraction of OVOL2+ZNF90 cells") +
  xlab("Cluster") +
  scale_y_continuous(labels = scales::percent)

cluster_counts <- coexpressing %>% count(ident)
cluster_counts
###########################
total_cluster_cells <- data %>%
  filter(species == "Human") %>%
  count(ident)

ovol2_cells <- data %>%
  filter(
    species == "Human",
    OVOL2 > 0
  ) %>%
  count(ident)
znf90_cells <- data %>%
  filter(
    species == "Human",
    ZNF90 > 0
  ) %>%
  count(ident)

coexpressing <- data %>%
  filter(
    species == "Human",
    OVOL2 > 0,
    ZNF90 > 0
  )

all_clusters <- total_cluster_cells$ident
ovol2_cells <- ovol2_cells %>%
  tidyr::complete(ident = all_clusters, fill = list(n = 0))
coexp_cluster_cells <- coexp_cluster_cells %>%
  tidyr::complete(ident = all_clusters, fill = list(n = 0))
znf90_cells <- znf90_cells %>%
  tidyr::complete(ident = all_clusters, fill = list(n = 0))
ovol2_props <- ovol2_cells %>%
  left_join(total_cluster_cells, by = "ident") %>%
  mutate(
    prop = n.x / n.y,
    group = "OVOL2"
  )
znf90_props <- znf90_cells %>%
  left_join(total_cluster_cells, by = "ident") %>%
  rename(
    expr_cells = n.x,
    total_cells = n.y
  ) %>%
  mutate(
    prop = expr_cells / total_cells,
    group = "ZNF90"
  )
double_props <- coexp_cluster_cells %>%
  left_join(total_cluster_cells, by = "ident") %>%
  rename(
    expr_cells = n.x,
    total_cells = n.y
  ) %>%
  mutate(
    prop = expr_cells / total_cells,
    group = "Double positive"
  )
combined_props <- bind_rows(
  ovol2_props,
  znf90_props,
  double_props
)
cluster_order <- total_cluster_cells %>%
  left_join(double_props, by = "ident") %>%
  mutate(prop = ifelse(is.na(prop), 0, prop)) %>%
  arrange(desc(prop)) %>%
  pull(ident)

all_clusters <- total_cluster_cells$ident

remaining_clusters <- setdiff(all_clusters, cluster_order)

cluster_order <- c(cluster_order, remaining_clusters)
combined_props$ident <- factor(combined_props$ident, levels = cluster_order)
ggplot(combined_props,
       aes(x = ident,
           y = prop,
           fill = group)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  scale_fill_manual(values = c(
    "OVOL2" = "#D55E00",
    "ZNF90" = "#0072B2",
    "Double positive" = "green"
  )) +
  theme_classic() +
  ylab("Fraction of expressing cells") +
  xlab("Cluster") +
  scale_y_continuous(labels = scales::percent) +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x  = element_text(size = 20, angle = 45, hjust = 1, vjust = 1),
    axis.text.y  = element_text(size = 20),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 18)
    
  )
#####now for cells expressing each gene individually
human_OVOL2 <- data %>%
  filter(species == "Human", OVOL2 > 0 & ZNF90 == 0) %>%
  mutate(Gene = "ZNF90- & OVOL2+")

human_ZNF90 <- data %>%
  filter(species == "Human", ZNF90 > 0 & OVOL2 == 0) %>%
  mutate(Gene = "ZNF90+ & OVOL2-")

human_double_pos <- data %>%
  filter(species == "Human", OVOL2 > 0, ZNF90 > 0) %>%
  mutate(Gene = "ZNF90+ & OVOL2+")

human_none <- data %>%
  filter(species == "Human", OVOL2 == 0, ZNF90 == 0) %>%
  mutate(Gene = "ZNF90- & OVOL2-")

combined <- bind_rows(
  human_OVOL2,
  human_ZNF90,
  human_double_pos,
  human_none   # remove this line if you don't want it
)

# plot_data <- combined %>%
#   group_by(ident, Gene, Phase) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(ident, Gene) %>%
#   mutate(prop = n / sum(n))
# 
# cluster_order <- combined %>%
#   count(ident) %>%
#   arrange(desc(n)) %>%
#   pull(ident)
# 
# plot_data$ident <- factor(plot_data$ident, levels = cluster_order)


selected_clusters <- c("aRG1","aRG3","bRG1")

plot_data <- combined %>%
  filter(ident %in% selected_clusters) %>%
  group_by(ident, Gene, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(ident, Gene) %>%
  mutate(prop = n / sum(n))

cluster_order <- plot_data %>%
  count(ident) %>%
  arrange(desc(n)) %>%
  pull(ident) %>%
  unique()


plot_data$ident <- factor(plot_data$ident, levels = cluster_order)
gene_order <- c(
  "ZNF90+ & OVOL2-",
  "ZNF90- & OVOL2+",
  "ZNF90- & OVOL2-",
  "ZNF90+ & OVOL2+"
)

plot_data$Gene <- factor(plot_data$Gene, levels = gene_order)


ggplot(plot_data, aes(x = ident, y = prop, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Gene) +
  theme_classic() +
  ylab("Fraction of expressing cells") +
  xlab("Cluster") +
  scale_y_continuous(labels = scales::percent)+theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x  = element_text(size = 18),
    axis.text.y  = element_text(size = 18), 
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )

##########################calculate enrichment in the coexpression of both genes!

data <- FetchData(
  seurat_obj_combined,
  vars = c("OVOL2", "ZNF90", "species", "ident")
) %>%
  filter(species == "Human")

cluster_stats <- data %>%
  group_by(ident) %>%
  summarise(
    total_cells = n(),
    OVOL2_pos = sum(OVOL2 > 0),
    ZNF90_pos = sum(ZNF90 > 0),
    double_pos = sum(OVOL2 > 0 & ZNF90 > 0),
    .groups = "drop"
  ) %>%
  mutate(
    prop_OVOL2 = OVOL2_pos / total_cells,
    prop_ZNF90 = ZNF90_pos / total_cells,
    prop_double = double_pos / total_cells
  )



plot_data <- cluster_stats %>%
  select(ident, prop_OVOL2, prop_ZNF90, prop_double) %>%
  pivot_longer(
    cols = starts_with("prop"),
    names_to = "type",
    values_to = "prop"
  )

ggplot(plot_data, aes(x = ident, y = prop, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  ylab("Fraction of cells") +
  xlab("Cluster") +
  scale_y_continuous(labels = scales::percent)

###around 20% cells expressing ovol2 also express ZNF90 but only 10% cells expressing ZNF90 express OVOL2

cluster_stats <- data %>%
  group_by(ident) %>%
  summarise(
    total = n(),
    OVOL2_pos = sum(OVOL2 > 0),
    ZNF90_pos = sum(ZNF90 > 0),
    double_pos = sum(OVOL2 > 0 & ZNF90 > 0),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    fisher = list(fisher.test(matrix(c(
      double_pos,
      OVOL2_pos - double_pos,
      ZNF90_pos - double_pos,
      total - OVOL2_pos - ZNF90_pos + double_pos
    ), nrow = 2))),
    p_value = fisher$p.value,
    odds_ratio = fisher$estimate,
    conf_low = fisher$conf.int[1],
    conf_high = fisher$conf.int[2]
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    log2_or = log2(odds_ratio)
  )

cluster_stats <- cluster_stats %>%
  arrange(p_adj) %>%
  mutate(ident = factor(ident, levels = rev(ident)))

ggplot(cluster_stats, aes(x = log2_or, y = ident)) +
  geom_point(aes(color = p_adj < 0.05), size = 3) +
  geom_errorbarh(aes(xmin = log2(conf_low), xmax = log2(conf_high)), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    values = c("grey", "red"),
    labels = c("", "FDR<0.05"),
    guide = guide_legend(override.aes = list(size = 3))
  ) +
  labs(
    x = "log2(Odds Ratio)",
    y = "Cluster",
    color = NULL
  ) +
  theme_classic() + 
  theme(
    axis.title = element_text(size = 14),   # increase axis label size
    axis.text = element_text(size = 12),    # increase tick label size
    legend.text = element_text(size = 12)   # increase legend text size
  )

cluster_stats_clean <- cluster_stats %>%
  select(-fisher)

write.table(
  cluster_stats_clean,
  "cluster_stats_OVOL2_ZNF90_coexpression.tsv",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE, sep = "\t"
)

######

library(monocle3)
library(dplyr)

genes_use <- c("ZNF90", "OVOL2")

# Expression threshold
# > 0 means detected expression
expr_threshold <- 0

# Get pseudotime
pt <- pseudotime(cds)
pt <- pt[colnames(cds)]

# Directly extract only the needed gene expression vectors
ZNF90_expr <- as.numeric(exprs(cds)["ZNF90", colnames(cds)])
OVOL2_expr <- as.numeric(exprs(cds)["OVOL2", colnames(cds)])

# Keep only cells with valid pseudotime
valid_cells <- is.finite(pt)

# Define cell groups
ZNF90_only_cells <- valid_cells & ZNF90_expr > expr_threshold & OVOL2_expr <= expr_threshold

OVOL2_only_cells <- valid_cells & OVOL2_expr > expr_threshold & ZNF90_expr <= expr_threshold

both_cells <- valid_cells & ZNF90_expr > expr_threshold & OVOL2_expr > expr_threshold
df_ZNF90_only <- data.frame(
  cell_id = colnames(cds)[ZNF90_only_cells],
  pseudotime = as.numeric(pt[ZNF90_only_cells])
)

df_OVOL2_only <- data.frame(
  cell_id = colnames(cds)[OVOL2_only_cells],
  pseudotime = as.numeric(pt[OVOL2_only_cells])
)

df_both <- data.frame(
  cell_id = colnames(cds)[both_cells],
  pseudotime = as.numeric(pt[both_cells])
)

nrow(df_ZNF90_only)
nrow(df_OVOL2_only)
nrow(df_both)
df_compare <- bind_rows(
  df_ZNF90_only %>% mutate(group = "ZNF90+"),
  df_OVOL2_only %>% mutate(group = "OVOL2+"),
  df_both %>% mutate(group = "ZNF90+ & OVOL2+")
)

df_compare$group<-factor((df_compare$group),levels = c("ZNF90+", "OVOL2+","ZNF90+ & OVOL2+"))

ggplot(df_compare, aes(x = group, y = pseudotime, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  theme_classic() +
  labs(
    x = "",
    y = "Pseudotime"
  ) +
  theme(
    legend.position = "none"
  )+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
    axis.text.y = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.title.x = element_text(size = 25)
  )


set.seed(123)

shapiro_results <- df_compare %>%
  group_by(group) %>%
  summarise(
    n = n(),
    shapiro_p = {
      x <- pseudotime[is.finite(pseudotime)]
      
      if (length(x) < 3) {
        NA_real_
      } else {
        if (length(x) > 5000) {
          x <- sample(x, 5000)
        }
        shapiro.test(x)$p.value
      }
    },
    .groups = "drop"
  )

shapiro_results
kruskal.test(pseudotime ~ group, data = df_compare)
pairwise.wilcox.test(
  x = df_compare$pseudotime,
  g = df_compare$group,
  p.adjust.method = "BH"
)


