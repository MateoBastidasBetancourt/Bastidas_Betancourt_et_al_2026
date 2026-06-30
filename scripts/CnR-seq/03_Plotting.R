library(ggplot2)
library(ggnewscale)
library(tidyverse)
library(tidyr)

##CnR plots: based on overall annotation stats


############################
# DATA
############################

df <- data.frame(
  Gene = c(
    "ZNF90","ZNF90","ZNF90","ZNF90","ZNF90",
    "OVOL2","OVOL2","OVOL2","OVOL2","OVOL2"
  ),
  
  Feature = c(
    "TTS","Exon","Intron","Intergenic","Promoter",
    "TTS","Exon","Intron","Intergenic","Promoter"
  ),
  
  Log2Ratio = c(
    0.337,0.777,0.006,-0.248,0.791,
    0.618,2.082,-0.060,-0.476,2.362
  ),
  
  LogP_FDR = c(
    11.22924,43.29409,0.97200,58.15510,64.25997,
    166.85306,7079.60500,60.17215,3771.26212,6689.08103
  )
)

############################
# CALCULATE OR
############################

df$OR <- 2^(df$Log2Ratio)

############################
# CALCULATE -log10(FDR)
############################
# HOMER LogP_FDR is already -log10(FDR)

df$neglog10FDR <- df$LogP_FDR

############################
# FEATURE ORDER
############################

df$Feature <- factor(
  df$Feature,
  levels = c(
    "Promoter",
    "Exon",
    "Intron",
    "TTS",
    "Intergenic"
  )
)
df$Gene <- factor(df$Gene, levels = c("ZNF90","OVOL2"))
df$LogP_plot <- pmin(df$LogP_FDR, 50)
############################
# PLOT
df$y <- NA

df$y[df$Feature == "TTS" & df$Gene == "OVOL2"] <- 0.7
df$y[df$Feature == "TTS" & df$Gene == "ZNF90"] <- 0.8

df$y[df$Feature == "Exon" & df$Gene == "OVOL2"] <- 1.2
df$y[df$Feature == "Exon" & df$Gene == "ZNF90"] <- 1.3

df$y[df$Feature == "Promoter" & df$Gene == "OVOL2"] <- 1.7
df$y[df$Feature == "Promoter" & df$Gene == "ZNF90"] <- 1.8

ggplot(df[df$Feature %in% c("TTS","Exon","Promoter"),],
       aes(x = OR,
           y = y,
           color = Gene)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey40"
  ) +
  scale_color_manual(
    values = c("OVOL2" = "#D55E00", "ZNF90" = "#0072B2")
  ) + guides(
    color = guide_legend(order = 2),
    size  = guide_legend(order = 1)
  ) +
  geom_segment(
    aes(x = 1,
        xend = OR,
        yend = y),
    linewidth = 0.8,
    alpha = 0.7
  ) +
  geom_point(
    aes(size = LogP_plot),
    alpha = 0.95
  ) +
  scale_y_continuous(
    breaks = c(0.75, 1.25, 1.75),
    labels = c("TTS", "Exon", "Promoter"),
    limits = NULL,
    expand = c(0.4, 0.)
  ) +
  labs(
  x = "Odds Ratio",
  y = ""
  )  + scale_size_continuous(
    name = expression(-log[10]("FDR")),
    range = c(3, 7)
  ) + 
  theme_classic(base_size = 14)+
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 19),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 21),
    legend.text = element_text(size = 21),
    legend.position = "right",
    axis.title.x = element_text(size = 21),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    
    axis.line = element_line(color = "black")
  )
  

############################

library(ggplot2)
library(ggnewscale)

############################
# TRUE HOMER ANNOTATION COUNTS
############################

ovol2_counts <- c(
  Promoter = 7762,
  Exon = 9917,
  Intron = 37037,
  TTS = 2009,
  Intergenic = 33123
)

znf90_counts <- c(
  Promoter = 462,
  Exon = 320,
  Intron = 3241,
  TTS = 339,
  Intergenic = 2194
)

############################
# PERCENTAGES
############################

ovol2_pct <- 100 * ovol2_counts / sum(ovol2_counts)
znf90_pct <- 100 * znf90_counts / sum(znf90_counts)

############################
# REPEAT CONTENT
############################

ovol2_repeat <- c(
  Repetitive = 5344,
  Non_repetitive = 91983 - 5344
)

znf90_repeat <- c(
  Repetitive = 708,
  Non_repetitive = 6599 - 708
)

ovol2_repeat_pct <- 100 * ovol2_repeat / sum(ovol2_repeat)
znf90_repeat_pct <- 100 * znf90_repeat / sum(znf90_repeat)

############################
# BUILD SEPARATE DATAFRAMES
############################

annotation_df <- data.frame(
  Gene = rep(c("OVOL2", "ZNF90"), each = 5),
  x = rep(c(0.85, 1.85), each = 5),
  Category = factor(
    rep(names(ovol2_pct), 2),
    levels = c("Promoter", "Exon", "Intron", "TTS", "Intergenic")
  ),
  Percent = c(as.numeric(ovol2_pct), as.numeric(znf90_pct))
)

repeat_df <- data.frame(
  Gene = rep(c("OVOL2", "ZNF90"), each = 2),
  x = rep(c(1.15, 2.15), each = 2),
  Category = factor(
    gsub("_", "-", rep(names(ovol2_repeat_pct), 2)),
    levels = c("Repetitive", "Non-repetitive")
  ),
  Percent = c(as.numeric(ovol2_repeat_pct), as.numeric(znf90_repeat_pct))
)

############################
# PLOT
############################

annotation_df$x <- ifelse(annotation_df$Gene == "OVOL2", 0.70, 1.70)
repeat_df$x     <- ifelse(repeat_df$Gene == "OVOL2", 1.30, 2.30)

annotation_df$x <- ifelse(annotation_df$Gene == "OVOL2", 0.60, 1.60)
repeat_df$x     <- ifelse(repeat_df$Gene == "OVOL2", 1.40, 2.40)
gene_centres <- c(
  OVOL2 = 1,
  ZNF90 = 2.8
)

bar_offset <- 0.25

annotation_df$x <- gene_centres[annotation_df$Gene] - bar_offset
repeat_df$x     <- gene_centres[repeat_df$Gene] + bar_offset

annotation_df$Gene <- factor(annotation_df$Gene,levels=c("ZNF90","OVOL2"))
repeat_df$Gene <- factor(repeat_df$Gene,levels=c("ZNF90","OVOL2"))

annotation_df$Gene <- factor(annotation_df$Gene, levels = c("ZNF90", "OVOL2"))
repeat_df$Gene <- factor(repeat_df$Gene, levels = c("ZNF90", "OVOL2"))

annotation_df$x <- as.numeric(annotation_df$Gene) - 0.25
repeat_df$x     <- as.numeric(repeat_df$Gene) + 0.25

gene_centres <- setNames(
  seq_along(levels(annotation_df$Gene)),
  levels(annotation_df$Gene)
)
bar_width <- 0.4
offset <- 0.4

ggplot() +
  geom_col(
    data = annotation_df,
    aes(x = x, y = Percent, fill = Category),
    width = 0.4,
    colour = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    name = "Genome annotations",
    values = c(
      Promoter = "#1b9e77",
      Exon = "#d95f02",
      Intron = "#7570b3",
      TTS = "#e7298a",
      Intergenic = "#66a61e"
    )
  ) +
  ggnewscale::new_scale_fill() +
  geom_col(
    data = repeat_df,
    aes(x = x, y = Percent, fill = Category),
    width = 0.4,
    colour = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    name = "Repetitive DNA content",
    values = c(
      Repetitive = "#e41a1c",
      `Non-repetitive` = "#377eb8"
    )
  ) +
  scale_x_continuous(
             breaks = gene_centres,
             labels = names(gene_centres),
             limits = c(
                   1 - offset - bar_width / 2,
                   n_genes + offset + bar_width / 2
               ),
             expand = c(0, 0)
    ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    x = "",
    y = "Percentage of peaks"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 25),
    axis.text.y = element_text(size=20),
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    legend.position = "right"
  )

#### now TE content:
ovol2_df <- data.frame(
  TF = "OVOL2",
  repeat_class = c("HERV", "STR", "ERV", "THE1", "AluY", "AluS"),
  OR = c(3.05889897, 1.86965802, 1.58450373, 1.52407327, 1.41755219, 1.10206261),
  FDR = c(1.099089e-13, 3.669678e-23, 5.805086e-14, 9.618219e-05, 3.761755e-11, 9.682997e-03)
)

znf90_df <- data.frame(
  TF = "ZNF90",
  repeat_class = c("HERV", "MST", "THE1", "AluS"),
  OR = c(3.3476267, 2.0317872, 1.9914589, 1.4172164),
  FDR = c(0.012535996, 0.034387260, 0.030051390, 0.002394863)
)
or_df <- bind_rows(ovol2_df, znf90_df) %>%
  mutate(
    neglog10FDR = -log10(FDR)
  )

repeat_order <- or_df %>%
  group_by(repeat_class) %>%
  summarise(max_OR = max(OR), .groups = "drop") %>%
  arrange(max_OR) %>%
  pull(repeat_class)

or_df <- or_df %>%
  mutate(
    repeat_class = factor(repeat_class, levels = repeat_order),
    y_base = as.numeric(repeat_class),
    y_pos = y_base + ifelse(TF == "OVOL2", -0.12, 0.12)
  )
or_df$TF <- factor(or_df$TF, levels = c("ZNF90","OVOL2"))

ggplot(or_df, aes(x = OR, y = y_pos, color = TF)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_segment(
    aes(x = 1, xend = OR, yend = y_pos),
    linewidth = 0.8,
    alpha = 0.7
  ) +
  geom_point(aes(size = neglog10FDR), alpha = 0.95) +
  scale_size_continuous(
    name = expression(-log[10]("FDR")),
    range = c(3, 7)
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(or_df$repeat_class)),
    labels = levels(or_df$repeat_class)
  ) +
  scale_x_log10(
    breaks = c(1, 1.25, 1.5, 2, 3, 4),
    labels = c("1", "1.25", "1.5", "2", "3", "4")
  ) +
  scale_color_manual(
    values = c("OVOL2" = "#D55E00", "ZNF90" = "#0072B2")
  ) +
  labs(
    x = "Odds Ratio",
    y = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25),
    legend.position = "right"
  )

###### now common targets:

library(VennDiagram)
library(grid)

znf90_genes
ovol2_genes

znf <- unique(na.omit(znf90_genes))
ovol <- unique(na.omit(ovol2_genes))

shared <- intersect(znf, ovol)

grid.newpage()

venn.plot <- draw.pairwise.venn(
  area1 = length(znf),
  area2 = length(ovol),
  cross.area = length(shared),
  category = c("", ""),
  fill = c("#0072B2", "#D55E00"),
  alpha = c(0.5, 0.5),
  lty = "blank",
  cex = 0,
  cat.cex = 0,
  scaled = TRUE
)

grid.draw(venn.plot)
# pct_znf_in_ovol
# # 47.8%
# 
# pct_ovol_in_znf
# # 7.8%

#################
library(ggplot2)
library(dplyr)
library(tidyr)

# Create dataframe
df <- data.frame(
  Feature = c("Promoter", "Exon", "Intron", "TTS", "Intergenic"),
  OVOL2 = c(8.64, 9.43, 54.72, 1.28, 25.93),
  ZNF90 = c(10.30, 5.67, 59.61, 5.59, 18.84)
)

# Set feature order
df$Feature <- factor(
  df$Feature,
  levels = c("Promoter", "Exon", "Intron", "TTS", "Intergenic")
)

# Convert to long format
df_long <- df %>%
  pivot_longer(
    cols = c(OVOL2, ZNF90),
    names_to = "Gene",
    values_to = "Percentage"
  )
df_long$Gene <- factor(df_long$Gene, levels = c("ZNF90", "OVOL2"))
# Plot stacked proportional bars
# Plot stacked proportional bars

ggplot(df_long, aes(x = Gene, y = Percentage, fill = Feature)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c(
    "Promoter"   = "#1B9E77",  # teal
    "Exon"       = "#D95F02",  # orange
    "Intron"     = "#7570B3",  # purple
    "TTS"        = "#E7298A",  # pink
    "Intergenic" = "#66A61E"   # green
  )) +
  labs(
    x = "",
    y = "Percentage of Peaks",
    fill = "Genome annotations"
  ) +
  theme_classic(base_size = 14) + geom_col(
    data = df_long,
    aes(x = Gene, y = Percentage, fill = Feature),
    width = 0.7,
    colour = "white",
    linewidth = 0.2
  ) + theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 22),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.position = "right"
  )

##########################
#plot genes with common peaks >1
df <- data.frame(
  Target = c(
    "DPP6", "ZNF423", "GABBR2", "SMOC2", "TPO", "CALN1",
    "DLGAP2", "HDAC4", "KCTD15", "SEPTIN9", "TCERG1L",
    "ZNF536", "ADCY1", "ADCY5", "APBA2", "CACNA1C",
    "CACNA1E", "CACNA2D2"
  ),
  OVOL2 = c(
    27, 17, 21, 20, 19, 19,
    23, 56, 18, 18, 24,
    27, 20, 16, 17, 15,
    17, 15
  ),
  ZNF90 = c(
    8, 5, 4, 4, 4, 3,
    3, 3, 3, 3, 3,
    3, 2, 2, 2, 2,
    2, 2
  )
)

# Order genes by total number of peaks
df$Total <- df$OVOL2 + df$ZNF90
df$Target <- factor(df$Target, levels = df$Target[order(df$Total, decreasing = TRUE)])

# Convert to long format, with ZNF90 first
df_long <- pivot_longer(
  df,
  cols = c("ZNF90", "OVOL2"),
  names_to = "Gene",
  values_to = "Peaks"
)

# Set TF order: ZNF90 first, OVOL2 second
df_long$Gene <- factor(df_long$Gene, levels = c("ZNF90", "OVOL2"))


ggplot(df_long, aes(x = Target, y = Peaks, fill = Gene)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "white"
  ) +
  geom_text(
    aes(label = Peaks),
    position = position_dodge(width = 0.8),
    vjust = -0.4,
    size = 7
  ) +
  scale_fill_manual(values = c(
    "ZNF90" = "#0072B2",
    "OVOL2" = "#D55E00"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    y = "Number of peaks",
    fill = "Target"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 25),
    axis.title.x = element_text(size = 0),
    legend.title = element_text(25),
    legend.text = element_text(size = 25)
  )

####plot enrichr data for the 100 promoters
df <- data.frame(
  term = c(
    "RUNX1", "WT1",
    "MYT1", "ARX", "CASZ1", "DPF1", "FOXG1",
    "NEUROD2", "PRDM16", "RAPGEF4", "RFX1",
    "KDM5B", "GATA3", "FOXA1"
  ),
  
  pvalue = c(
    0.000002794, 0.000004173,
    0.00001914, 0.0001286, 0.0001286, 0.0001286, 0.0001286,
    0.0001286, 0.0001286, 0.0001286, 0.00001914,
    0.00001578, 0.0002900, 0.001085
  ),
  
  padj = c(
    0.007201, 0.007201,
    0.01740, 0.01740, 0.01740, 0.01740, 0.01740,
    0.01740, 0.01740, 0.01740, 0.01740,
    0.001973, 0.01355, 0.02713
  ),
  
  odds_ratio = c(
    12.84, 5.31,
    6.69, 5.86, 5.86, 5.86, 5.86,
    5.86, 5.86, 5.86, 6.69,
    18.01, 26.73, 50.74
  ),
  
  combined_score = c(
    164.16, 65.73,
    72.65, 52.50, 52.50, 52.50, 52.50,
    52.50, 52.50, 52.50, 72.65,
    199.08, 217.71, 346.39
  )
)
df$logOR <- log2(df$odds_ratio)
df$neglogFDR <- -log10(df$padj)

df$term <- factor(df$term, levels = df$term)
ggplot(df, aes(x = logOR, y = term)) +
  geom_point(aes(size = neglogFDR, color = combined_score)) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(
    x = expression(log[2]~"Odds Ratio"),
    size = expression(-log[10]~"FDR"),
    color = "Combined score"
  ) +
  theme(
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    axis.title.y = element_text(size = 0),
    axis.title.x = element_text(size = 25),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 19)
  ) +
  scale_size_continuous(
    name = expression(-log[10]("FDR")),
    range = c(3, 10)
  )

###only OVOL2
df <- data.frame(
  term = c(
    "MYBL2", "KMT2A",
    "CHD1", "ZNF592", "RBM27", "ZNF644",
    "CUL5", "BRD9", "SP3"
  ),
  
  pvalue = c(
    2.44E-10, 2.53E-09,
    1.03E-09, 1.03E-09, 4.23E-09, 4.23E-09,
    8.42E-09, 1.65E-08, 1.65E-08
  ),
  
  padj = c(
    0.000001766, 0.000009145,
    8.83E-07, 8.83E-07, 0.00000182, 0.00000182,
    0.000002899, 0.000004067, 0.000004067
  ),
  
  odds_ratio = c(
    1.66, 1.33,
    2.03, 2.03, 1.98, 1.98,
    1.95, 1.93, 1.93
  ),
  
  combined_score = c(
    36.84, 26.26,
    42.08, 42.08, 38.14, 38.14,
    36.28, 34.49, 34.49
  )
)
df$logOR <- log2(df$odds_ratio)
df$neglogFDR <- -log10(df$padj)
df$term <- factor(df$term, levels = rev(df$term))
ggplot(df, aes(x = logOR, y = term)) +
  geom_point(aes(size = neglogFDR, color = combined_score)) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(
    x = expression(log[2]~"Odds Ratio"),
    size = expression(-log[10]~"FDR"),
    color = "Combined score"
  ) +
  theme(
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    axis.title.y = element_text(size = 0),
    axis.title.x = element_text(size = 25),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 19)
  ) +
  scale_size_continuous(
    name = expression(-log[10]("FDR")),
    range = c(3, 10)
  )

### ZNF90 has PRDM16 and RUNX1 as well!

