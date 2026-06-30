

#setwd("/mnt/vast-standard/home/cbastid/u15756/macs3_test/ZNF90_output_HOMER_RepeatMasker")
### now dataframe with all TEs:
df_TEs <- data.frame(
  class = c("Alu", "LTR_ERV", "L1", "L2", "MIR", "DNA_transposon", "Simple_repeat", "L3", "Other"),
  
  a = c(225, 159, 157, 93, 96, 80, 52, 8, 21),
  
  b = c(193468, 110972, 149201, 86080, 109258, 90417, 129908, 10656, 22110),
  
  c = c(666, 732, 734, 798, 795, 811, 839, 883, 870),
  
  d = c(708602, 791098, 752869, 815990, 792812, 811653, 772162, 891414, 879960)
)

df_TEs$p_value <- NA
df_TEs$odds_ratio <- NA

for (i in 1:nrow(df_TEs)) {
  mat <- matrix(c(df_TEs$a[i], df_TEs$b[i], df_TEs$c[i], df_TEs$d[i]), nrow = 2)
  test <- fisher.test(mat)
  
  df_TEs$p_value[i] <- test$p.value
  df_TEs$odds_ratio[i] <- test$estimate
}
df_TEs$log2_OR <- log2(df_TEs$odds_ratio)
df_TEs$FDR <- p.adjust(df_TEs$p_value, method = "BH")

###Alu and LTR_ERVs are enriched but simple_repeats are depleted. Let's analyze subcategories of Alu and ERV/LTRs
# -------------------------
# INPUT DATA
# -------------------------

dataset <- c(
  "LTR_ERV|MST"=16,
  "LTR_ERV|LTR_misc"=19,
  "LTR_ERV|THE1"=18,
  "LTR_ERV|MLT"=50,
  "LTR_ERV|LTR"=37,
  "LTR_ERV|ERV_other"=8,
  "LTR_ERV|HERV_other"=9,
  "LTR_ERV|HERVK"=2,
  
  "Alu|FLAM"=4,
  "Alu|AluJ"=46,
  "Alu|AluY"=27,
  "Alu|FRAM"=1,
  "Alu|AluS"=147
)

genome <- c(
  "LTR_ERV|MST"=8054,
  "LTR_ERV|LTR_misc"=15584,
  "LTR_ERV|THE1"=9253,
  "LTR_ERV|MLT"=44567,
  "LTR_ERV|LTR"=27040,
  "LTR_ERV|ERV_other"=3565,
  "LTR_ERV|HERV_other"=2744,
  "LTR_ERV|HERVK"=324,
  
  "Alu|FLAM"=6711,
  "Alu|AluJ"=51715,
  "Alu|AluY"=22927,
  "Alu|FRAM"=1159,
  "Alu|AluS"=110483
)

# totals
dataset_total <- 891   # from your previous sum
genome_total  <- sum(genome) + (902961 - sum(genome))  # OR just use 902961 if that's total TE count

genome_total <- 902961  # safer since you computed it already

# -------------------------
# BUILD TABLE
# -------------------------

results <- data.frame(
  subclass = names(dataset),
  a = dataset,
  c = genome[names(dataset)]
)

results$b <- dataset_total - results$a
results$d <- genome_total - results$c

# -------------------------
# OR + Fisher
# -------------------------

results$OR <- (results$a * results$d) / (results$b * results$c)
results$log2OR <- log2(results$OR)

results$pvalue <- mapply(function(a,b,c,d) {
  fisher.test(matrix(c(a,b,c,d), nrow=2))$p.value
}, results$a, results$b, results$c, results$d)

# -------------------------
# FDR correction
# -------------------------

results$FDR <- p.adjust(results$pvalue, method="BH")

# -------------------------
# SORT (most enriched first)
# -------------------------

results <- results[order(-results$log2OR), ]

print(results) #enrichment in AluS, THE1A/B and MSTA/B                    
                          # OR  FDR
# LTR_ERV|HERVK      6.2675291 0.090301007
# LTR_ERV|HERV_other 3.3476267 0.012535996
# LTR_ERV|ERV_other  2.2857078 0.070891343
# LTR_ERV|MST        2.0317872 0.034387260
# LTR_ERV|THE1       1.9914589 0.030051390
# Alu|AluS           1.4172164 0.002394863
# LTR_ERV|LTR        1.4034667 0.090301007
# LTR_ERV|LTR_misc   1.2406987 0.474450534
# Alu|AluY           1.1995055 0.474450534
# LTR_ERV|MLT        1.1451102 0.474450534
# Alu|AluJ           0.8960653 0.559231390
# Alu|FRAM           0.8742542 1.000000000
# Alu|FLAM           0.6022521 0.510000241

###now same abundance test for OVOL2:
# Genome abundance
genome <- c(
  # LTR / ERV
  MST = 5199,
  LTR_misc = 10186,
  THE1 = 5957,
  MLT = 28600,
  ERV_other = 2368,
  LTR = 18924,
  HERV_other = 1823,
  HERVK = 346,
  
  # Alu
  FLAM = 5524,
  AluJ = 33574,
  AluY = 28363,
  FRAM = 2465,
  Alu_other = 743,
  AluS = 81368,
  
  # Simple repeats
  di_repeat = 23659,
  mono_repeat = 11959,
  penta_repeat = 7597,
  tri_repeat = 8397,
  tetra_repeat = 18656,
  long_repeat = 8316,
  hexa_repeat = 16112,
  
  # DNA transposons
  PABL = 158,
  Zaphod = 645,
  Tigger = 9828,
  MER = 42152,
  Charlie = 8279
)

# =========================
# PEAK OVERLAPS
# =========================

peaks <- c(
  # LTR / ERV
  LTR_misc = 46,
  MST = 70,
  THE1 = 106,
  MLT = 308,
  ERV_other = 41,
  LTR = 343,
  HERV_other = 65,
  HERVK = 4,
  
  # Alu
  FLAM = 29,
  AluJ = 305,
  AluY = 459,
  FRAM = 14,
  Alu_other = 5,
  AluS = 1033,
  
  # Simple repeats
  di_repeat = 149,
  mono_repeat = 6,
  penta_repeat = 135,
  tri_repeat = 263,
  tetra_repeat = 156,
  long_repeat = 193,
  hexa_repeat = 342,
  
  # DNA transposons
  PABL = 2,
  Zaphod = 4,
  Tigger = 36,
  MER = 351,
  Charlie = 22
)

total_genome <- sum(genome)
total_peaks <- sum(peaks)

results <- data.frame()

pvals <- c()

for(te in names(genome)) {
  
  mat <- matrix(c(
    peaks[te],
    total_peaks - peaks[te],
    genome[te],
    total_genome - genome[te]
  ),
  nrow = 2,
  byrow = TRUE)
  
  ft <- fisher.test(mat)
  
  expected <- total_peaks * (genome[te] / total_genome)
  fold <- peaks[te] / expected
  
  pvals <- c(pvals, ft$p.value)
  
  results <- rbind(results,
                   data.frame(
                     TE = te,
                     observed = peaks[te],
                     expected = round(expected, 2),
                     fold_enrichment = round(fold, 3),
                     odds_ratio = as.numeric(ft$estimate)
                   )
  )
}

results$FDR <- p.adjust(pvals, method = "BH")

results <- results[order(-results$fold_enrichment), ]

results[,c(4,5,6)]

# odds_ratio          FDR
# HERV_other   3.05889897 1.099089e-13
# tri_repeat   2.76426463 3.838049e-42
# long_repeat  2.01535293 4.823662e-17
# hexa_repeat  1.86965802 3.669678e-23
# LTR          1.58450373 5.805086e-14
# THE1         1.52407327 9.618219e-05
# penta_repeat 1.52546298 1.012760e-05
# ERV_other    1.47525224 2.208688e-02
# AluY         1.41755219 3.761755e-11
# MST          1.14613783 3.061120e-01
# AluS         1.10206261 9.682997e-03
# PABL         1.07541145 7.388175e-01
# HERVK        0.98213551 1.000000e+00
# MLT          0.90863693 1.364927e-01
# AluJ         0.75516120 2.181953e-06
# tetra_repeat 0.69995689 9.762359e-06
# MER          0.68260204 1.612003e-12
# Alu_other    0.57123264 3.260144e-01
# di_repeat    0.51904476 2.215015e-17
# Zaphod       0.52643786 3.061120e-01
# FRAM         0.48086901 4.889466e-03
# FLAM         0.44240094 1.545764e-06
# LTR_misc     0.37727850 5.166648e-14
# Tigger       0.30562395 2.215015e-17
# Charlie      0.22192508 2.669619e-19
# mono_repeat  0.04134514 9.485539e-51
###OVOL2 binds HERVs, LTR, trinucleotide repeats and AluY (young)
##### now let's annotate the peaks to gene symbols and find potential CREs in the peaks_in_repeats

#setwd("~/Downloads/ZNF90_output_HOMER_RepeatMasker")

##### now let's annotate the peaks to gene symbols and find potential CREs in the peaks_in_repeats
library(data.table)
library(dplyr)
library(biomaRt)

############################################################
# 1. LOAD HOMER annotatePeaks OUTPUT
############################################################

annot <- fread("ZNF90_homer_annotation.txt")
annot_ZNF90 <- annot
# Check column names
colnames(annot)

############################################################
# 2. EXTRACT ENSEMBL IDs
############################################################

# Your IDs seem to look like:
# ENST00000498319.2

annot$ensembl_transcript <- annot$`Nearest Unigene`

# remove transcript version number
annot$ensembl_transcript <- sub("\\..*", "", annot$ensembl_transcript)

############################################################
# 3. CONNECT TO ENSEMBL
############################################################

mart <- useMart(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl"
)

############################################################
# 4. RETRIEVE GENE SYMBOLS
############################################################

gene_map <- getBM(
  attributes = c(
    "ensembl_transcript_id",
    "external_gene_name",
    "description",
    "gene_biotype"
  ),
  filters = "ensembl_transcript_id",
  values = unique(annot$ensembl_transcript),
  mart = mart
)

############################################################
# 5. MERGE BACK INTO ANNOTATION TABLE
############################################################

annot2 <- annot %>%
  left_join(
    gene_map,
    by = c("ensembl_transcript" = "ensembl_transcript_id")
  )

write.table(annot2,
            "annotated_peaks_symbols.txt",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

promoter_peaks <- annot2 %>%
  filter(grepl(c("promoter"), Annotation, ignore.case = TRUE))

############################################################
# 6. LOAD REPEAT-OVERLAPPING PEAKS
############################################################

# Example structure:
# chr1 1168000 1168178 chr1 1167999 1168144 G-rich 39 +

repeats <- fread(
  "peaks_in_repeats_ZNF90.bed",
  header = FALSE
)

colnames(repeats) <- c(
  "chr_peak",
  "start_peak",
  "end_peak",
  "chr_repeat",
  "start_repeat",
  "end_repeat",
  "repeat_name",
  "repeat_score",
  "strand"
)

############################################################
# 7. MATCH REPEAT PEAKS TO HOMER-ANNOTATED PEAKS
############################################################

# create coordinate keys
annot2$peak_key <- paste(
  annot2$Chr,
  annot2$Start,
  annot2$End,
  sep = "_"
)

repeats$peak_key <- paste(
  repeats$chr_peak,
  repeats$start_peak + 1,
  repeats$end_peak,
  sep = "_"
)

repeat_annot <- repeats %>%
  left_join(annot2, by = "peak_key")

############################################################
# 9. FILTER FOR TSS DISTANCE <= 2 kb
############################################################

repeat_CRE_candidates <- repeat_annot %>%
  filter(abs(`Distance to TSS`) <= 2000)
##a lot are AluS for ZNF90

############################################################
# 10. SAVE OUTPUT
############################################################

write.csv(
  repeat_CRE_candidates,
  "repeat_CRE_candidates_2kb.csv",
  row.names = FALSE
)

#########promoter - TE analysis for OVOL2

library(dplyr, lib.loc = "/usr/lib64/R/library")
library(data.table, lib.loc = "/usr/lib64/R/library")

#setwd("/mnt/vast-standard/home/cbastid/u15756/macs3_test/OVOL2_output_HOMER_RepeatMasker")
annot <- fread("OVOL2_homer_annotation.txt")
annot_OVOL2 <- annot

repeats <- fread(
  "peaks_in_repeats_OVOL2.bed",
  header = FALSE
)

colnames(repeats) <- c(
  "chr_peak",
  "start_peak",
  "end_peak",
  "chr_repeat",
  "start_repeat",
  "end_repeat",
  "repeat_name",
  "repeat_score",
  "strand"
)


# create coordinate keys
annot$peak_key <- paste(
  annot$Chr,
  annot$Start,
  annot$End,
  sep = "_"
)

repeats$peak_key <- paste(
  repeats$chr_peak,
  repeats$start_peak + 1,
  repeats$end_peak,
  sep = "_"
)

repeat_annot <- repeats %>%
  left_join(annot, by = "peak_key")



repeat_CRE_candidates_OVOL2 <- repeat_annot %>%
  filter(abs(`Distance to TSS`) <= 2000)


write.csv(
  repeat_CRE_candidates_OVOL2,
  "repeat_CRE_candidates_2kb_OVOL2.csv",
  row.names = FALSE
)



########binding overlap for both genes' annotated peaks
head(annot_OVOL2)
head(annot_ZNF90)
znf90_genes <- unique(
  annot_ZNF90$external_gene_name[
    !is.na(annot_ZNF90$external_gene_name) &
      annot_ZNF90$external_gene_name != ""
  ]
)

ovol2_genes <- unique(
  annot_OVOL2$`Gene Name`[
    !is.na(annot_OVOL2$`Gene Name`) &
      annot_OVOL2$`Gene Name` != ""
  ]
)

# -----------------------------
# BASIC COUNTS
# -----------------------------

length(znf90_genes) #4211
length(ovol2_genes) #25747

# -----------------------------
# SHARED GENES
# -----------------------------

shared_genes <- intersect(znf90_genes, ovol2_genes)

length(shared_genes)#2011

head(shared_genes)

head(annot_OVOL2$`Gene Name`)
head(annot_ZNF90$external_gene_name)

# Filter annot_OVOL2
annot_OVOL2_shared <- annot_OVOL2[
  annot_OVOL2$`Gene Name` %in% shared_genes,
]

# Filter annot_ZNF90
annot_ZNF90_shared <- annot_ZNF90[
  annot_ZNF90$external_gene_name %in% shared_genes,
]

#####
gene_counts_ZNF90 <- sort(
  table(annot_ZNF90_shared$external_gene_name),
  decreasing = TRUE
)

# View top genes
head(gene_counts_ZNF90, 100)
## Count occurrences of each gene
gene_counts_OVOL2 <- sort(
  table(annot_OVOL2_shared$`Gene Name`),
  decreasing = TRUE
)

# View top genes
head(gene_counts_OVOL2, 100)
#####overlap between genes present more than once in BOTH datasets
# | Gene     | OVOL2 peaks | ZNF90 peaks |
#   | -------- | ----------: | ----------: |
#   | DPP6     |          27 |           8 |
#   | ZNF423   |          17 |           5 |
#   | GABBR2   |          21 |           4 |
#   | SMOC2    |          20 |           4 |
#   | TPO      |          19 |           4 |
#   | CALN1    |          19 |           3 |
#   | DLGAP2   |          23 |           3 |
#   | HDAC4    |          56 |           3 |
#   | KCTD15   |          18 |           3 |
#   | SEPTIN9  |          18 |           3 |
#   | TCERG1L  |          24 |           3 |
#   | ZNF536   |          27 |           3 |
#   | ADCY1    |          20 |           2 |
#   | ADCY5    |          16 |           2 |
#   | APBA2    |          17 |           2 |
#   | CACNA1C  |          15 |           2 |
#   | CACNA1E  |          17 |           2 |
#   | CACNA2D2 |          15 |           2 |
#   


######overall category counts for overlap peaks between both genes:
annot_ZNF90_shared$SimpleAnnotation <- ifelse(
  grepl("^intron", annot_ZNF90_shared$Annotation), "Intron",
  
  ifelse(
    grepl("^exon", annot_ZNF90_shared$Annotation), "Exon",
    
    ifelse(
      grepl("^promoter-TSS", annot_ZNF90_shared$Annotation), "Promoter",
      
      ifelse(
        grepl("^TTS", annot_ZNF90_shared$Annotation), "TTS",
        
        ifelse(
          grepl("^Intergenic", annot_ZNF90_shared$Annotation), "Intergenic",
          "Other"
        )
      )
    )
  )
)

# Count categories
category_counts_ZNF90 <- sort(
  table(annot_ZNF90_shared$SimpleAnnotation),
  decreasing = TRUE
)

category_counts_ZNF90
# category_counts_ZNF90
# 
# Intron               Intergenic           Promoter            Exon               TTS 
# 1430(0,5960)        452(0,1884)        247 (0,102)        136(0,0566)        134(0,0558)

annot_OVOL2_shared$SimpleAnnotation <- ifelse(
  grepl("^intron", annot_OVOL2_shared$Annotation), "Intron",
  
  ifelse(
    grepl("^exon", annot_OVOL2_shared$Annotation), "Exon",
    
    ifelse(
      grepl("^promoter-TSS", annot_OVOL2_shared$Annotation), "Promoter",
      
      ifelse(
        grepl("^TTS", annot_OVOL2_shared$Annotation), "TTS",
        
        ifelse(
          grepl("^Intergenic", annot_OVOL2_shared$Annotation), "Intergenic",
          "Other"
        )
      )
    )
  )
)

# Count categories
category_counts_OVOL2 <- sort(
  table(annot_OVOL2_shared$SimpleAnnotation),
  decreasing = TRUE
)

category_counts_OVOL2
# Intron               Intergenic       Exon                 Promoter        TTS 
# 5366(0,547)       2543(0,259)        925(0,094)        847(0,0943)        126(0,0128) 

########binding overlap for both genes' in promoters
head(annot_OVOL2)
head(annot_ZNF90)

znf90_genes <- annot_ZNF90 %>%
  filter(grepl(c("promoter"), Annotation, ignore.case = TRUE)) ##462 promoter peaks

ovol2_genes <- annot_OVOL2 %>%
  filter(grepl(c("promoter"), Annotation, ignore.case = TRUE)) ##7767 promoter peaks

znf90_genes <- unique(
  znf90_genes$external_gene_name[
    !is.na(znf90_genes$external_gene_name) &
      znf90_genes$external_gene_name != ""
  ]
)

ovol2_genes <- unique(
  ovol2_genes$`Gene Name`[
    !is.na(ovol2_genes$`Gene Name`) &
      ovol2_genes$`Gene Name` != ""
  ]
)

length(znf90_genes)#393
length(ovol2_genes)#7036

shared_genes <- intersect(znf90_genes, ovol2_genes)

length(shared_genes) #=100
shared_genes
# [1] "ATPAF2"   "DUOX1"    "POLR3C"   "ZNF513"   "ARHGAP23" "PEX10"    "FAM53A"  
# [8] "FBN3"     "DMXL1"    "MTSS1"    "ATP6V1A"  "HLF"      "CACNA2D2" "AR"      
# [15] "SLC1A6"   "CDK2AP1"  "STK11"    "IL13RA1"  "NPL"      "PLXNB1"   "KCTD16"  
# [22] "ACSF3"    "MFSD5"    "HSF2"     "TP53AIP1" "EMC1"     "SF3B2"    "PLEKHG4" 
# [29] "STT3B"    "CNGB3"    "ZNF775"   "CCNE2"    "MORC2"    "ZBTB25"   "VIPR1"   
# [36] "ADARB1"   "TGFA"     "RPAP1"    "RRM1"     "DLX1"     "TBCD"     "FER"     
# [43] "TMEM38A"  "CAMK2B"   "ROR2"     "HDAC4"    "RBM15B"   "EPN1"     "CUL7"    
# [50] "TMEM92"   "MOGAT1"   "KCNJ3"    "NCAPD3"   "DYNC1LI2" "PRKAR2B"  "BANP"    
# [57] "GJA3"     "DYRK1A"   "POLR1D"   "PIK3R5"   "CACNA1A"  "OSM"      "ESRP2"   
# [64] "WDR82"    "GFPT1"    "LDHA"     "FDPS"     "APTX"     "RIOX2"    "GCN1"    
# [71] "TXNRD2"   "ZNF598"   "HSBP1"    "OGA"      "SNAP25"   "CSRP1"    "PTBP1"   
# [78] "DEPDC1"   "MEIS2"    "C1QTNF5"  "MPP2"     "OGT"      "ISYNA1"   "CREBZF"  
# [85] "HDAC5"    "HYDIN"    "GDE1"     "SPART"    "PRRT3"    "NCOA7"    "FANCA"   
# [92] "KIAA1755" "GAD1"     "MSI2"     "GSE1"     "ATR"      "NFKBID"   "SEPTIN9" 
# [99] "ZMYM3"    "SLC38A5"

###now for promoters with TEs
znf90_TE_promoters <- unique(
  repeat_CRE_candidates_2kb_ZNF90$external_gene_name[
    !is.na(repeat_CRE_candidates_2kb_ZNF90$external_gene_name) &
      repeat_CRE_candidates_2kb_ZNF90$external_gene_name != ""
  ]
)

ovol2_TE_promoters <- unique(
  repeat_CRE_candidates_OVOL2$`Gene Name`[
    !is.na(repeat_CRE_candidates_OVOL2$`Gene Name`) &
      repeat_CRE_candidates_OVOL2$`Gene Name` != ""
  ]
)

length(znf90_TE_promoters)#82
length(ovol2_TE_promoters)#621

shared_promoters_TE <- intersect(znf90_TE_promoters, ovol2_TE_promoters)
# shared_promoters_TE
# "MAN1C1"  "DEPDC1"  "SRSF11"  "S100A11" "LGR6"    "DNAH14" 

