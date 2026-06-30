
#setwd("~/Downloads/ZNF90_output_HOMER_RepeatMasker")

##### now let's annotate the peaks to gene symbols and find potential CREs in the peaks_in_repeats
library(data.table)
library(dplyr)
library(biomaRt)

############################################################
# 1. LOAD HOMER annotatePeaks OUTPUT
############################################################

annot <- fread("ZNF90_homer_annotation.txt")

# Check column names
colnames(annot)

############################################################
# 2. EXTRACT ENSEMBL IDs
############################################################

# Your IDs seem to look like:
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
# 9. FILTER FOR TSS DISTANCE <= 100 kb
############################################################

repeat_CRE_candidates <- repeat_annot %>%
  filter(abs(`Distance to TSS`) <= 2000)
##a lot are AluS for ZNF90

############################################################
# 10. SAVE OUTPUT
############################################################

write.csv(
  repeat_CRE_candidates,
  "repeat_CRE_candidates_2kb_ZNF90.csv",
  row.names = FALSE
)
