library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(org.Mmu.eg.db)
library(dplyr)
library(GenomicRanges)
library(tidyr)
library(AnnotationDbi)

#data mining reilly et al 2015-ChIP seq h3k27ac and h3k4me2 from human and rhesus macaque
# -------------------------------
# 1. Set paths
# -------------------------------
bed_folder <- "~/Downloads/GSE63649_RAW_ChIP_seq_HS_MM/Human/" #or human
out_folder <- "~/Downloads/GSE63649_RAW_ChIP_seq_HS_MM/Human/CUTRUN_gene_scores_by_mark/"
dir.create(out_folder, showWarnings = FALSE)

# -------------------------------
# 2. Define histone marks / TF marks
# -------------------------------
marks <- c("H3K4me2", "H3K27ac")

# Create one folder per mark
for(mark in marks){
  dir.create(file.path(out_folder, mark), showWarnings = FALSE)
}

# -------------------------------
# 3. Load TxDb object
# -------------------------------
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#or human TxDb.Hsapiens.UCSC.hg19.knownGene

# Create promoter regions once (±3kb)
promoters_gr <- promoters(genes(txdb), upstream=3000, downstream=3000)
# -------------------------------
# 4. Function to detect mark from filename
# -------------------------------
get_mark <- function(filename){
  filename <- tolower(filename)
  if(grepl("h3k4me2", filename)) return("H3K4me2")
  if(grepl("h3k27ac", filename)) return("H3K27ac")
  return("Other")
}

# -------------------------------
# 5. Annotate and merge peaks per mark
# -------------------------------
bed_files <- list.files(bed_folder, pattern = "\\.bed$", full.names = TRUE)
goi_bed_df <- data.frame()
ovol2_coordinates_df <- data.frame()
genes_of_interest <- c("OVOL2", "SOX2", "GAPDH", "HBB")
######
for(mark in marks){
  
  mark_files <- bed_files[sapply(bed_files, function(f) get_mark(basename(f)) == mark)]
  if(length(mark_files) == 0) next
  
  combined_df <- NULL  # will store gene‑level peak counts across samples
  
  for(bed_file in mark_files){
    
    message("Processing: ", basename(bed_file))
    sample_name <- gsub(".bed.gz|.bed", "", basename(bed_file))
    
    bed_data <- read.table(bed_file, header = FALSE)
    peaks_gr <- GRanges(
      seqnames = bed_data$V1,
      ranges   = IRanges(bed_data$V2 + 1, bed_data$V3)
    )
    
    hits <- findOverlaps(peaks_gr, promoters_gr)
    if(length(hits) == 0) next
    
    ## ---------------------------
    ## Overlap‑level gene mapping
    ## ---------------------------
    gene_ids <- names(promoters_gr)[subjectHits(hits)]
    
    gene_symbols <- mapIds(
      org.Hs.eg.db,
      keys = gene_ids,
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    )
    
    ## ===========================
    ## 1️⃣ Gene‑level peak counts
    ## ===========================
    gene_counts <- table(gene_symbols)
    
    df_score <- data.frame(
      gene_symbol = names(gene_counts),
      score = as.integer(gene_counts),
      sample = sample_name,
      stringsAsFactors = FALSE
    )
    
    if(is.null(combined_df)){
      combined_df <- df_score
    } else {
      combined_df <- bind_rows(combined_df, df_score)
    }
    
    ## ===========================
    ## 2️⃣ GOI peak coordinates
    ## ===========================
    goi_hits <- which(gene_symbols %in% genes_of_interest)
    if(length(goi_hits) > 0){
      
      selected_peaks <- peaks_gr[unique(queryHits(hits)[goi_hits])]
      
      bed_df <- data.frame(
        chr    = as.character(seqnames(selected_peaks)),
        start  = start(selected_peaks) - 1,
        end    = end(selected_peaks),
        sample = sample_name,
        stringsAsFactors = FALSE
      )
      
      goi_bed_df <- bind_rows(goi_bed_df, bed_df)
    }
    
    rm(bed_data, peaks_gr, hits, gene_ids, gene_symbols)
    gc()
  }
  
  ## ---------------------------
  ## Save gene-level matrix
  ## ---------------------------
  gene_matrix <- combined_df %>%
    pivot_wider(
      names_from = sample,
      values_from = score,
      values_fill = NA
    ) %>%
    as.data.frame()
  
  rownames(gene_matrix) <- gene_matrix$gene_symbol
  gene_matrix$gene_symbol <- NULL
  
  write.table(
    gene_matrix,
    file = file.path(out_folder,
                     paste0(mark, "_gene_promoter_peakcount_matrix.txt")),
    sep = "\t",
    quote = FALSE
  )
  
  ## ---------------------------
  ## Save GOI peak coordinates
  ## ---------------------------
  write.table(
    goi_bed_df,
    file = file.path(out_folder, mark,
                     paste0(mark, "_GOI_promoter_peaks.bed")),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  out_bed <- file.path(
    out_folder,
    mark,
    paste0(mark, "_GOI_promoter_peaks.bed")
  )
  
  dir.create(file.path(out_folder, mark), showWarnings = FALSE)
  
  write.table(
    goi_bed_df,
    out_bed,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
}

goi_bed_df_human <- goi_bed_df

###################
# 1. Set paths
# -------------------------------
bed_folder <- "~/Downloads/GSE63649_RAW_ChIP_seq_HS_MM/Rhesus/" #or human
out_folder <- "~/Downloads/GSE63649_RAW_ChIP_seq_HS_MM/Rhesus/CUTRUN_gene_scores_by_mark/"
dir.create(out_folder, showWarnings = FALSE)

# -------------------------------
# 2. Define histone marks / TF marks
# -------------------------------
marks <- c("H3K4me2", "H3K27ac")

# Create one folder per mark
for(mark in marks){
  dir.create(file.path(out_folder, mark), showWarnings = FALSE)
}

# -------------------------------
# 3. Load genomic ranges object
# -------------------------------
ncbi_dataset_ <- read.delim("~/Downloads/ncbi_dataset.tsv") #https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_000002255.3/ from RheMac2: https://academic.oup.com/nar/article/33/suppl_1/D501/2505241?login=false

#colnames(df) <- c("start","end","chr","strand","description","gene_name")

# add chr prefix
ncbi_dataset_$Chromosome <- paste0("chr", ncbi_dataset_$Chromosome)

library(GenomicRanges)

genes_gr <- GRanges(
  seqnames = ncbi_dataset_$Chromosome,
  ranges = IRanges(start=ncbi_dataset_$Begin, end=ncbi_dataset_$End),
  strand = ifelse(ncbi_dataset_$Orientation=="plus","+","-"),
  gene_name = ncbi_dataset_$Symbol
)

promoters_gr <- promoters(
  genes_gr,
  upstream = 3000,
  downstream = 3000
)

# 4. Function to detect mark from filename
# -------------------------------
get_mark <- function(filename){
  filename <- tolower(filename)
  if(grepl("h3k4me2", filename)) return("H3K4me2")
  if(grepl("h3k27ac", filename)) return("H3K27ac")
  return("Other")
}

# -------------------------------
# 5. Annotate and merge peaks per mark
# -------------------------------
bed_files <- list.files(bed_folder, pattern = "\\.bed$", full.names = TRUE)
# Pre-map OVOL2 ENTREZ ID once
# Initialize OVOL2 coordinates dataframe
ovol2_coordinates_df <- data.frame()
goi_bed_df <- data.frame()
genes_of_interest <- c("OVOL2", "SOX2", "GAPDH", "HBB")

for(mark in marks){
  
  mark_files <- bed_files[sapply(bed_files, function(f) get_mark(basename(f)) == mark)]
  if(length(mark_files) == 0) next
  
  combined_df <- NULL
  
  for(bed_file in mark_files){
    
    message("Processing: ", basename(bed_file))
    
    # Read peaks
    bed_data <- read.table(bed_file, header = FALSE)
    
    # Convert to GRanges
    peaks_gr <- GRanges(
      seqnames = bed_data$V1,
      ranges   = IRanges(start = bed_data$V2 + 1, end = bed_data$V3)
    )
    
    # -----------------------------
    # Normalize chromosome names
    # -----------------------------
    seqlevels(peaks_gr)    <- sub("^chr0*", "chr", seqlevels(peaks_gr))
    seqlevels(promoters_gr) <- sub("^chr0*", "chr", seqlevels(promoters_gr))
    
    # Find overlaps with promoters
    hits <- findOverlaps(peaks_gr, promoters_gr)
    hit_promoters <- promoters_gr[subjectHits(hits)]
    hit_gene_names <- hit_promoters$gene_name
    
    # Count peaks per gene
    gene_counts <- table(hit_gene_names)
    df_score <- data.frame(
      gene_name = names(gene_counts),
      score     = as.numeric(gene_counts)
    )
    
    sample_name <- gsub(".bed.gz|.bed", "", basename(bed_file))
    
    # Ensure row exists
    if(!(sample_name %in% rownames(ovol2_coordinates_df))){
      ovol2_coordinates_df[sample_name, genes_of_interest] <- NA
    }
    
    # -------------------------------
    # Extract peaks overlapping genes of interest (BED-style)
    # -------------------------------
    
    subject_gene_names <- hit_promoters$gene_name
    
    # Find overlaps where promoter gene is in genes of interest
    goi_hit_indices <- which(subject_gene_names %in% genes_of_interest)
    
    if(length(goi_hit_indices) > 0){
      # Get the peak indices (from peaks_gr) corresponding to those overlaps
      peak_indices <- unique(queryHits(hits)[goi_hit_indices])
      selected_peaks <- peaks_gr[peak_indices]
      
      sample_name <- gsub(".bed.gz|.bed", "", basename(bed_file))
      
      bed_df <- data.frame(
        chr   = as.character(seqnames(selected_peaks)),
        start = start(selected_peaks) - 1,  # BED = 0-based
        end   = end(selected_peaks),
        sample = sample_name
      )
      
      # Append
      goi_bed_df <- rbind(goi_bed_df, bed_df)
    }
    
    if(is.null(combined_df)){
      combined_df <- df_score
    } else {
      combined_df <- full_join(combined_df, df_score, by="gene_name")
    }
    
    rm(bed_data, peaks_gr, hits, gene_counts, df_score)
    gc()
  }
  
  out_file <- file.path(out_folder, paste0(mark, "_combined_promoter_peakcounts.txt"))
  write.table(combined_df, out_file,
              row.names=FALSE, col.names=TRUE,
              quote=FALSE, sep="\t")
  out_bed <- file.path(out_folder, paste0(mark, "_GOI_peaks.bed"))
  
  write.table(goi_bed_df,
              out_bed,
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE,
              sep = "\t")
  rm(combined_df)
  gc()
}
goi_bed_rhesus <- goi_bed_df
