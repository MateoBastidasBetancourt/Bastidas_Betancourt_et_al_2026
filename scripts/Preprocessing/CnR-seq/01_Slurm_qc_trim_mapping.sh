#!/bin/bash
#SBATCH --job-name=CnR
#SBATCH --mail-user=c.bastidasbetancou@stud.uni-goettingen.de
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=300gb
#SBATCH --time=48:00:00
#SBATCH --partition=medium
#SBATCH --output=CnR.out
#SBATCH --error=CnR.err

set -euo pipefail
#
set -x
module load gcc/14.2.0
module load bowtie2/2.5.4
module load samtools/1.21
module load apptainer
BASE_DIR="$(pwd)/p1827"

TRIMMED_DIR="$(pwd)/trimmed"
FASTQC_DIR="$(pwd)/fastqc"



mkdir -p "$TRIMMED_DIR" "$FASTQC_DIR"

# Containers
CUTADAPT_SIF="/sw/container/bioinformatics/cutadapt-5.0.sif"
FASTQC_SIF="/sw/container/bioinformatics/fastqc-0.12.1.sif"
MULTIQC_SIF="/sw/container/bioinformatics/multiqc-1.27.1.sif"

# Loop over all sample directories
for sample_dir in "$BASE_DIR"/*/; do
    sample_group="$(basename "$sample_dir")"

    # Find R1 and R2 files
    r1=$(find "$sample_dir" -maxdepth 1 -name "*_R1_001.fastq.gz" | head -n 1)
    r2=$(find "$sample_dir" -maxdepth 1 -name "*_R2_001.fastq.gz" | head -n 1)

    if [[ -z "${r1:-}" || -z "${r2:-}" ]]; then
        echo "Skipping $sample_group: could not find both R1 and R2"
        continue
    fi

    mkdir -p "${TRIMMED_DIR}/${sample_group}" "${FASTQC_DIR}/${sample_group}"

    sample=$(basename "$r1")
    sample=${sample%%_R1_*}
    echo "Processing $sample_group / $sample"
    apptainer run \
    --bind "$(pwd):$(pwd)" \
    "$CUTADAPT_SIF" \
    cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -q 20 \
    -m 20 \
    --pair-filter=any \
    -o "${TRIMMED_DIR}/${sample_group}/${sample}_R1_trimmed.fastq.gz" \
    -p "${TRIMMED_DIR}/${sample_group}/${sample}_R2_trimmed.fastq.gz" \
    "$r1" "$r2" \
    > "${TRIMMED_DIR}/${sample_group}/${sample}_cutadapt.log"

    # Run FastQC
    apptainer run \
        --bind "$(pwd):$(pwd)" \
        "$FASTQC_SIF" \
        fastqc \
        --outdir "${FASTQC_DIR}/${sample_group}" \
        "${TRIMMED_DIR}/${sample_group}/${sample}_R1_trimmed.fastq.gz" \
        "${TRIMMED_DIR}/${sample_group}/${sample}_R2_trimmed.fastq.gz"

done
#

Run MultiQC on all FastQC outputs
apptainer run \
    --bind "$(pwd):$(pwd)" \
    "$MULTIQC_SIF" \
    multiqc "$FASTQC_DIR" -o "$(pwd)/multiqc_report"
BARCODE_DIR="/scratch/users/u15756/barcode_counts"
mkdir -p "$BARCODE_DIR"
bowtie2-build /scratch/users/u15756/hg38_ecoli.fa /scratch/users/u15756/hg38_ecoli_index
bowtie2-build /scratch/users/u15756/rhemac10_ecoli.fa /scratch/users/u15756/rhemac10_ecoli_index


OUTPUT_DIR="/scratch/users/u15756/output_bowtie/"
INPUT_DIR=""
INDEX="/scratch/users/u15756/hg38_ecoli_index/hg38_ecoli_index"

mkdir -p "$OUTPUT_DIR"

for sample_dir in "$INPUT_DIR"/*/; do
    sample_group=$(basename "$sample_dir")

    r1=$(find "$sample_dir" -name "*R1_trimmed.fastq.gz" | head -n 1)
    r2=$(find "$sample_dir" -name "*R2_trimmed.fastq.gz" | head -n 1)

    if [[ -z "${r1:-}" || -z "${r2:-}" ]]; then
        echo "Skipping $sample_group: missing reads"
        continue
    fi

    sample=$(basename "$r1")
    sample=${sample%%_R1_*}
    rm -f "${OUTPUT_DIR}/${sample}_sorted.bam"*
    rm -f "${OUTPUT_DIR}/${sample}_dedup.bam"*
    rm -f "${OUTPUT_DIR}/${sample}_tmp"*
    echo "Aligning $sample_group / $sample"

    bowtie2 \
        -p 8 \
        --very-sensitive \
        --dovetail \
        -x "$INDEX" \
        -1 "$r1" \
        -2 "$r2" \
        2> "${OUTPUT_DIR}/${sample}_bowtie2.log" \
    | samtools view -bS - \
    | samtools collate -u -O - \
    | samtools fixmate -m - - \
    | samtools sort -T "${OUTPUT_DIR}/${sample}_tmp" -o "${OUTPUT_DIR}/${sample}_sorted.bam" -

    samtools markdup \
        "${OUTPUT_DIR}/${sample}_sorted.bam" \
        "${OUTPUT_DIR}/${sample}_dedup.bam"

    samtools index "${OUTPUT_DIR}/${sample}_dedup.bam"
    echo "Counting barcodes for $sample_group / $sample"

    # R1
    echo "R1 barcode counts for $sample" > "${BARCODE_DIR}/${sample}_R1_counts.txt"
    for barcode in TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA
    do
        count=$(zgrep -c "$barcode" "${INPUT_DIR}/${sample_group}/${sample}_R1_trimmed.fastq.gz")
        echo -e "$barcode\t$count" >> "${BARCODE_DIR}/${sample}_R1_counts.txt"
    done

    # R2
    echo "R2 barcode counts for $sample" > "${BARCODE_DIR}/${sample}_R2_counts.txt"
    for barcode in TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA
    do
        count=$(zgrep -c "$barcode" "${INPUT_DIR}/${sample_group}/${sample}_R2_trimmed.fastq.gz")
        echo -e "$barcode\t$count" >> "${BARCODE_DIR}/${sample}_R2_counts.txt"
    done
done


INPUT_DIR="/scratch/users/u15756/p1841/p1841"
INDEX="/scratch/users/u15756/rhemac10_ecoli_index/rhemac10_ecoli_index"
BARCODE_DIR="/scratch/users/u15756/barcode_counts"

for sample_dir in "$INPUT_DIR"/*/; do
    sample_group=$(basename "$sample_dir")

    r1=$(find "$sample_dir" -name "*R1_001.fastq.gz" | head -n 1)
    r2=$(find "$sample_dir" -name "*R2_001.fastq.gz" | head -n 1)

    if [[ -z "${r1:-}" || -z "${r2:-}" ]]; then
        echo "Skipping $sample_group: missing reads"
        continue
    fi

    sample=$(basename "$r1")
    sample=${sample%%_R1_*}
    rm -f "${OUTPUT_DIR}/${sample}_sorted.bam"*
    rm -f "${OUTPUT_DIR}/${sample}_dedup.bam"*
    rm -f "${OUTPUT_DIR}/${sample}_tmp"*
    echo "Aligning $sample_group / $sample"
    bowtie2 \
        --very-sensitive \
        -p 8 \
        --dovetail \
        -x "$INDEX" \
        -1 "$r1" \
        -2 "$r2" \
        2> "${OUTPUT_DIR}/${sample}_bowtie2.log" \
    | samtools view -bS - \
    | samtools collate -u -O - \
    | samtools fixmate -m - - \
    | samtools sort -T "${OUTPUT_DIR}/${sample}_tmp" -o "${OUTPUT_DIR}/${sample}_sorted.bam" -

    samtools markdup \
        "${OUTPUT_DIR}/${sample}_sorted.bam" \
        "${OUTPUT_DIR}/${sample}_dedup.bam"

    samtools index "${OUTPUT_DIR}/${sample}_dedup.bam"

    echo "Counting barcodes for $sample_group / $sample"

#     # R1
    echo "R1 barcode counts for $sample" > "${BARCODE_DIR}/${sample}_R1_counts.txt"
    for barcode in TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA
    do
        count=$(zgrep -c "$barcode" "${INPUT_DIR}/${sample_group}/${sample}_R1_001.fastq.gz" || true)
        echo -e "$barcode\t$count" >> "${BARCODE_DIR}/${sample}_R1_counts.txt"
    done

    # R2
    echo "R2 barcode counts for $sample" > "${BARCODE_DIR}/${sample}_R2_counts.txt"
    for barcode in TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA
    do
        count=$(zgrep -c "$barcode" "${INPUT_DIR}/${sample_group}/${sample}_R2_001.fastq.gz" || true)
        echo -e "$barcode\t$count" >> "${BARCODE_DIR}/${sample}_R2_counts.txt"
    done
done

OUTPUT_DIR="/scratch/users/u15756/output_bowtie/"
shopt -s nullglob
shopt -s nullglob

for bam in "$OUTPUT_DIR"/*_dedup.bam; do
    sample=$(basename "$bam" _dedup.bam)
    sample_dir="${OUTPUT_DIR}/${sample}"

    mkdir -p "$sample_dir"

    mv "$bam" "$sample_dir/"
    mv "${bam}.bai" "$sample_dir/" 2>/dev/null || true

    bam="${sample_dir}/$(basename "$bam")"

    # idxstats
    idx="${sample_dir}/${sample}_idxstats.txt"
    samtools idxstats "$bam" > "$idx"
    
    # Extract E. coli contig
    ecoli_contig=$(awk '$1 ~ /U00096/ {print $1}' "$idx")
    if [[ -z "$ecoli_contig" ]]; then
        echo "No E. coli contig found for $sample"
        continue  
    fi

    # Split BAM
    samtools view -b "$bam" "$ecoli_contig" > "${sample_dir}/${sample}_ecoli.bam"

    samtools idxstats "$bam" \
    | awk '$1 != "*" && $1 !~ /U00096/ {print $1}' \
    | xargs samtools view -b "$bam" > "${sample_dir}/${sample}_host.bam"

    samtools index "${sample_dir}/${sample}_ecoli.bam"
    samtools index "${sample_dir}/${sample}_host.bam"

    # QC
    samtools flagstat "$bam" > "${sample_dir}/${sample}_flagstat.txt"
    
    # Summary
    awk '
    $1 != "*" { total += $3 }
    $1 ~ /U00096/ { bact = $3 }
    END {
        if (total > 0) {
            print "total_mapped=" total
            print "bacteria=" bact
            print "fraction=" bact/total
        } else {
            print "No reads found"
        }
    }' "$idx"

done



