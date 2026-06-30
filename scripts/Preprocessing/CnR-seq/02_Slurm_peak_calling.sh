#!/bin/bash
#SBATCH --job-name=CnR_peak_calling
#SBATCH --mail-user=c.bastidasbetancou@stud.uni-goettingen.de
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --time=48:00:00
#SBATCH --partition=medium
#SBATCH --output=CnR_peak_calling.out
#SBATCH --error=CnR_peak_calling.err


###350305 would be the median to use for E.coli spike normalization
# scaling factors:
# p1827sOVOL2-HA-plus_Heide_S49 1.15
# p1827sOVOL2-minus_Heide_S51 0.43
# p1827sOVOL2-plus_Heide_S50 0.99
# 
# p1827sZNF90-HA-minus_Heide_S46 1.18
# p1827sZNF90-HA-plus_Heide_S45 2.24
# p1827sZNF90-minus_Heide_S48 1.00
# p1827sZNF90-plus_Heide_S47 0.44


module load apptainer
module load gcc/14.2.0
module load bedtools2/2.31.1
module load samtools

###350305 would be the median to use for E.coli spike normalization
# scaling factors:
# p1841s0V0L2-HA-minus_Heide_S2	0.0146
# p1841s0V0L2-HA-plus_Heide_S1 0.0083
# p1841s0V0L2-minus_Heide_S4 0.0394
# p1841s0V0L2-plus_Heide_S3	0.0181

# p1827sZNF90-HA-minus_Heide_S46 1.18
# p1827sZNF90-HA-plus_Heide_S45 2.24
# p1827sZNF90-minus_Heide_S48 1.00
# p1827sZNF90-plus_Heide_S47 0.44

cd /scratch/users/u15756/output_bowtie/ovol2/

apptainer exec \
  --bind /scratch \
  /sw/container/bioinformatics/deeptools-3.5.6.sif \
  bamCoverage \
  -b /scratch/users/u15756/output_bowtie/ovol2/p1841s0V0L2-HA-minus_Heide_S2/p1841s0V0L2-HA-minus_Heide_S2_host.bam \
  -o /scratch/users/u15756/output_bowtie/ovol2/p1841s0V0L2-HA-minus_Heide_S2/p1841s0V0L2-HA-minus_Heide_S2_host.bw \
  --binSize 10 \
  --scaleFactor 1.37

apptainer exec \
  --bind /scratch \
  /sw/container/bioinformatics/deeptools-3.5.6.sif \
  bamCoverage \
  -b /scratch/users/u15756/output_bowtie/ovol2/p1841s0V0L2-HA-plus_Heide_S1/p1841s0V0L2-HA-plus_Heide_S1_host.bam \
  -o /scratch/users/u15756/output_bowtie/ovol2/p1841s0V0L2-HA-plus_Heide_S1/p1841s0V0L2-HA-plus_Heide_S1_host.bw \
  --binSize 10 \
  --scaleFactor 2.42

apptainer exec \
  --bind /scratch \
  /sw/container/bioinformatics/deeptools-3.5.6.sif \
  bamCoverage \
  -b /scratch/users/u15756/output_bowtie/ovol2/p1841s0V0L2-minus_Heide_S4/p1841s0V0L2-minus_Heide_S4_host.bam \
  -o /scratch/users/u15756/output_bowtie/ovol2/p1841s0V0L2-minus_Heide_S4/p1841s0V0L2-minus_Heide_S4_host.bw \
  --binSize 10 \
  --scaleFactor 0.51

apptainer exec \
  --bind /scratch \
  /sw/container/bioinformatics/deeptools-3.5.6.sif \
  bamCoverage \
  -b /scratch/users/u15756/output_bowtie/ovol2/p1841s0V0L2-plus_Heide_S3/p1841s0V0L2-plus_Heide_S3_host.bam \
  -o /scratch/users/u15756/output_bowtie/ovol2/p1841s0V0L2-plus_Heide_S3/p1841s0V0L2-plus_Heide_S3_host.bw \
  --binSize 10 \
  --scaleFactor 1.11
########################
apptainer exec \
  --bind /scratch \
  /sw/container/bioinformatics/deeptools-3.5.6.sif \
  bamCoverage \
  -b /scratch/users/u15756/output_bowtie/p1827sZNF90-HA-minus_Heide_S46/p1827sZNF90-HA-minus_Heide_S46_host.bam \
  -o  /scratch/users/u15756/output_bowtie/p1827sZNF90-HA-minus_Heide_S46/p1827sZNF90-HA-minus_Heide_S46_host.bw \
  --binSize 10 \
  --scaleFactor 1.18
  
apptainer exec \
  --bind /scratch \
  /sw/container/bioinformatics/deeptools-3.5.6.sif \
  bamCoverage \
  -b /scratch/users/u15756/output_bowtie/p1827sZNF90-HA-plus_Heide_S45/p1827sZNF90-HA-plus_Heide_S45_host.bam \
  -o /scratch/users/u15756/output_bowtie/p1827sZNF90-HA-plus_Heide_S45/p1827sZNF90-HA-plus_Heide_S45_host.bw \
  --binSize 10 \
  --scaleFactor 2.24
  
apptainer exec \
  --bind /scratch \
  /sw/container/bioinformatics/deeptools-3.5.6.sif \
  bamCoverage \
  -b /scratch/users/u15756/output_bowtie/p1827sZNF90-minus_Heide_S48/p1827sZNF90-minus_Heide_S48_host.bam \
  -o /scratch/users/u15756/output_bowtie/p1827sZNF90-minus_Heide_S48/p1827sZNF90-minus_Heide_S48_host.bw \
  --binSize 10 \
  --scaleFactor 1.00

apptainer exec \
  --bind /scratch \
  /sw/container/bioinformatics/deeptools-3.5.6.sif \
  bamCoverage \
  -b /scratch/users/u15756/output_bowtie/p1827sZNF90-plus_Heide_S47/p1827sZNF90-plus_Heide_S47_host.bam \
  -o /scratch/users/u15756/output_bowtie/p1827sZNF90-plus_Heide_S47/p1827sZNF90-plus_Heide_S47_host.bw \
  --binSize 10 \
  --scaleFactor 0.44

###HEATMAP PLOTS

awk 'BEGIN{OFS="\t"} $3=="transcript" {

    match($0, /gene_id "([^"]+)"/, a)
    gene=a[1]

    if (gene=="") next
    if ($7!="+" && $7!="-") next

    if ($7=="+") {
        tss=$4
    } else {
        tss=$5
    }

    if (!(gene in best)) {
        best[gene]=tss
        chr[gene]=$1
        strand[gene]=$7
    }

    # keep most upstream / 5-prime TSS per gene
    if ($7=="+" && tss < best[gene]) {
        best[gene]=tss
        chr[gene]=$1
        strand[gene]=$7
    }

    if ($7=="-" && tss > best[gene]) {
        best[gene]=tss
        chr[gene]=$1
        strand[gene]=$7
    }
}

END {
    for (g in best) {
        # Convert GTF 1-based TSS to BED 0-based half-open interval
        start=best[g]-1
        end=best[g]

        if (start<0) start=0

        print chr[g], start, end, g, ".", strand[g]
    }
}' rheMac10.ncbiRefSeq.gtf > promoters_rhesus_unique.bed

awk 'BEGIN{OFS="\t"} $3=="transcript" {

    match($0, /gene_id "([^"]+)"/, a)
    gene=a[1]

    if (gene=="") next
    if ($7!="+" && $7!="-") next

    if ($7=="+") {
        tss=$4
    } else {
        tss=$5
    }

    if (!(gene in best)) {
        best[gene]=tss
        chr[gene]=$1
        strand[gene]=$7
    }

    # keep most upstream / 5-prime TSS per gene
    if ($7=="+" && tss < best[gene]) {
        best[gene]=tss
        chr[gene]=$1
        strand[gene]=$7
    }

    if ($7=="-" && tss > best[gene]) {
        best[gene]=tss
        chr[gene]=$1
        strand[gene]=$7
    }
}

END {
    for (g in best) {
        # Convert GTF 1-based TSS to BED 0-based half-open interval
        start=best[g]-1
        end=best[g]

        if (start<0) start=0

        print chr[g], start, end, g, ".", strand[g]
    }
}' ./Downloads/annotation_human.gtf > promoters_human_unique.bed


computeMatrix reference-point \
  -S p1827sZNF90-HA-plus_Heide_S45_host.bw p1827sZNF90-minus_Heide_S48_host.bw \
  -R ./promoters_human_unique.bed \
  --referencePoint TSS \
  -b 2000 -a 2000 \
  --skipZeros \
  --missingDataAsZero \
  --binSize 50 \
  -out promoter_matrix_human.gz
  
computeMatrix reference-point \
  -S p1841s0V0L2-HA-plus_Heide_S1_host.bw p1841s0V0L2-HA-minus_Heide_S2_host.bw \
  -R ./promoters_rhesus_unique.bed \
  --referencePoint TSS \
  -b 2000 -a 2000 \
  --skipZeros \
  --missingDataAsZero \
  --binSize 50 \
  -out promoter_matrix_rhesus.gz
  
plotHeatmap \
  -m promoter_matrix_rhesus.gz \
  -out promoter_heatmap_rhesus.png
  
plotHeatmap \
  -m promoter_matrix_human.gz \
  -out promoter_heatmap_human.png
  
plotProfile \
  -m promoter_matrix_human.gz \
  -out promoter_profile_human.png
plotProfile \
  -m promoter_matrix_rhesus.gz \
  -out promoter_profile_rhesus.png
  
#####merging untagged assays to run another MACS with those
samtools merge -o merged_control_ZNF90.bam ./p1827sZNF90-minus_Heide_S48/p1827sZNF90-minus_Heide_S48_host.bam ./p1827sZNF90-plus_Heide_S47/p1827sZNF90-plus_Heide_S47_host.bam

#cd /scratch/users/u15756/output_bowtie/ovol2/
samtools merge -o merged_control_OVOL2.bam ./p1841s0V0L2-minus_Heide_S4/p1841s0V0L2-minus_Heide_S4_host.bam ./p1841s0V0L2-plus_Heide_S3/p1841s0V0L2-plus_Heide_S3_host.bam

############now MACS3 with IgG controls
#####or rhesys macaque genome and put control sample!
apptainer exec \
  --bind /scratch \
  /sw/container/bioinformatics/macs3-3.0.3.sif \
  macs3 callpeak \
  -t /scratch/users/u15756/output_bowtie/p1827sZNF90-HA-plus_Heide_S45/p1827sZNF90-HA-plus_Heide_S45_host.bam \
  -c /scratch/users/u15756/output_bowtie/p1827sZNF90-HA-minus_Heide_S46/p1827sZNF90-HA-minus_Heide_S46_host.bam \
  -f BAMPE \
  -g hs \
  -n ZNF90_test \
  --outdir /scratch/users/u15756/output_bowtie/macs3_test/ZNF90 \
  --nomodel \
  --keep-dup all \
  -q 0.05
  
apptainer exec \
  --bind /scratch \
  /sw/container/bioinformatics/macs3-3.0.3.sif \
  macs3 callpeak \
  -t /scratch/users/u15756/output_bowtie/p1827sZNF90-HA-plus_Heide_S45/p1827sZNF90-HA-plus_Heide_S45_host.bam \
  -c /scratch/users/u15756/output_bowtie/merged_control_ZNF90.bam \
  -f BAMPE \
  -g hs \
  -n ZNF90_test_untagged \
  --outdir /scratch/users/u15756/output_bowtie/macs3_test/ZNF90_untagged \
  --nomodel \
  --keep-dup all \
  -q 0.05
    
  
######################### OVOL2-wait until all samples have been sequenced  
  
apptainer exec \
  --bind /scratch \
  /sw/container/bioinformatics/macs3-3.0.3.sif \
  macs3 callpeak \
  -t /scratch/users/u15756/output_bowtie/ovol2/p1841s0V0L2-HA-plus_Heide_S1/p1841s0V0L2-HA-plus_Heide_S1_host.bam \
  -c /scratch/users/u15756/output_bowtie/ovol2/p1841s0V0L2-HA-minus_Heide_S2/p1841s0V0L2-HA-minus_Heide_S2_host.bam \
  -f BAMPE \
  -g 3.1e9 \
  -n OVOL2_negative \
  --outdir /scratch/users/u15756/output_bowtie/macs3_test/OVOL2 \
  --nomodel \
  --keep-dup all \
  -q 0.05
  
apptainer exec \
  --bind /scratch \
  /sw/container/bioinformatics/macs3-3.0.3.sif \
  macs3 callpeak \
  -t /scratch/users/u15756/output_bowtie/ovol2/p1841s0V0L2-HA-plus_Heide_S1/p1841s0V0L2-HA-plus_Heide_S1_host.bam \
  -c /scratch/users/u15756/output_bowtie/ovol2/merged_control_OVOL2.bam \
  -f BAMPE \
  -g 3.1e9 \
  -n OVOL2_merge\
  --outdir /scratch/users/u15756/output_bowtie/macs3_test/OVOL2 \
  --nomodel \
  --keep-dup all \
  -q 0.05
  
  
  
################# estimate overlap between controls, ideally should be high
###ZNF90
bedtools intersect -u -a ./ZNF90/ZNF90_test_peaks.narrowPeak -b ./ZNF90_untagged/ZNF90_test_untagged_peaks.narrowPeak > ZNF90_shared_peaks_bedtools_intersect.narrowPeak

bedtools intersect -v -a ./ZNF90/ZNF90_test_peaks.narrowPeak -b ./ZNF90_untagged/ZNF90_test_untagged_peaks.narrowPeak > ZNF90_HA_vs_tagged_IgG.narrowPeak
bedtools intersect -v -a ./ZNF90_untagged/ZNF90_test_untagged_peaks.narrowPeak -b ./ZNF90/ZNF90_test_peaks.narrowPeak > ZNF90_HA_vs_untagged.narrowPeak
####OVOL2
cd /scratch/users/u15756/output_bowtie/macs3_test/
bedtools intersect -u -a ./OVOL2/OVOL2_negative_peaks.narrowPeak -b ./OVOL2/OVOL2_merge_peaks.narrowPeak > OVOL2_shared_peaks_bedtools_intersect.narrowPeak

bedtools intersect -v -a ./OVOL2/OVOL2_negative_peaks.narrowPeak -b ./OVOL2/OVOL2_merge_peaks.narrowPeak > OVOL2_HA_vs_tagged_IgG.narrowPeak
bedtools intersect -v -a ./OVOL2/OVOL2_merge_peaks.narrowPeak -b ./OVOL2/OVOL2_negative_peaks.narrowPeak > OVOL2_HA_vs_untagged.narrowPeak




