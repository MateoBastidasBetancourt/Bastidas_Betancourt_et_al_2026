#!/bin/bash
#SBATCH --job-name=GenomeIndexing
#SBATCH --mail-user=c.bastidasbetancou@stud.uni-goettingen.de
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=2
#SBATCH --ntasks=48
#SBATCH --mem=320gb   # Adjust memory as needed based on your dataset
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --account=all
#SBATCH --output=genome_indexing_rhesus_liftoff.out
#SBATCH --error=genome_indexing_rhesus_liftoff.err

echo Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID
echo Running on host `hostname`
echo Job started at `date +"%T %a %d %b %Y"`
echo Directory is `pwd`
echo Using $SLURM_NTASKS processors across $SLURM_NNODES nodes

module load star/2.7.3a
GENOME_RHESUS=~/RNAseq/STAR/rheMac10.fa.gz

ANNOTATION_RHESUS=~/RNAseq/STAR/rheMac10_liftoff_GRCh38.p13_RefSeq.gtf.gz


## Create separate output directories for each species
mkdir -p /scratch/users/cbastid/index_rhesus_liftoff

# Decompress genome reference files for human and marmoset
gunzip -c $GENOME_RHESUS > /scratch/users/cbastid/genomes_annotations/rheMac10.fa
# Decompress annotation files for human and marmoset
gunzip -c $ANNOTATION_RHESUS > /scratch/users/cbastid/genomes_annotations/rheMac10_liftoff_GRCh38.p13_RefSeq.gtf


# Generate genome index for human
STAR \
  --runMode genomeGenerate \
  --genomeDir /scratch/users/cbastid/index_rhesus \
  --genomeFastaFiles /scratch/users/cbastid/genomes_annotations/rheMac10.fa \
  --sjdbGTFfile /scratch/users/cbastid/genomes_annotations/rheMac10_liftoff_GRCh38.p13_RefSeq.gtf \
  --sjdbOverhang 100 \
  --runThreadN $SLURM_NTASKS \
  --limitGenomeGenerateRAM 31000000000 \
  --outTmpDir /scratch/users/cbastid/STAR \

# Delete the decompressed files to save space
rm /scratch/users/cbastid/genomes_annotations/rheMac10.fa
rm /scratch/users/cbastid/genomes_annotations/rheMac10_liftoff_GRCh38.p13_RefSeq.gtf

GENOME_HUMAN=~/RNAseq/STAR/fasta/reference_human.fa.gz

## PATH to the GTF annotation files for human and marmoset
ANNOTATION_HUMAN=~/RNAseq/STAR/fasta/annotation_human.gtf.gz

rm -rf /scratch/users/cbastid/STAR/

## Create separate output directories for each species
mkdir -p /scratch/users/cbastid/index_star_human_149

# Decompress genome reference files for human and marmoset
gunzip -c $GENOME_HUMAN > ~/RNAseq/STAR/fasta/reference_human.fa
# Decompress annotation files for human and marmoset
gunzip -c $ANNOTATION_HUMAN > ~/RNAseq/STAR/fasta/annotation_human.gtf


# Generate genome index for human
STAR \
  --runMode genomeGenerate \
  --genomeDir /scratch/users/cbastid/index_star_human_149 \
  --genomeFastaFiles ~/RNAseq/STAR/fasta/reference_human.fa \
  --sjdbGTFfile ~/RNAseq/STAR/fasta/annotation_human.gtf \
  --sjdbOverhang 149 \
  --runThreadN $SLURM_NTASKS \
  --limitGenomeGenerateRAM 39500000000 \
  --outTmpDir /scratch/users/cbastid/STAR \

# Delete the decompressed files to save space
rm ~/RNAseq/STAR/fasta/reference_human.fa
rm ~/RNAseq/STAR/fasta/annotation_human.gtf


rm -rf /scratch/users/cbastid/STAR/

