#!/bin/bash
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=trinity_run
#SBATCH --mem=350G
#SBATCH --mail-user=c.bastidasbetancou@stud.uni-goettingen.de
#SBATCH --mail-type=BEGIN,END
#SBATCH --account=all
#SBATCH --time=48:00:00
#SBATCH --output=qctrim_run_liftoff_human.out
#SBATCH --error=qctrim_run_liftoff_human.err

# Load Trinity module and dependencies (replace 'trinity/2.11.0' with the correct version if needed)

cd /scratch/users/cbastid/reads
module load rev/21.12
module load multiqc/1.13
module load trimgalore/0.6.10
module load cutadapt
module load fastqc
module load star
parent_folder="/scratch/users/cbastid/rhesus"

cd /scratch/users/cbastid/reads


extract_srr() {
    filename="$1"
    # Extract SRR identifier from the filename
    srr=$(basename "$filename" | cut -d '_' -f 1)
    echo "$srr"
}

process_folder() {
    folder="$1"
    echo "Processing folder: $folder"
    cd "$folder" || return
    mkdir -p /scratch/users/cbastid/STAR_output/rhesus/"$subfolder"
    # Iterate over each FASTQ file
    if [ "$subfolder" = "johnson" ]; then
	for file in *_1.fastq.gz; do
        	[ -e "$file" ] || continue  # Skip if no matching file found
        	srr=$(extract_srr "$file")
		#gunzip "$srr"_1.fastq.gz
		#gunzip  "$srr"_2.fastq.gz
        	# Run Trim Galore! for read trimming
		#cutadapt -q 20 --minimum-length 1 -o "$srr"_1_trimmed.fastq.gz -p "$srr"_2_trimmed.fastq.gz "$srr"_1.fastq.gz "$srr"_2.fastq.gz
        	trim_galore --paired -q 20 --fastqc --gzip "$srr"_1.fastq.gz "$srr"_2.fastq.gz
		#fastqc --outdir ../results "$srr"_1_trimmed.fastq.gz "$srr"_2_trimmed.fastq.gz



		gunzip -c "$srr"_1_val_1.fq.gz > /scratch/users/cbastid/rhesus/"$subfolder"/"$srr"_1_val_1.fq
		gunzip -c "$srr"_2_val_2.fq.gz > /scratch/users/cbastid/rhesus/"$subfolder"/"$srr"_2_val_2.fq
		rm -rf /scratch/users/cbastid/star_tmp/

		 Run STAR with the decompressed paired-end read files
 		STAR \
 		--genomeDir /scratch/users/cbastid/index_star_human \
 		--readFilesIn /scratch/users/cbastid/rhesus/"$subfolder"/"$srr"_1_val_1.fq /scratch/users/cbastid/rhesus/"$subfolder"/"$srr"_2_val_2.fq \
 		--quantMode GeneCounts \
 		--outSAMtype None \
 		--runThreadN $SLURM_NTASKS \
 		--outFileNamePrefix /scratch/users/cbastid/STAR_output/rhesus/"$subfolder"/"$srr"_output_ \
 		--outTmpDir /scratch/users/cbastid/star_tmp

		rm /scratch/users/cbastid/rhesus/"$subfolder"/"$srr"_1_val_1.fq
		rm /scratch/users/cbastid/rhesus/"$subfolder"/"$srr"_2_val_2.fq
	    done
    else
	for file in *.fastq.gz; do
		#echo "$folder"
		#extract_SRR_SE "$file"
        	#cutadapt -q 20 -minimum-length 1 -o "$file"_trimmed.fastq.gz "$file"
		trim_galore -q 20 --fastqc --gzip "$file"
		#fastqc --outdir ../results "$file"_trimmed.fastq.gz
		SRR=$(basename "$file" | cut -d'.' -f1)
		gunzip -c "$SRR"_trimmed.fq.gz > /scratch/users/cbastid/rhesus/"$subfolder"/"$SRR"_trimmed.fq
		STAR \
                --genomeDir /scratch/users/cbastid/index_star_human \
                --readFilesIn /scratch/users/cbastid/rhesus/"$subfolder"/"$SRR"_trimmed.fq \
                --quantMode GeneCounts \
                --outSAMtype None \
                --runThreadN $SLURM_NTASKS \
                --outFileNamePrefix /scratch/users/cbastid/STAR_output/rhesus/"$subfolder"/"$SRR"_output_ \
                --outTmpDir /scratch/users/cbastid/star_tmp
		
		rm -r /scratch/users/cbastid/star_tmp/
		rm /scratch/users/cbastid/rhesus/"$subfolder"/"$SRR"_trimmed.fastq
	done
    fi
    echo "Finished processing folder: $folder"
}



# Iterate over each folder
for folder in "$parent_folder"/{johnson,florio,fietz}; do
    if [ -d "$folder" ]; then
        subfolder=$(basename "$folder")
	process_folder "$folder"

	#subfolder=$(basename "$folder")
    else
        echo "Folder $folder does not exist."
    fi
done
