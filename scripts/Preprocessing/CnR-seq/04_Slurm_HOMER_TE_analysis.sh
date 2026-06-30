#!/bin/bash
#SBATCH --job-name=CnR_HOMER
#SBATCH --mail-user=c.bastidasbetancou@stud.uni-goettingen.de
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --time=48:00:00
#SBATCH --partition=medium
#SBATCH --output=CnR_HOMER.out
#SBATCH --error=CnR_HOMER.err
module load gcc/14.2.0
module load apptainer
#module load gcc/11.5.0 
module load openmpi/4.1.7
module load repeatmasker

module load bedtools2/2.31.1

################# installing genomes

apptainer exec \
--bind $HOME:$HOME \
--bind /scratch/users/u15756/homer_data:/opt/conda/share/homer/data \
/sw/container/bioinformatics/homer-latest.sif \
perl /opt/conda/share/homer/configureHomer.pl -install rheMac10

apptainer exec \
--bind $HOME:$HOME \
--bind /scratch/users/u15756/homer_data:/opt/conda/share/homer/data \
/sw/container/bioinformatics/homer-latest.sif \
perl /opt/conda/share/homer/configureHomer.pl -install hg38


####convert input files to bed3 format and deal with chromosome naming
awk 'BEGIN{OFS="\t"} {print "chr"$1,$2,$3}' ~/macs3_test/ZNF90_shared_peaks_bedtools_intersect.narrowPeak > ~/macs3_test/ZNF90_shared_peaks_bedtools_intersect.bed
awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' \
~/macs3_test/OVOL2_shared_peaks_bedtools_intersect.narrowPeak \
> ~/macs3_test/OVOL2_shared_peaks_bedtools_intersect.bed

 
###ZNF90 HOMER
apptainer exec \
  --bind $HOME:$HOME \
  --bind /scratch/users/u15756:/scratch/users/u15756 \
  /sw/container/bioinformatics/homer-latest.sif \
  /opt/conda/share/homer/bin/findMotifsGenome.pl \
  ~/macs3_test/ZNF90_shared_peaks_bedtools_intersect.bed \
  /scratch/users/u15756/homer_data/genomes/hg38/genome.fa \
  ~/macs3_test/homer_motif_output_short \
  -size 200 \
  -len 8,10,12 
  

apptainer exec \
  --bind $HOME:$HOME \
  --bind /scratch/users/u15756:/scratch/users/u15756 \
  /sw/container/bioinformatics/homer-latest.sif \
  /opt/conda/share/homer/bin/annotatePeaks.pl \
  ~/macs3_test/ZNF90_shared_peaks_bedtools_intersect.bed \
  /scratch/users/u15756/homer_data/genomes/hg38/genome.fa \
  -gtf ~/hg38.knownGene.gtf -gid \
  > ~/macs3_test/ZNF90_homer_annotation.txt 
  
####OVOL2 HOMER
apptainer exec \
  --bind $HOME:$HOME \
  --bind /scratch/users/u15756:/scratch/users/u15756 \
  /sw/container/bioinformatics/homer-latest.sif \
  /opt/conda/share/homer/bin/findMotifsGenome.pl \
  ~/macs3_test/OVOL2_shared_peaks_bedtools_intersect.bed \
  /scratch/users/u15756/homer_data/genomes/rheMac10/genome.fa \
  ~/macs3_test/homer_motif_output_OVOL2 \
  -size 200 \
  -len 8,10,12 
  

apptainer exec \
  --bind $HOME:$HOME \
  --bind /scratch/users/u15756:/scratch/users/u15756 \
  /sw/container/bioinformatics/homer-latest.sif \
  /opt/conda/share/homer/bin/annotatePeaks.pl \
  ~/macs3_test/OVOL2_shared_peaks_bedtools_intersect.bed \
  /scratch/users/u15756/homer_data/genomes/rheMac10/genome.fa \
  -gtf ~/rheMac10.ncbiRefSeq.gtf -gid \
  > ~/macs3_test/OVOL2_homer_annotation.txt 
#### repeatmasker after uploading the UCSC repeat files for hg38 and rheMac10

bedtools intersect -a ~/macs3_test/ZNF90_shared_peaks_bedtools_intersect.bed -b ~/repeats_human.bed -wa -wb > ~/macs3_test/peaks_in_repeats_ZNF90.bed
bedtools intersect \
-a ~/macs3_test/OVOL2_shared_peaks_bedtools_intersect.bed \
-b ~/repeats_rhesus.bed \
-wa -wb \
> ~/macs3_test/OVOL2_output_HOMER_RepeatMasker/peaks_in_repeats_OVOL2.bed

####
##count
cut -f1-3 ~/macs3_test/peaks_in_repeats_ZNF90.bed | sort | uniq | wc -l
#708/6599
cut -f1-3 ~/macs3_test/peaks_in_repeats_OVOL2.bed | sort | uniq | wc -l
#5344/91983

#ZNF90 repetitive element abundance
cut -f7 ~/macs3_test/ZNF90_output_HOMER_RepeatMasker/peaks_in_repeats_ZNF90.bed | \
awk '
{
  # -------------------------
  # NORMALIZATION
  # -------------------------
  gsub(/_int$/, "")
  gsub(/-int$/, "")
  gsub(/_I-int$/, "")
  gsub(/\r/, "")

  print
}
' | sort | uniq -c | \
awk '
{
  count=$1
  name=$2

  # -------------------------
  # CLASSIFICATION
  # -------------------------
  if (name ~ /^Alu/ || name ~ /^FLAM/ || name ~ /^FRAM/)
    class="Alu"

  else if (name ~ /^MIR/)
    class="MIR"

  else if (name ~ /^L1/)
    class="L1"

  else if (name ~ /^L2/)
    class="L2"

  else if (name ~ /^L3/)
    class="L3"

  else if (name ~ /^HERV/ || name ~ /^HUERS/ || name ~ /^ERV/ || name ~ /^ERVL/ ||
           name ~ /^LTR/ || name ~ /^MLT/ || name ~ /^MST/ || name ~ /^THE1/ ||
           name ~ /^HAL1/ || name ~ /^LOR1/ ||
           name ~ /^MamGyp/ || name ~ /^MamRep/ || name ~ /^MamSINE/ ||
           name ~ /^MADE/ || name ~ /^PRIMAX/)
    class="LTR_ERV"

  else if (name ~ /^Charlie/ || name ~ /^Tigger/ || name ~ /^Zaphod/ ||
           name ~ /^MER/ || name ~ /^PABL/)
    class="DNA_transposon"

  else if (name ~ /^\(/ || name ~ /rich/)
    class="Simple_repeat"

  else
    class="Other"

  # -------------------------
  # AGE GROUPING
  # -------------------------
  if (name ~ /^SVA/ || name ~ /^HERVK/ || name ~ /^L1HS/ || name ~ /^L1PA1/)
    age="young"

  else if (name ~ /^AluY/ || name ~ /^L1PA/)
    age="young"

  else if (name ~ /^AluS/ || name ~ /^AluJ/ || name ~ /^L1PB/)
    age="intermediate"

  else if (name ~ /^MIR/ || name ~ /^L2/ || name ~ /^L3/ ||
           name ~ /^L1MA/ || name ~ /^L1MB/ || name ~ /^L1ME/ ||
           name ~ /^L1M/ || name ~ /^L1P/ ||
           name ~ /^MER/ || name ~ /^THE1/ || name ~ /^MST/ ||
           name ~ /^MLT/ || name ~ /^ERV/ || name ~ /^ERVL/ ||
           name ~ /^HERV/ || name ~ /^HUERS/ ||
           name ~ /^FLAM/ || name ~ /^FRAM/ || name ~ /^HAL1/)
    age="old"

  else if (name ~ /^\(/ || name ~ /rich/ || name ~ /^ALR/)
    age="old"

  else
    age="unknown"

  # -------------------------
  # SUMMARIES
  # -------------------------
  class_sum[class] += count
  age_sum[age] += count
  combo_sum[class "|" age] += count
}

END {
  print "=== Class summary ==="
  for (c in class_sum)
    print c, class_sum[c]

  print "\n=== Age summary ==="
  for (a in age_sum)
    print a, age_sum[a]

  print "\n=== Class + Age breakdown ==="
  for (k in combo_sum)
    print k, combo_sum[k]
}
'
# === Class summary ===
# LTR_ERV 159
# Alu 225
# MIR 96
# L1 157
# L2 93
# L3 8
# Other 21
# Simple_repeat 52
# DNA_transposon 80
# 
# === Age summary ===
# young 86
# intermediate 199
# old 527
# unknown 79
# 
# === Class + Age breakdown ===
# Other|young 4
# L2|old 93
# L1|young 53
# Alu|young 27
# Other|old 7
# Alu|old 5
# LTR_ERV|old 106
# DNA_transposon|old 62
# L3|old 8
# L1|intermediate 6
# L1|old 98
# Alu|intermediate 193
# LTR_ERV|unknown 51
# Other|unknown 10
# DNA_transposon|unknown 18
# LTR_ERV|young 2
# MIR|old 96
# Simple_repeat|old 52

cut -f4 repeats_human.bed | sort | uniq -c > genome_TE_counts.txt ### specific TEs
cut -f4 repeats_human.bed | \
awk '
{
  gsub(/_int$/, "")
  gsub(/-int$/, "")
  gsub(/_I-int$/, "")
  gsub(/\r/, "")
  print
}
' | sort | uniq -c | \
awk '
{
  count=$1
  name=$2

  # CLASSIFICATION
  if (name ~ /^Alu/ || name ~ /^FLAM/ || name ~ /^FRAM/)
    class="Alu"
  else if (name ~ /^MIR/)
    class="MIR"
  else if (name ~ /^L1/)
    class="L1"
  else if (name ~ /^L2/)
    class="L2"
  else if (name ~ /^L3/)
    class="L3"
  else if (name ~ /^HERV/ || name ~ /^HUERS/ || name ~ /^ERV/ || name ~ /^ERVL/ ||
           name ~ /^LTR/ || name ~ /^MLT/ || name ~ /^MST/ || name ~ /^THE1/ ||
           name ~ /^HAL1/ || name ~ /^LOR1/ ||
           name ~ /^MamGyp/ || name ~ /^MamRep/ || name ~ /^MamSINE/ ||
           name ~ /^MADE/ || name ~ /^PRIMAX/)
    class="LTR_ERV"
  else if (name ~ /^Charlie/ || name ~ /^Tigger/ || name ~ /^Zaphod/ ||
           name ~ /^MER/ || name ~ /^PABL/)
    class="DNA_transposon"
  else if (name ~ /^\(/ || name ~ /rich/)
    class="Simple_repeat"
  else
    class="Other"

  class_sum[class] += count
}

END {
  for (c in class_sum)
    print c, class_sum[c]
}
'
#numbers from repeatMasker database
# LTR_ERV 111131
# Alu 193693
# MIR 109354
# L1 149358
# L2 86173
# L3 10664
# Other 22131
# Simple_repeat 129960
# DNA_transposon 90497


######Alu and HERV/LTR are significantly enriched: now look inside those groups:
cut -f4 repeats_human.bed | \
awk '
{
  # -------------------------
  # CLEANING
  # -------------------------
  gsub(/_int$/, "")
  gsub(/-int$/, "")
  gsub(/_I-int$/, "")
  gsub(/\r/, "")
  name=$0

  # -------------------------
  # MAIN CLASS
  # -------------------------
  if (name ~ /^Alu/ || name ~ /^FLAM/ || name ~ /^FRAM/)
    class="Alu"

  else if (name ~ /^MIR/)
    class="MIR"

  else if (name ~ /^L1/)
    class="L1"

  else if (name ~ /^L2/)
    class="L2"

  else if (name ~ /^L3/)
    class="L3"

  else if (name ~ /^HERV/ || name ~ /^HUERS/ || name ~ /^ERV/ || name ~ /^ERVL/ ||
           name ~ /^LTR/ || name ~ /^MLT/ || name ~ /^MST/ || name ~ /^THE1/ ||
           name ~ /^HAL1/ || name ~ /^LOR1/ ||
           name ~ /^MamGyp/ || name ~ /^MamRep/ || name ~ /^MamSINE/ ||
           name ~ /^MADE/ || name ~ /^PRIMAX/)
    class="LTR_ERV"

  else if (name ~ /^Charlie/ || name ~ /^Tigger/ || name ~ /^Zaphod/ ||
           name ~ /^MER/ || name ~ /^PABL/)
    class="DNA_transposon"

  else if (name ~ /^\(/ || name ~ /rich/)
    class="Simple_repeat"

  else
    class="Other"

  # -------------------------
  # SUBCLASS (KEY PART 🔥)
  # -------------------------

  subclass="Other"

  if (class=="Alu") {
    if (name ~ /^AluY/)
      subclass="AluY"
    else if (name ~ /^AluS/)
      subclass="AluS"
    else if (name ~ /^AluJ/)
      subclass="AluJ"
    else if (name ~ /^FLAM/)
      subclass="FLAM"
    else if (name ~ /^FRAM/)
      subclass="FRAM"
    else
      subclass="Alu_other"
  }

  else if (class=="LTR_ERV") {
    if (name ~ /^HERVK/)
      subclass="HERVK"
    else if (name ~ /^HERV/)
      subclass="HERV_other"
    else if (name ~ /^LTR/)
      subclass="LTR"
    else if (name ~ /^THE1/)
      subclass="THE1"
    else if (name ~ /^MLT/)
      subclass="MLT"
    else if (name ~ /^MST/)
      subclass="MST"
    else if (name ~ /^ERV/)
      subclass="ERV_other"
    else
      subclass="LTR_misc"
  }

  else if (class=="DNA_transposon") {
    if (name ~ /^MER/)
      subclass="MER"
    else if (name ~ /^Charlie/)
      subclass="Charlie"
    else if (name ~ /^Tigger/)
      subclass="Tigger"
    else if (name ~ /^Zaphod/)
      subclass="Zaphod"
    else if (name ~ /^PABL/)
      subclass="PABL"
    else
      subclass="DNA_other"
  }

  else {
    subclass=class
  }

  # -------------------------
  # SUMMARIES
  # -------------------------
  class_sum[class]++
  subclass_sum[class "|" subclass]++
}

END {
  print "=== Class summary ==="
  for (c in class_sum)
    print c, class_sum[c]

  print "\n=== Subclass breakdown ==="
  for (k in subclass_sum)
    print k, subclass_sum[k]
}
'

##
# LTR_ERV|MST 8054
# LTR_ERV|LTR_misc 15584
# DNA_transposon|PABL 243
# Alu|FLAM 6711
# Simple_repeat|Simple_repeat 129960
# Alu|AluJ 51715
# L2|L2 86173
# DNA_transposon|Zaphod 782
# Alu|AluY 22927
# DNA_transposon|Tigger 14587
# L3|L3 10664
# LTR_ERV|THE1 9253
# Alu|FRAM 1159
# Alu|Alu_other 698
# Other|Other 22131
# LTR_ERV|MLT 44567
# LTR_ERV|ERV_other 3565
# LTR_ERV|LTR 27040
# LTR_ERV|HERV_other 2744
# MIR|MIR 109354
# Alu|AluS 110483
# LTR_ERV|HERVK 324
# DNA_transposon|MER 63274
# L1|L1 149358
# DNA_transposon|Charlie 11611

######now subcategories in my dataset:
cut -f7 ~/macs3_test/ZNF90_output_HOMER_RepeatMasker/peaks_in_repeats_ZNF90.bed | \
awk '
{
  # -------------------------
  # NORMALIZATION
  # -------------------------
  gsub(/_int$/, "")
  gsub(/-int$/, "")
  gsub(/_I-int$/, "")
  gsub(/\r/, "")

  print
}
' | sort | uniq -c | \
awk '
{
  count=$1
  name=$2

  # -------------------------
  # CLASSIFICATION
  # -------------------------
  if (name ~ /^Alu/ || name ~ /^FLAM/ || name ~ /^FRAM/)
    class="Alu"

  else if (name ~ /^MIR/)
    class="MIR"

  else if (name ~ /^L1/)
    class="L1"

  else if (name ~ /^L2/)
    class="L2"

  else if (name ~ /^L3/)
    class="L3"

  else if (name ~ /^HERV/ || name ~ /^HUERS/ || name ~ /^ERV/ || name ~ /^ERVL/ ||
           name ~ /^LTR/ || name ~ /^MLT/ || name ~ /^MST/ || name ~ /^THE1/ ||
           name ~ /^HAL1/ || name ~ /^LOR1/ ||
           name ~ /^MamGyp/ || name ~ /^MamRep/ || name ~ /^MamSINE/ ||
           name ~ /^MADE/ || name ~ /^PRIMAX/)
    class="LTR_ERV"

  else if (name ~ /^Charlie/ || name ~ /^Tigger/ || name ~ /^Zaphod/ ||
           name ~ /^MER/ || name ~ /^PABL/)
    class="DNA_transposon"

  else if (name ~ /^\(/ || name ~ /rich/)
    class="Simple_repeat"

  else
    class="Other"

  # -------------------------
  # SUBCLASS (NEW 🔥)
  # -------------------------
  if (class=="Alu") {
    if (name ~ /^AluY/)
      subclass="AluY"
    else if (name ~ /^AluS/)
      subclass="AluS"
    else if (name ~ /^AluJ/)
      subclass="AluJ"
    else if (name ~ /^FLAM/)
      subclass="FLAM"
    else if (name ~ /^FRAM/)
      subclass="FRAM"
    else
      subclass="Alu_other"
  }

  else if (class=="LTR_ERV") {
    if (name ~ /^HERVK/)
      subclass="HERVK"
    else if (name ~ /^HERV/)
      subclass="HERV_other"
    else if (name ~ /^LTR/)
      subclass="LTR"
    else if (name ~ /^THE1/)
      subclass="THE1"
    else if (name ~ /^MLT/)
      subclass="MLT"
    else if (name ~ /^MST/)
      subclass="MST"
    else if (name ~ /^ERV/)
      subclass="ERV_other"
    else
      subclass="LTR_misc"
  }

  else if (class=="DNA_transposon") {
    if (name ~ /^MER/)
      subclass="MER"
    else if (name ~ /^Charlie/)
      subclass="Charlie"
    else if (name ~ /^Tigger/)
      subclass="Tigger"
    else if (name ~ /^Zaphod/)
      subclass="Zaphod"
    else if (name ~ /^PABL/)
      subclass="PABL"
    else
      subclass="DNA_other"
  }

  else {
    subclass=class
  }

  # -------------------------
  # AGE GROUPING (UNCHANGED)
  # -------------------------
  if (name ~ /^SVA/ || name ~ /^HERVK/ || name ~ /^L1HS/ || name ~ /^L1PA1/)
    age="young"

  else if (name ~ /^AluY/ || name ~ /^L1PA/)
    age="young"

  else if (name ~ /^AluS/ || name ~ /^AluJ/ || name ~ /^L1PB/)
    age="intermediate"

  else if (name ~ /^MIR/ || name ~ /^L2/ || name ~ /^L3/ ||
           name ~ /^L1MA/ || name ~ /^L1MB/ || name ~ /^L1ME/ ||
           name ~ /^L1M/ || name ~ /^L1P/ ||
           name ~ /^MER/ || name ~ /^THE1/ || name ~ /^MST/ ||
           name ~ /^MLT/ || name ~ /^ERV/ || name ~ /^ERVL/ ||
           name ~ /^HERV/ || name ~ /^HUERS/ ||
           name ~ /^FLAM/ || name ~ /^FRAM/ || name ~ /^HAL1/)
    age="old"

  else if (name ~ /^\(/ || name ~ /rich/ || name ~ /^ALR/)
    age="old"

  else
    age="unknown"

  # -------------------------
  # SUMMARIES
  # -------------------------
  class_sum[class] += count
  age_sum[age] += count
  combo_sum[class "|" age] += count
  subclass_sum[class "|" subclass] += count
}

END {
  print "=== Class summary ==="
  for (c in class_sum)
    print c, class_sum[c]

  print "\n=== Age summary ==="
  for (a in age_sum)
    print a, age_sum[a]

  print "\n=== Class + Age breakdown ==="
  for (k in combo_sum)
    print k, combo_sum[k]

  print "\n=== Subclass breakdown ==="
  for (k in subclass_sum)
    print k, subclass_sum[k]
}
'
# === Class summary ===
# LTR_ERV 159
# Alu 225
# MIR 96
# L1 157
# L2 93
# L3 8
# Other 21
# Simple_repeat 52
# DNA_transposon 80
# 
# === Age summary ===
# young 86
# intermediate 199
# old 527
# unknown 79
# 
# === Class + Age breakdown ===
# Other|young 4
# L2|old 93
# L1|young 53
# Alu|young 27
# Other|old 7
# Alu|old 5
# LTR_ERV|old 106
# DNA_transposon|old 62
# L3|old 8
# L1|intermediate 6
# L1|old 98
# Alu|intermediate 193
# LTR_ERV|unknown 51
# Other|unknown 10
# DNA_transposon|unknown 18
# LTR_ERV|young 2
# MIR|old 96
# Simple_repeat|old 52
# 
# === Subclass breakdown ===
# LTR_ERV|MST 16
# LTR_ERV|LTR_misc 19
# DNA_transposon|PABL 3
# L2|L2 93
# Alu|FLAM 4
# Alu|AluJ 46
# Simple_repeat|Simple_repeat 52
# DNA_transposon|Zaphod 2
# Alu|AluY 27
# DNA_transposon|Tigger 8
# LTR_ERV|THE1 18
# L3|L3 8
# Alu|FRAM 1
# LTR_ERV|MLT 50
# Other|Other 21
# LTR_ERV|LTR 37
# LTR_ERV|ERV_other 8
# LTR_ERV|HERV_other 9
# MIR|MIR 96
# Alu|AluS 147
# DNA_transposon|MER 62
# L1|L1 157
# LTR_ERV|HERVK 2
# DNA_transposon|Charlie 5

############################### OVOL2 abundance analysis

#OVOL2 repetitive element abundance
cut -f7 ~/macs3_test/OVOL2_output_HOMER_RepeatMasker/peaks_in_repeats_OVOL2.bed | \
awk '
{
  # -------------------------
  # NORMALIZATION
  # -------------------------
  gsub(/_int$/, "")
  gsub(/-int$/, "")
  gsub(/_I-int$/, "")
  gsub(/\r/, "")

  print
}
' | sort | uniq -c | \
awk '
{
  count=$1
  name=$2

  # -------------------------
  # CLASSIFICATION
  # -------------------------
  if (name ~ /^Alu/ || name ~ /^FLAM/ || name ~ /^FRAM/)
    class="Alu"

  else if (name ~ /^MIR/)
    class="MIR"

  else if (name ~ /^L1/)
    class="L1"

  else if (name ~ /^L2/)
    class="L2"

  else if (name ~ /^L3/)
    class="L3"

  else if (name ~ /^HERV/ || name ~ /^HUERS/ || name ~ /^ERV/ || name ~ /^ERVL/ ||
           name ~ /^LTR/ || name ~ /^MLT/ || name ~ /^MST/ || name ~ /^THE1/ ||
           name ~ /^HAL1/ || name ~ /^LOR1/ ||
           name ~ /^MamGyp/ || name ~ /^MamRep/ || name ~ /^MamSINE/ ||
           name ~ /^MADE/ || name ~ /^PRIMAX/)
    class="LTR_ERV"

  else if (name ~ /^Charlie/ || name ~ /^Tigger/ || name ~ /^Zaphod/ ||
           name ~ /^MER/ || name ~ /^PABL/)
    class="DNA_transposon"

  else if (name ~ /^\(/ || name ~ /rich/)
    class="Simple_repeat"

  else
    class="Other"

  # -------------------------
  # AGE GROUPING
  # -------------------------
  if (name ~ /^SVA/ || name ~ /^HERVK/ || name ~ /^L1HS/ || name ~ /^L1PA1/)
    age="young"

  else if (name ~ /^AluY/ || name ~ /^L1PA/)
    age="young"

  else if (name ~ /^AluS/ || name ~ /^AluJ/ || name ~ /^L1PB/)
    age="intermediate"

  else if (name ~ /^MIR/ || name ~ /^L2/ || name ~ /^L3/ ||
           name ~ /^L1MA/ || name ~ /^L1MB/ || name ~ /^L1ME/ ||
           name ~ /^L1M/ || name ~ /^L1P/ ||
           name ~ /^MER/ || name ~ /^THE1/ || name ~ /^MST/ ||
           name ~ /^MLT/ || name ~ /^ERV/ || name ~ /^ERVL/ ||
           name ~ /^HERV/ || name ~ /^HUERS/ ||
           name ~ /^FLAM/ || name ~ /^FRAM/ || name ~ /^HAL1/)
    age="old"

  else if (name ~ /^\(/ || name ~ /rich/ || name ~ /^ALR/)
    age="old"

  else
    age="unknown"

  # -------------------------
  # SUMMARIES
  # -------------------------
  class_sum[class] += count
  age_sum[age] += count
  combo_sum[class "|" age] += count
}

END {
  print "=== Class summary ==="
  for (c in class_sum)
    print c, class_sum[c]

  print "\n=== Age summary ==="
  for (a in age_sum)
    print a, age_sum[a]

  print "\n=== Class + Age breakdown ==="
  for (k in combo_sum)
    print k, combo_sum[k]
}
'
# == Class summary ===
# LTR_ERV 983
# Alu 1845
# MIR 776
# L1 632
# L2 738
# L3 32
# Other 147
# Simple_repeat 1244
# DNA_transposon 415
# 
# === Age summary ===
# young 612
# intermediate 1378
# old 4119
# unknown 703
# 
# === Class + Age breakdown ===
# L2|old 738
# L1|young 149
# Alu|young 459
# Other|old 8
# Alu|old 43
# LTR_ERV|old 607
# DNA_transposon|old 351
# L1|unknown 123
# L3|old 32
# L1|intermediate 40
# L1|old 320
# Alu|intermediate 1338
# LTR_ERV|unknown 372
# Other|unknown 139
# DNA_transposon|unknown 64
# LTR_ERV|young 4
# MIR|old 776
# Alu|unknown 5
# Simple_repeat|old 1244

cut -f4 repeats_rhesus.bed | \
awk '
{
  gsub(/_int$/, "")
  gsub(/-int$/, "")
  gsub(/_I-int$/, "")
  gsub(/\r/, "")
  print
}
' | sort | uniq -c | \
awk '
{
  count=$1
  name=$2

  # CLASSIFICATION
  if (name ~ /^Alu/ || name ~ /^FLAM/ || name ~ /^FRAM/)
    class="Alu"
  else if (name ~ /^MIR/)
    class="MIR"
  else if (name ~ /^L1/)
    class="L1"
  else if (name ~ /^L2/)
    class="L2"
  else if (name ~ /^L3/)
    class="L3"
  else if (name ~ /^HERV/ || name ~ /^HUERS/ || name ~ /^ERV/ || name ~ /^ERVL/ ||
           name ~ /^LTR/ || name ~ /^MLT/ || name ~ /^MST/ || name ~ /^THE1/ ||
           name ~ /^HAL1/ || name ~ /^LOR1/ ||
           name ~ /^MamGyp/ || name ~ /^MamRep/ || name ~ /^MamSINE/ ||
           name ~ /^MADE/ || name ~ /^PRIMAX/)
    class="LTR_ERV"
  else if (name ~ /^Charlie/ || name ~ /^Tigger/ || name ~ /^Zaphod/ ||
           name ~ /^MER/ || name ~ /^PABL/)
    class="DNA_transposon"
  else if (name ~ /^\(/ || name ~ /rich/)
    class="Simple_repeat"
  else
    class="Other"

  class_sum[class] += count
}

END {
  for (c in class_sum)
    print c, class_sum[c]
}
'

# LTR_ERV 73403
# Alu 152037
# MIR 77388
# L1 96607
# L2 66822
# L3 7010
# Other 18362
# Simple_repeat 94696
# DNA_transposon 61062
#####enriched in Alu, ERVs and simple repeats: Adjust subcategories accordingly:

cut -f4 repeats_rhesus.bed | \
awk '
{
  # -------------------------
  # CLEANING
  # -------------------------
  gsub(/_int$/, "")
  gsub(/-int$/, "")
  gsub(/_I-int$/, "")
  gsub(/\r/, "")
  name=$0

  # -------------------------
  # MAIN CLASS
  # -------------------------
  if (name ~ /^Alu/ || name ~ /^FLAM/ || name ~ /^FRAM/)
    class="Alu"

  else if (name ~ /^MIR/)
    class="MIR"

  else if (name ~ /^L1/)
    class="L1"

  else if (name ~ /^L2/)
    class="L2"

  else if (name ~ /^L3/)
    class="L3"

  else if (name ~ /^HERV/ || name ~ /^HUERS/ || name ~ /^ERV/ || name ~ /^ERVL/ ||
           name ~ /^LTR/ || name ~ /^MLT/ || name ~ /^MST/ || name ~ /^THE1/ ||
           name ~ /^HAL1/ || name ~ /^LOR1/ ||
           name ~ /^MamGyp/ || name ~ /^MamRep/ || name ~ /^MamSINE/ ||
           name ~ /^MADE/ || name ~ /^PRIMAX/)
    class="LTR_ERV"

  else if (name ~ /^Charlie/ || name ~ /^Tigger/ || name ~ /^Zaphod/ ||
           name ~ /^MER/ || name ~ /^PABL/)
    class="DNA_transposon"

  else if (name ~ /^\(/ || name ~ /rich/)
    class="Simple_repeat"

  else
    class="Other"

  # -------------------------
  # SUBCLASS (KEY PART 🔥)
  # -------------------------

  subclass="Other"

  if (class=="Alu") {
    if (name ~ /^AluY/)
      subclass="AluY"
    else if (name ~ /^AluS/)
      subclass="AluS"
    else if (name ~ /^AluJ/)
      subclass="AluJ"
    else if (name ~ /^FLAM/)
      subclass="FLAM"
    else if (name ~ /^FRAM/)
      subclass="FRAM"
    else
      subclass="Alu_other"
  }

  else if (class=="LTR_ERV") {
    if (name ~ /^HERVK/)
      subclass="HERVK"
    else if (name ~ /^HERV/)
      subclass="HERV_other"
    else if (name ~ /^LTR/)
      subclass="LTR"
    else if (name ~ /^THE1/)
      subclass="THE1"
    else if (name ~ /^MLT/)
      subclass="MLT"
    else if (name ~ /^MST/)
      subclass="MST"
    else if (name ~ /^ERV/)
      subclass="ERV_other"
    else
      subclass="LTR_misc"
  }

  else if (class=="DNA_transposon") {
    if (name ~ /^MER/)
      subclass="MER"
    else if (name ~ /^Charlie/)
      subclass="Charlie"
    else if (name ~ /^Tigger/)
      subclass="Tigger"
    else if (name ~ /^Zaphod/)
      subclass="Zaphod"
    else if (name ~ /^PABL/)
      subclass="PABL"
    else
      subclass="DNA_other"
  }

  else if (class=="Simple_repeat") {

    # Extract motif inside parentheses
    motif=name

    gsub(/^\(/, "", motif)
    gsub(/\)n$/, "", motif)

    motif_len=length(motif)

    if (motif_len==1)
      subclass="mono_repeat"

    else if (motif_len==2)
      subclass="di_repeat"

    else if (motif_len==3)
      subclass="tri_repeat"

    else if (motif_len==4)
      subclass="tetra_repeat"

    else if (motif_len==5)
      subclass="penta_repeat"

    else if (motif_len==6)
      subclass="hexa_repeat"

    else
      subclass="long_repeat"
}

else {
    subclass=class
}

  # -------------------------
  # SUMMARIES
  # -------------------------
  class_sum[class]++
  subclass_sum[class "|" subclass]++
}

END {
  print "=== Class summary ==="
  for (c in class_sum)
    print c, class_sum[c]

  print "\n=== Subclass breakdown ==="
  for (k in subclass_sum)
    print k, subclass_sum[k]
}
'
# === Class summary ===
# LTR_ERV 73403
# Alu 152037
# MIR 77388
# L1 96607
# L2 66822
# Other 18362
# L3 7010
# Simple_repeat 94696
# DNA_transposon 61062
# 
# === Subclass breakdown ===
# LTR_ERV|MST 5199
# LTR_ERV|LTR_misc 10186
# DNA_transposon|PABL 158
# Alu|FLAM 5524
# Simple_repeat|di_repeat 23659
# L2|L2 66822
# Alu|AluJ 33574
# DNA_transposon|Zaphod 645
# Simple_repeat|mono_repeat 11959
# Alu|AluY 28363
# Simple_repeat|penta_repeat 7597
# Simple_repeat|tri_repeat 8397
# Simple_repeat|tetra_repeat 18656
# LTR_ERV|THE1 5957
# DNA_transposon|Tigger 9828
# L3|L3 7010
# Alu|FRAM 2465
# Alu|Alu_other 743
# Simple_repeat|long_repeat 8316
# Other|Other 18362
# LTR_ERV|MLT 28600
# LTR_ERV|ERV_other 2368
# LTR_ERV|LTR 18924
# Simple_repeat|hexa_repeat 16112
# LTR_ERV|HERV_other 1823
# MIR|MIR 77388
# Alu|AluS 81368
# LTR_ERV|HERVK 346
# L1|L1 96607
# DNA_transposon|MER 42152
# DNA_transposon|Charlie 8279


######now subcategories in my dataset:
cut -f7 ~/macs3_test/OVOL2_output_HOMER_RepeatMasker/peaks_in_repeats_OVOL2.bed | \
awk '
{
  # -------------------------
  # CLEANING
  # -------------------------
  gsub(/_int$/, "")
  gsub(/-int$/, "")
  gsub(/_I-int$/, "")
  gsub(/\r/, "")
  name=$0

  # -------------------------
  # MAIN CLASS
  # -------------------------
  if (name ~ /^Alu/ || name ~ /^FLAM/ || name ~ /^FRAM/)
    class="Alu"

  else if (name ~ /^MIR/)
    class="MIR"

  else if (name ~ /^L1/)
    class="L1"

  else if (name ~ /^L2/)
    class="L2"

  else if (name ~ /^L3/)
    class="L3"

  else if (name ~ /^HERV/ || name ~ /^HUERS/ || name ~ /^ERV/ || name ~ /^ERVL/ ||
           name ~ /^LTR/ || name ~ /^MLT/ || name ~ /^MST/ || name ~ /^THE1/ ||
           name ~ /^HAL1/ || name ~ /^LOR1/ ||
           name ~ /^MamGyp/ || name ~ /^MamRep/ || name ~ /^MamSINE/ ||
           name ~ /^MADE/ || name ~ /^PRIMAX/)
    class="LTR_ERV"

  else if (name ~ /^Charlie/ || name ~ /^Tigger/ || name ~ /^Zaphod/ ||
           name ~ /^MER/ || name ~ /^PABL/)
    class="DNA_transposon"

  else if (name ~ /^\(/ || name ~ /rich/)
    class="Simple_repeat"

  else
    class="Other"

  # -------------------------
  # SUBCLASS (KEY PART 🔥)
  # -------------------------

  subclass="Other"

  if (class=="Alu") {
    if (name ~ /^AluY/)
      subclass="AluY"
    else if (name ~ /^AluS/)
      subclass="AluS"
    else if (name ~ /^AluJ/)
      subclass="AluJ"
    else if (name ~ /^FLAM/)
      subclass="FLAM"
    else if (name ~ /^FRAM/)
      subclass="FRAM"
    else
      subclass="Alu_other"
  }

  else if (class=="LTR_ERV") {
    if (name ~ /^HERVK/)
      subclass="HERVK"
    else if (name ~ /^HERV/)
      subclass="HERV_other"
    else if (name ~ /^LTR/)
      subclass="LTR"
    else if (name ~ /^THE1/)
      subclass="THE1"
    else if (name ~ /^MLT/)
      subclass="MLT"
    else if (name ~ /^MST/)
      subclass="MST"
    else if (name ~ /^ERV/)
      subclass="ERV_other"
    else
      subclass="LTR_misc"
  }

  else if (class=="DNA_transposon") {
    if (name ~ /^MER/)
      subclass="MER"
    else if (name ~ /^Charlie/)
      subclass="Charlie"
    else if (name ~ /^Tigger/)
      subclass="Tigger"
    else if (name ~ /^Zaphod/)
      subclass="Zaphod"
    else if (name ~ /^PABL/)
      subclass="PABL"
    else
      subclass="DNA_other"
  }

  else if (class=="Simple_repeat") {

    # Extract motif inside parentheses
    motif=name

    gsub(/^\(/, "", motif)
    gsub(/\)n$/, "", motif)

    motif_len=length(motif)

    if (motif_len==1)
      subclass="mono_repeat"

    else if (motif_len==2)
      subclass="di_repeat"

    else if (motif_len==3)
      subclass="tri_repeat"

    else if (motif_len==4)
      subclass="tetra_repeat"

    else if (motif_len==5)
      subclass="penta_repeat"

    else if (motif_len==6)
      subclass="hexa_repeat"

    else
      subclass="long_repeat"
}

else {
    subclass=class
}

  # -------------------------
  # SUMMARIES
  # -------------------------
  class_sum[class]++
  subclass_sum[class "|" subclass]++
}

END {
  print "=== Class summary ==="
  for (c in class_sum)
    print c, class_sum[c]

  print "\n=== Subclass breakdown ==="
  for (k in subclass_sum)
    print k, subclass_sum[k]
}
'
# === Class summary ===
# LTR_ERV 983
# Alu 1845
# MIR 776
# L1 632
# L2 738
# L3 32
# Other 147
# Simple_repeat 1244
# DNA_transposon 415
# 
# === Subclass breakdown ===
# LTR_ERV|LTR_misc 46
# LTR_ERV|MST 70
# DNA_transposon|PABL 2
# Alu|FLAM 29
# Alu|AluJ 305
# Simple_repeat|di_repeat 149
# L2|L2 738
# DNA_transposon|Zaphod 4
# Simple_repeat|mono_repeat 6
# Simple_repeat|penta_repeat 135
# Alu|AluY 459
# L3|L3 32
# Simple_repeat|tetra_repeat 156
# DNA_transposon|Tigger 36
# Simple_repeat|tri_repeat 263
# LTR_ERV|THE1 106
# Alu|FRAM 14
# Alu|Alu_other 5
# Simple_repeat|long_repeat 193
# Other|Other 147
# LTR_ERV|MLT 308
# LTR_ERV|ERV_other 41
# LTR_ERV|LTR 343
# Simple_repeat|hexa_repeat 342
# LTR_ERV|HERV_other 65
# Alu|AluS 1033
# MIR|MIR 776
# LTR_ERV|HERVK 4
# DNA_transposon|MER 351
# L1|L1 632
# DNA_transposon|Charlie 22






