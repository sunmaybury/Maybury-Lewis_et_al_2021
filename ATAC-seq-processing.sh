#!/bin/bash

#Requested resources
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH -J atac
#SBATCH -o atac%J.out
#SBATCH-e atac%J.err

# Data processing started July 20, 2016
# Organized scripts into pipeline January 5, 2019 

# SML ATAC-Seq pipeline
# From fastq files to narrowPeak files (insertsite_blacklistcleared_narrowPeak)


# Run on Oscar

samples=(
	1_TAAGGCGA_L005
	2_CGTACTAG_L005
	3_AGGCAGAA_L005
	4_TCCTGAGC_L005
	5_GGACTCCT_L005
	6_TAGGCATG_L005
	7_TCCTGAGC_L003
	8_GGACTCCT_L003
	9_TAGGCATG_L003
	10_TAAGGCGA_L003
	11_CGTACTAG_L003
	12_AGGCAGAA_L003
	13_TAAGGCGA_L004
	14_CGTACTAG_L004
	15_AGGCAGAA_L004
	16_TCCTGAGC_L004
	17_GGACTCCT_L004
	18_TAGGCATG_L005
	19_TAAGGCGA_L005
	20_CGTACTAG_L005
	21_AGGCAGAA_L005
	22_TCCTGAGC_L005
	23_GGACTCCT_L008
	24_TAGGCATG_L008
	25_TAAGGCGA_L008
	26_CGTACTAG_L008
	27_AGGCAGAA_L008)


############################
# Trimming
############################

# Old module
# module load trimgalore/0.4.0

# New modules
module load trimgalore/0.5.0
module load fastqc/0.11.5
module load cutadapt/1.14

for read1 in *_R1_001.fastq.gz
do
	read2=$(echo $read1| sed 's/R1_001.fastq.gz/R2_001.fastq.gz/')
	trim_galore $read1 $read2 --paired --fastqc
done

######################################
# Mapping with Bowtie2
######################################

# Old module
# module load bowtie2/2.2.5

# New module
module load bowtie2/2.3.0

for i in "${samples[@]}"; do
	bowtie2 -p8 -q --no-discordant -x /gpfs/data/shared/bowtie/mm10 -1 "$i"_R1_001_val_1.fq.gz -2 "$i"_R2_001_val_2.fq.gz -S "$i".sam
done

############################################
# Filter and remove duplicates
############################################

for i in "${samples[@]}"; do
echo Filter chrM $i
	sed '/chrM/d' < "$i".sam > "$i"_filtered.sam 		# filter mitochondrial chromosome
done

module load samtools/1.3.1

for i in "${samples[@]}"; do
echo Convert to bam $i
	samtools view -b -S "$i"_filtered.sam > "$i"_filtered.bam 		# convert to bam and sort
	samtools sort "$i"_filtered.bam -o "$i"_sorted.bam
echo Post alignment filter $i
	samtools view -F 1804 -b "$i"_sorted.bam > "$i"_filtered_FILT.bam 	# flag 1804
done

# Old modules
# module load java/7u5
# module load picard-tools/1.88  

# New modules
module load java/8u111
module load picard-tools/2.9.2


######### mark duplicates with picard


# MARKDUP="/gpfs/runtime/opt/picard-tools/1.88/picard-tools-1.88/MarkDuplicates.jar"
MARKDUP="/gpfs/runtime/opt/picard-tools/2.9.2/picard.jar"

for i in ${samples[@]}; do
	echo Mark duplicates chrM filtered $i 			
		java -Xmx4g -jar ${MARKDUP} MarkDuplicates \
			I=${i}_filtered_FILT.bam \
			O=${i}_filtered_markdup.bam \
			M=${i}_filtered_log.txt \
			REMOVE_DUPLICATES=false \
			VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp
done


######### filter marked duplicates with samtools

for i in "${samples[@]}"; do
echo Remove duplicates chrM filtered $i
	mv "$i"_filtered_markdup.bam "$i"_filtered_FILT.bam 
	samtools view -F 1804 -b "$i"_filtered_FILT.bam > "$i"_filtered_FINAL.bam
echo Index final BAM file $i
	samtools index "$i"_filtered_FINAL.bam "$i"_filtered_FINAL_INDEX.bai
	samtools flagstat "$i"_filtered_FINAL.bam > "$i"_filtered_FINAL_flagstat
done


### Save the final bam files locally
### Following steps are done locally, but you need to install: samtools, bedtools, MACS (optional)
### These tools are good to have on your computer if possible. If not, you can just run on Oscar (all three modules are available on Oscar)


############################# Start with filtered and dupsremoved bam files ##############################################

## Run locally

# FILE_PATH="/Volumes/External/Sequencing/TruSeq_ATAC-Seq"

###########################
# Get insert size metrics
###########################

# # Run locally
# cd $FILE_PATH/20160814_dupsremoved_finalbam 		# Start with filtered and dupsremoved bam files
# for i in "${samples[@]}"; do
# 	echo Get insert sizes $i
# 		java -Xms1g -Xmx16g -jar /Users/syk3/Softwares/picard-2.7.1/picard-2.7.1.jar CollectInsertSizeMetrics METRIC_ACCUMULATION_LEVEL=ALL_READS OUTPUT="${i}_insertSizes.txt" HISTOGRAM_FILE="{i}_insertSizes.pdf" INPUT="${i}_filtered_FINAL.bam"
# done

########### Or run on Oscar
module load java/8u111
module load picard-tools/2.9.2
module load R/3.5.2

PICARD="/gpfs/runtime/opt/picard-tools/2.9.2/picard.jar"

for i in ${samples[@]}; do
	echo Get insert sizes $i
		java -Xms1g -Xms16g -jar ${PICARD} CollectInsertSizeMetrics \
			I=${i}_filtered_FINAL.bam \
			O=${i}_insertSizes.txt \
			H=${i}_insertSizes.pdf \
			METRIC_ACCUMULATION_LEVEL=ALL_READS
done



#############################################
# Get insert sites
# Also adjust for Tn5 binding footprint
#############################################

# Run locally

# Convert bam to bed
# Then adjust 4 or 5bp for the Tn5 binding site, and trim to the 5' end of the read

cd $FILE_PATH/20160814_dupsremoved_finalbam
for i in "${samples[@]}"; do
	echo get insert sites $i
		bamToBed -i "${i}_filtered_FINAL.bam" > temp.bed
		samtools view "${i}_filtered_FINAL.bam" | perl -lane 'print $F[0]."\t".$F[8];' > insert.temp

		paste temp.bed insert.temp | \
		perl -lane '$name=(split/\//,$F[3])[0]; if ($name =~ /$F[6]/){ print $F[0]."\t".$F[1]."\t".$F[2]."\t".$F[3]."\t".$F[7]."\t".$F[5]; }else{ print STDERR "NAMES DONT MATCH for $name $F[6]!!!!";}' | \
		awk -F $'\t' 'BEGIN {OFS = FS}{if ($6 == "+") { $2 = $2 + 4 } else if ($6 == "-") {$3 = $3 - 5} print $0}' | \
		awk 'BEGIN{OFS="\t"}{if($6 == "-") $2=$3-1; print $1, $2, $2+1, $4, $5, $6}' | \
		sort -k1,1 -k2,2n - >> "${i}_insertSites.bed"

		rm temp.bed insert.temp

	echo convert back to bam $i
		bedToBam -i "${i}_insertSites.bed" -g mm10.chrom.sizes > temp.bam
		samtools sort temp.bam -o "${i}_sort.bam"
		rm temp.bam
		samtools index "${i}_sort.bam" "${i}_sort.bam.bai"
done

# # Clean up and move files to new folder
# cd $FILE_PATH
# mkdir 20161118_insertsite_bed
# mkdir 20161118_insertsite_bam
# cd $FILE_PATH/20160814_dupsremoved_finalbam
# mv *_insert_site.bed $FILE_PATH/20161118_insertsite_bed
# mv *_sort.bam $FILE_PATH/20161118_insertsite_bam
# mv *_sort.bam.bai $FILE_PATH/20161118_insertsite_bam


# Get insertsites bed files ready for calling peaks (with MACS - extending windows from the insertsites)
# insertsite.bed files were zipped 


cd $FILE_PATH/20161118_insertsite_bed 
for i in "${samples[@]}"; do
	echo get ready to call peaks $i
		gunzip -c "${i}_insertSites.bed.gz" | slopBed -i - -g mm10.chrom.sizes -l 75 -r -75 -s > "${i}_insertSites_75bpShiftedForMacs.bed"
done 

# Clean up
mkdir 75bpShiftedForMacs
mv *_insertSites_75bpShiftedForMacs.bed $FILE_PATH/20161118_insertsite_bed_75bpShiftedForMacs


################################################################
# Call peaks
################################################################

# module load macs/2.1.1
# module load bedtools/2.25.0

# Ran locally

# Call peaks
for i in "${samples[@]}"; do
	echo Decompress file
		gzip -d $i.bed.gz 
	echo Call peaks
		macs2 callpeak -t ${i}_insertSites_75bpShiftedForMacs.bed -f BED -n "${i}" -g mm -q 0.05 --nomodel --extsize 150 -B --keep-dup all --call-summits
done

# Split subpeaks
for i in "${samples[@]}"; do
	echo Split peaks
		perl splitMACS2SubPeaks.pl ${i}_peaks.narrowPeak > ${i}_splitSubPeaks.narrowPeak
done

# Clear blacklist peaks mm10 using Bedtools
BLACKLIST_PEAKS="mm10.blacklist.bed"
BLACKLIST_CLEARED_PEAKS="splitSubPeaks.blacklistCleared.narrowPeak"
BLACKLIST_LOST_PEAKS="splitSubPeaks.blacklistLost.narrowPeak"

for i in "${samples[@]}"; do
	echo intersect blacklist peaks
		intersectBed -v -a "${i}_splitSubPeaks.narrowPeak" -b "${BLACKLIST_PEAKS}" > "${i}_${BLACKLIST_CLEARED_PEAKS}"
		intersectBed -u -a "${i}_splitSubPeaks.narrowPeak" -b "${BLACKLIST_PEAKS}" > "${i}_${BLACKLIST_LOST_PEAKS}"
done





################################################################
# Make UCSC bedgraph files
################################################################

# cd $FILE_PATH/20190109_UCSCfile

BAM_PATH="/Volumes/research/BM_WebbLab/SunML/AB_ATAC-seq/20190108_insertsite_bam"

# Make tag directory
for i in "${samples[@]}"; do
	echo making tag directory $i
		makeTagDirectory "$i"-ATAC-seq/ $BAM_PATH/"$i"_sort.bam
done

# Remove reads outside UCSC chromosome lengths
for i in "${samples[@]}"; do
	echo remove out of bound reads $i
		removeOutofBoundsReads.pl "$i"-ATAC-Seq/ mm10
done

# Make genome browser files
for i in "${samples[@]}"; do
	makeUCSCfile "$i"-ATAC-Seq/ -o auto -fsize 1e8
done





