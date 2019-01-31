#################################################
##### Whole Exome Sequencing data treatment #####

#!/bin/bash

#Directory with demultiplexed fastq.gz files
input="../../datashare/Spruce/exome_capture/fastq/cleaned_reads/"

#Directory for outputs
output="./output"

############ STEP 1 ##############
### Sorting and filtering data ###

# numer of reads in individuals
rm ${output}/NBreads_1P.txt
rm ${output}/NBreads_2P.txt
rm ${output}/NBreads_names.txt
for forward in ${input}*_fastq_1P.gz
do
	reverse=${forward/_fastq_1P.gz/_fastq_2P.gz}
	f=${forward/_fastq_1P.gz/}
	name=`basename ${f}`
	echo ${name} >> ${output}/NBreads_names.txt
	echo $(zcat ${forward} | wc -l)/4|bc >> ${output}/NBreads_1P.txt
	echo $(zcat ${reverse} | wc -l)/4|bc >> ${output}/NBreads_2P.txt
done

R --vanilla <<EOF

	reads_R1 = read.table("${output}/NBreads_1P.txt",header = FALSE)
	reads_R2 = read.table("${output}/NBreads_2P.txt",header = FALSE)
	ind_names = read.table("${output}/NBreads_names.txt", header = FALSE)
	bwa_res = cbind (ind_names, reads_R1, reads_R2)
	colnames(bwa_res) = c("Individuals", "Total_reads", "Total_mapped", "Total_mapq10", "Total_mapq30")
	write.table(bwa_res, "${output}/bwa_results.txt")
EOF

################## STEP 2 ###################
### Mapping the reads with bwa and stampy ###

#Reference for mapping
ref="./WES_RedSpurce_reference.fna"

#number of CPU used
t=20

### First round of mapping with bwa ###

# Indexing the genome
#bwa index ${ref}

# Aligning individual sequences to the reference
# Running in parallel to speed up the mapping

for forward in ${input}O*_fastq_1P.gz
do
	f=${forward/_fastq_1P.gz/}
	name=`basename ${f}`
	echo "@ Aligning $name..."
	reverse=${forward/_fastq_1P.gz/_fastq_2P.gz}
    bwa mem -t ${t} -M -a ${ref} ${forward} ${reverse} > ${output}/aligned/BWA/${name}.sam
done

ls ${output}/aligned/BWA/*.sam | parallel -j${t} "samtools view -bS {} | samtools sort {} -n -T {.} -o {.}.bam" #; samtools index {.}.bam"

# Stats on bwa alignments
rm ${output}/aligned/BWA/res.aln.reads.out
rm ${output}/aligned/BWA/names.txt
for INDIV in ${output}/aligned/BWA/*.bam
do
	f=${INDIV/.bam/}
	name=`basename ${f}`
	echo ${name} >> ${output}/aligned/BWA/names.txt
	samtools flagstat ${INDIV} | awk 'NR>=6&&NR<=13 {print $1}' | column -x
done >> ${output}/aligned/BWA/res.aln.reads.out

R --vanilla <<EOF

	align_res = read.table("${output}/aligned/BWA/res.aln.reads.out",header = FALSE)
	ind_names = read.table("${output}/aligned/BWA/names.txt", header = FALSE)
	bwa_res = cbind (ind_names, align_res)
	colnames(bwa_res) = c("Individuals", "Total_reads", "R1_reads", "R2_reads", "Properly_paired", "ReadMapped_MateUnmapped", "ReadUnmapped_MateMapped", "MateDifferentContig", "MateDifferentContig_Q5")
	bwa_res$Percent_properlypaired = as.numeric(bwa_res$Properly_paired) / as.numeric(bwa_res$Total_reads)
	write.table(bwa_res, "${output}/bwa_results.txt")
EOF

# Reads mapping quality scores after bwa alignment
rm ${output}/aligned/BWA/reads_mapping_Qscores.txt
for file in ${output}/aligned/BWA/*.bam
do
	samtools view ${file} | awk '$5>0{c1++}; $5>29{c29++}; $5>=0{c0++}; $5>19{c19++}; $5>9{c9++} END {print c0 " " c1 " " c9 " " c19 " " c29}'
done >> ${output}/aligned/BWA/reads_mapping_Qscores.txt

# Nucleotide coverage on bwa .bam files
rm ${output}/aligned/BWA/mean_coverage.txt
for file in ${output}/aligned/BWA/*.bam
do
	samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'
done >> ${output}/aligned/BWA/mean_coverage.txt

R --vanilla <<EOF

	reads_Q = read.table("${output}/aligned/BWA/reads_mapping_Qscores.txt",header = FALSE)
	mean_cov = read.table("${output}/aligned/BWA/mean_coverage.txt",header = FALSE)
	ind_names = read.table("${output}/aligned/BWA/names.txt", header = FALSE)
	qual_res = cbind (ind_names, reads_Q, mean_cov)
	colnames(qual_res) = c("Individuals", "Total_reads", "Total_mapped", "Total_mapq10", "Total_mapq20", "Total_mapq30", "Mean_cov")
	write.table(qual_res, "${output}/bwa_quality_results.txt")
EOF

# Deleting all the sam files to save space
rm ${output}/aligned/BWA/*.sam

### Remove PCR duplicates
ls ${output}/aligned/BWA/*.bam | parallel -j${t} "sambamba-0.6.8 markdup -r -t 2 {} {.}_rmdup.bam"

# Stats after PCR duplicates removal
rm ${output}/aligned/BWA/rmdup.res.aln.reads.out
rm ${output}/aligned/BWA/rmdup.names.txt
for INDIV in ${output}/aligned/BWA/*_rmdup.bam
do
	f=${INDIV/.bam/}
	name=`basename ${f}`
	echo ${name} >> ${output}/aligned/BWA/rmdup.names.txt
	samtools flagstat ${INDIV} | awk 'NR>=6&&NR<=13 {print $1}' | column -x
done >> ${output}/aligned/BWA/rmdup.res.aln.reads.out

R --vanilla <<EOF

	align_res = read.table("${output}/aligned/BWA/rmdup.res.aln.reads.out",header = FALSE)
	ind_names = read.table("${output}/aligned/BWA/rmdup.names.txt", header = FALSE)
	bwa_res = cbind (ind_names, align_res)
	colnames(bwa_res) = c("Individuals", "Total_reads", "R1_reads", "R2_reads", "Properly_paired", "ReadMapped_MateUnmapped", "ReadUnmapped_MateMapped", "MateDifferentContig", "MateDifferentContig_Q5")
	Percent_properlypaired = as.numeric(bwa_res[,5]) / as.numeric(bwa_res[,2])
	bwa_res = data.frame(bwa_res, Percent_properlypaired=Percent_properlypaired)
	write.table(bwa_res, "${output}/bwa_rmdup_results.txt")
EOF

# Reads mapping quality scores
rm ${output}/aligned/BWA/rmdup_reads_mapping_Qscores.txt
for file in ${output}/aligned/BWA/*_rmdup.bam
do
	samtools view ${file} | awk '$5>0{c1++}; $5>29{c29++}; $5>=0{c0++}; $5>19{c19++}; $5>9{c9++} END {print c0 " " c1 " " c9 " " c19 " " c29}'
done >> ${output}/aligned/BWA/rmdup_reads_mapping_Qscores.txt

# Nucleotide coverage
rm ${output}/aligned/BWA/rmdup_mean_coverage.txt
for file in ${output}/aligned/BWA/*.bam
do
	samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'
done >> ${output}/aligned/BWA/rmdup_mean_coverage.txt

R --vanilla <<EOF

	reads_Q = read.table("${output}/aligned/BWA/rmdup_reads_mapping_Qscores.txt",header = FALSE)
	mean_cov = read.table("${output}/aligned/BWA/rmdup_mean_coverage.txt",header = FALSE)
	ind_names = read.table("${output}/aligned/BWA/rmdup.names.txt", header = FALSE)
	qual_res = cbind (ind_names, reads_Q, mean_cov)
	colnames(qual_res) = c("Individuals", "Total_reads", "Total_mapped", "Total_mapq10", "Total_mapq20", "Total_mapq30", "Mean_cov")
	write.table(qual_res, "${output}/bwa_rmdup_quality_results.txt")
EOF


#####################################
### Mapping the reads with Stampy ###

### Re-do the mapping with Stampy from the bam files and keeping the well-mapped reads ###
mkdir ${output}/aligned/stampy
stampy=~/TOOLS/stampy-1.0.32/stampy.py

#Build a genome (.stidx) file:
${stampy} -t ${t} -G WES_RedSpurce_reference ${ref}

#Build a hash (.sthash) file:
${stampy} -t ${t} -g WES_RedSpurce_reference -H WES_RedSpurce_reference

for file in ${output}/aligned/BWA/*.bam
do
	f=${file/.bam/}
	name=`basename ${f}`
	echo "@ Aligning $name..."
	${stampy} -t${t} -g ${ref/.fna/} -h ${ref/.fna/} --bamkeepgoodreads -M ${output}/aligned/${name}.bam > ${output}/aligned/stampy/${name}.stampy.sam
done

ls ${output}/aligned/stampy/*.stampy.sam | parallel -j${t} "samtools view -bS {} | samtools sort - {.}" #; samtools index {.}.bam"

# Deleting all the stampy.sam files to save space
rm ${output}/aligned/stampy/*.stampy.sam

# Stats on stampy alignments
rm ${output}/aligned/stampy/res.aln.reads.stampy.out
rm ${output}/aligned/stampy/stampy_names.txt
for INDIV in ${output}/aligned/stampy/*.stampy.bam
do
	f=${INDIV/.stampy.bam/}
	name=`basename ${f}`
	echo ${name} >> ${output}/aligned/stampy/stampy_names.txt
	samtools flagstat ${INDIV} | awk 'NR>=6&&NR<=13 {print $1}' | column -x
done >> ${output}/aligned/stampy/res.aln.reads.stampy.out

R --vanilla <<EOF

	align_res = read.table("${output}/aligned/stampy/res.aln.reads.stampy.out",header = FALSE)
	ind_names = read.table("${output}/aligned/stampy/stampy_names.txt", header = FALSE)
	stampy_res = cbind(ind_names, align_res)
	colnames(stampy_res) = c("Individuals", "Total_reads", "R1_reads", "R2_reads", "Properly_paired", "ReadMapped_MateUnmapped", "ReadUnmapped_MateMapped", "MateDifferentContig", "MateDifferentContig_Q5")
	Percent_properlypaired = as.numeric(stampy_res[,5]) / as.numeric(stampy_res[,2])
	stampy_res = data.frame(stampy_res, Percent_properlypaired=Percent_properlypaired)
	write.table(stampy_res, "${output}/stampy_results.txt")
EOF

# Reads mapping quality scores after stampy alignment
rm ${output}/aligned/stampy/reads_mapping_Qscores.txt
for file in ${output}/aligned/stampy/*.stampy.bam
do
	samtools view ${file} | awk '$5>0{c1++}; $5>29{c29++}; $5>=0{c0++}; $5>19{c19++}; $5>9{c9++} END {print c0 " " c1 " " c9 " " c19 " " c29}'
done >> ${output}/aligned/stampy/reads_mapping_Qscores.txt

# Nucleotide coverage on stampy .bam files
rm ${output}/aligned/stampy/mean_coverage.txt
for file in ${output}/aligned/stampy/*.stampy.bam
do
	samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'
done >> ${output}/aligned/stampy/mean_coverage.txt

R --vanilla <<EOF

	reads_Q = read.table("${output}/aligned/stampy/reads_mapping_Qscores.txt",header = FALSE)
	mean_cov = read.table("${output}/aligned/stampy/mean_coverage.txt",header = FALSE)
	ind_names = read.table("${output}/aligned/stampy/stampy_names.txt", header = FALSE)
	qual_res = cbind (ind_names, reads_Q, mean_cov)
	colnames(qual_res) = c("Individuals", "Total_reads", "Total_mapped", "Total_mapq10", "Total_mapq20", "Total_mapq30", "Mean_cov")
	write.table(qual_res, "${output}/stampy_quality_results.txt")
EOF


#########################################
### Mapping the reads with NextGenMap ###

mkdir ${output}/aligned/NGM
for forward in ${input}*_fastq_1P.gz
do
	f=${forward/_fastq_1P.gz/}
	name=`basename ${f}`
	echo "@ Aligning $name..."
	reverse=${forward/_fastq_1P.gz/_fastq_2P.gz}
    ngm -r ${ref} -1 ${forward} -2 ${reverse} -o ${output}/aligned/NGM/${name}.ngm.sam -t ${t}
done

ls ${output}/aligned/NGM/*.ngm.sam | parallel -j${t} "samtools view -bS {} | samtools sort - {.}" #; samtools index {.}.bam"

# Deleting all the ngm.sam files to save space
rm ${output}/aligned/NGM/*.ngm.sam

# Stats on NGM alignments
rm ${output}/aligned/NGM/res.aln.reads.ngm.out
rm ${output}/aligned/NGM/ngm_names.txt
for INDIV in ${output}/aligned/NGM/*.ngm.bam
do
	f=${INDIV/.ngm.bam/}
	name=`basename ${f}`
	echo ${name} >> ${output}/aligned/NGM/ngm_names.txt
	samtools flagstat ${INDIV} | awk 'NR>=6&&NR<=13 {print $1}' | column -x
done >> ${output}/aligned/NGM/res.aln.reads.ngm.out

R --vanilla <<EOF

	align_res = read.table("${output}/aligned/NGM/res.aln.reads.ngm.out",header = FALSE)
	ind_names = read.table("${output}/aligned/NGM/ngm_names.txt", header = FALSE)
	ngm_res = cbind(ind_names, align_res)
	colnames(ngm_res) = c("Individuals", "Total_reads", "R1_reads", "R2_reads", "Properly_paired", "ReadMapped_MateUnmapped", "ReadUnmapped_MateMapped", "MateDifferentContig", "MateDifferentContig_Q5")
	Percent_properlypaired = as.numeric(ngm_res[,5]) / as.numeric(ngm_res[,2])
	stampy_res = data.frame(ngm_res, Percent_properlypaired=Percent_properlypaired)
	write.table(ngm_res, "${output}/ngm_results.txt")
EOF

# Reads mapping quality scores after NGM alignment
rm ${output}/aligned/NGM/reads_mapping_Qscores.txt
for file in ${output}/aligned/NGM/*.ngm.bam
do
	samtools view ${file} | awk '$5>0{c1++}; $5>29{c29++}; $5>=0{c0++}; $5>19{c19++}; $5>9{c9++} END {print c0 " " c1 " " c9 " " c19 " " c29}'
done >> ${output}/aligned/NGM/reads_mapping_Qscores.txt

# Nucleotide coverage on NGM .bam files
rm ${output}/aligned/NGM/mean_coverage.txt
for file in ${output}/aligned/NGM/*.ngm.bam
do
	samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'
done >> ${output}/aligned/NGM/mean_coverage.txt

R --vanilla <<EOF

	reads_Q = read.table("${output}/aligned/NGM/reads_mapping_Qscores.txt",header = FALSE)
	mean_cov = read.table("${output}/aligned/NGM/mean_coverage.txt",header = FALSE)
	ind_names = read.table("${output}/aligned/NGM/ngm_names.txt", header = FALSE)
	qual_res = cbind (ind_names, reads_Q, mean_cov)
	colnames(qual_res) = c("Individuals", "Total_reads", "Total_mapped", "Total_mapq10", "Total_mapq20", "Total_mapq30", "Mean_cov")
	write.table(qual_res, "${output}/ngm_quality_results.txt")
EOF
