###########################
#### MAPPING THE READS ####
###########################

#!/bin/bash

#Directory with demultiplexed fastq.gz files
input=$1

#Directory for outputs
output=$2

#Reference for mapping (fasta file)
ref=$3

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


########### STEP 2 ###############
### Mapping the reads with bwa ###

#number of CPU used
t=10

# Indexing the genome
#bwa index ${ref}

# Aligning individual sequences to the reference
for forward in ${input}*_fastq_1P.gz
do
	f=${forward/_fastq_1P.gz/}
	name=`basename ${f}`
	echo "@ Aligning $name..."
	reverse=${forward/_fastq_1P.gz/_fastq_2P.gz}
    bwa mem -t ${t} -M -a ${ref} ${forward} ${reverse} > ${output}/BWA/${name}.sam
done

### Sorting SAM files and converting to BAM files
for f in ${output}/BWA/*.sam
do
	out=${f/.sam/}
	samtools view -bS ${f} -o ${out}.bam
	samtools sort ${out}.bam -o ${out}.sorted.bam
	rm ${out}.bam
done

# Stats on bwa alignments
rm ${output}/BWA/res.aln.reads.out
rm ${output}/BWA/names.txt
for INDIV in ${output}/BWA/*.sorted.bam
do
	f=${INDIV/.bam/}
	name=`basename ${f}`
	echo ${name} >> ${output}/BWA/names.txt
	samtools flagstat ${INDIV} | awk 'NR>=6&&NR<=13 {print $1}' | column -x
done >> ${output}/BWA/res.aln.reads.out

R --vanilla <<EOF

	align_res = read.table("${output}/BWA/res.aln.reads.out",header = FALSE)
	ind_names = read.table("${output}/BWA/names.txt", header = FALSE)
	bwa_res = cbind (ind_names, align_res)
	colnames(bwa_res) = c("Individuals", "Total_reads", "R1_reads", "R2_reads", "Properly_paired", "ReadMapped_MateUnmapped", "ReadUnmapped_MateMapped", "MateDifferentContig", "MateDifferentContig_Q5")
	Percent_properlypaired = as.numeric(bwa_res[,5]) / as.numeric(bwa_res[,2])
	bwa_res = data.frame(bwa_res, Percent_properlypaired=Percent_properlypaired)
	write.table(bwa_res, "${output}/BWA/bwa_results.txt")
EOF

# Reads mapping quality scores after bwa alignment
rm ${output}/BWA/reads_mapping_Qscores.txt
for file in ${output}/BWA/*.sorted.bam
do
	samtools view ${file} | awk '$5>0{c1++}; $5>29{c29++}; $5>=0{c0++}; $5>19{c19++}; $5>9{c9++} END {print c0 " " c1 " " c9 " " c19 " " c29}'
done >> ${output}/BWA/reads_mapping_Qscores.txt

# Nucleotide coverage on bwa .bam files
rm ${output}/BWA/mean_coverage.txt
for file in ${output}/BWA/*.sorted.bam
do
	samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'
done >> ${output}/BWA/mean_coverage.txt

R --vanilla <<EOF

	reads_Q = read.table("${output}/BWA/reads_mapping_Qscores.txt",header = FALSE)
	mean_cov = read.table("${output}/BWA/mean_coverage.txt",header = FALSE)
	ind_names = read.table("${output}/BWA/names.txt", header = FALSE)
	qual_res = cbind (ind_names, reads_Q, mean_cov)
	colnames(qual_res) = c("Individuals", "Total_reads", "Total_mapped", "Total_mapq10", "Total_mapq20", "Total_mapq30", "Mean_cov")
	write.table(qual_res, "${output}/BWA/bwa_quality_results.txt")
EOF


### Removing PCR duplicates
for file in ${output}/BWA/*.sorted.bam
do
	sambamba-0.6.8 markdup -r -t 2 ${file} ${file}.rmdup.bam
	samtools sort ${file}.rmdup.bam -o ${file}.rmdup.bam.sorted.bam
done

for file in ${output}/BWA/*.sorted.bam.rmdup.bam.sorted.bam
do
	f=${file/.sorted.bam.rmdup.bam.sorted.bam/}
	mv ${file} ${f}.final.bam
done

for file in ${output}/BWA/*.final.bam
do
	samtools index ${file}
done

# Deleting all the sam files and other temporary files to save space
rm ${output}/BWA/*.sam
rm ${output}/BWA/*.sorted.bam
rm ${output}/BWA/*.rmdup.bam
rm ${output}/BWA/*.rmdup.bam.bai

# Stats after PCR duplicates removal
rm ${output}/BWA/rmdup.res.aln.reads.out
rm ${output}/BWA/rmdup.names.txt
for INDIV in ${output}/BWA/*.final.bam
do
	f=${INDIV/final.bam/}
	name=`basename ${f}`
	echo ${name} >> ${output}/BWA/rmdup.names.txt
	samtools flagstat ${INDIV} | awk 'NR>=6&&NR<=13 {print $1}' | column -x
done >> ${output}/BWA/rmdup.res.aln.reads.out

R --vanilla <<EOF

	align_res = read.table("${output}/BWA/rmdup.res.aln.reads.out",header = FALSE)
	ind_names = read.table("${output}/BWA/rmdup.names.txt", header = FALSE)
	bwa_res = cbind (ind_names, align_res)
	colnames(bwa_res) = c("Individuals", "Total_reads", "R1_reads", "R2_reads", "Properly_paired", "ReadMapped_MateUnmapped", "ReadUnmapped_MateMapped", "MateDifferentContig", "MateDifferentContig_Q5")
	Percent_properlypaired = as.numeric(bwa_res[,5]) / as.numeric(bwa_res[,2])
	bwa_res = data.frame(bwa_res, Percent_properlypaired=Percent_properlypaired)
	write.table(bwa_res, "${output}/BWA/bwa_rmdup_results.txt")
EOF

# Reads mapping quality scores
rm ${output}/BWA/rmdup_reads_mapping_Qscores.txt
for file in ${output}/BWA/*.final.bam
do
	samtools view ${file} | awk '$5>0{c1++}; $5>29{c29++}; $5>=0{c0++}; $5>19{c19++}; $5>9{c9++} END {print c0 " " c1 " " c9 " " c19 " " c29}'
done >> ${output}/BWA/rmdup_reads_mapping_Qscores.txt

# Nucleotide coverage
rm ${output}/BWA/rmdup_mean_coverage.txt
for file in ${output}/BWA/*.final.bam
do
	samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'
done >> ${output}/BWA/rmdup_mean_coverage.txt

R --vanilla <<EOF

	reads_Q = read.table("${output}/BWA/rmdup_reads_mapping_Qscores.txt",header = FALSE)
	mean_cov = read.table("${output}/BWA/rmdup_mean_coverage.txt",header = FALSE)
	ind_names = read.table("${output}/BWA/rmdup.names.txt", header = FALSE)
	qual_res = cbind (ind_names, reads_Q, mean_cov)
	colnames(qual_res) = c("Individuals", "Total_reads", "Total_mapped", "Total_mapq10", "Total_mapq20", "Total_mapq30", "Mean_cov")
	write.table(qual_res, "${output}/BWA/bwa_rmdup_quality_results.txt")
EOF
