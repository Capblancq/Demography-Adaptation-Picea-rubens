##########################
####### ANNOTATIONS ######
##########################

#!/bin/bash

###########################################
### Annotation of the polymorphic sites ###

### Finding the match between White and Norway spruce references ###
# Scaffold and positions of all the variants 
zcat ~/mydata/datashare/Spruce/exome_capture/WES_mapping/ANGSD/ref_Pglauca/Full_Sampling_intersect_poly.vcf.gz | grep -v "#" | cut -f -2 - > temp.bed
zcat ~/mydata/datashare/Spruce/exome_capture/WES_mapping/ANGSD/ref_Pglauca/Full_Sampling_intersect_poly.vcf.gz | grep -v "#" | cut -f 2 - | paste temp.bed - > chrompos_vcf_WS.bed
rm temp.bed

# Creating the new scaffolds with 300 bp on each side of every variants
cut -f 1,2 ~/mydata/datashare/Spruce/exome_capture/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna.fai > chrom.sizes
bedtools slop -i chrompos_vcf_WS.bed -g chrom.sizes -b 300 > near_variants.bed
awk '$1 != "NA"' FS=' ' near_variants.bed > near_variants_filtered.bed
bedtools getfasta -fi ~/mydata/datashare/Spruce/exome_capture/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -bed ./near_variants.bed > near_variants.fasta

# blastn white spruce scaffolds against norway spruce reference
/popgen/blast/ncbi-blast-2.7.1+/bin/blastn -query near_variants.fasta -db ~/mydata/datashare/Spruce/exome_capture/WES_mapping/ReferenceGenomes/Pabies1.0-genome.fa -task blastn -num_threads 15 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" -out blast_results_variants.out -max_target_seqs 1


### Changing the name of the scaffold and the position according to the match between White and Norway spruce references ###
# Splitting the blastn result file to run my R code
split -l 50000 blast_results_variants.out blast_results_variants_

# First step -> producing a file with matches between Norway spruce chromo:pos and our sites aligned to white spruce reference
R --vanilla <<EOF

## Blastn output with correspondance between white and norway pruce genomes for each site
data <- lapply(list.files("./", pattern = "blast_results_variants_"), read.table)

## Sites from the vcf file and correspondances with the newly produced ~600bp scaffolds
sites <- read.table("chrompos_vcf_WS.bed")[,1:2]
corresp <- read.table("near_variants.bed")
chrom_site_scaf <- data.frame(Chrom = sites[,1], Site = sites[,2], Scaffold = paste(corresp[,2], corresp[,3], sep = "-"), Crom_scaf = paste(corresp[,1], paste(corresp[,2], corresp[,3], sep = "-"), sep = ":"))
first_pos <- as.integer(unlist(lapply(strsplit(as.character(chrom_site_scaf$Scaffold), split = "-"), function(x) x[1])))
chrom_site_scaf$new_pos <- chrom_site_scaf$Site - first_pos

for (d in 1:length(data)){
  blastn <- data[[d]]
  colnames(blastn) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore", "qseq", "sseq")
  # The function getfasta (Bedtools) would have needed range like 0-500 and I gave it 1-500 meaning that 1 is 2 and so on...
  blastn$qstart <- blastn$qstart+1
  blastn$qend <- blastn$qend+1
  ## Ordering by scaffold and blasting score
  blastn_filtered <- blastn[-which(blastn$pident<80),]
  ## Finding the correspondance between norway and white spruce 
  TAB <- chrom_site_scaf
  TAB[,1] <- as.character(TAB[,1])
  TAB[,2] <- as.integer(TAB[,2])
  TAB[,3] <- as.character(TAB[,3])
  TAB[,4] <- as.character(TAB[,4])
  TAB[,5] <- as.integer(TAB[,5])
  for (i in 1:nrow(TAB)){
    hit <- which(TAB$Crom_scaf[i]==blastn_filtered$qseqid & TAB$new_pos[i]>=blastn_filtered$qstart & TAB$new_pos[i]<=blastn_filtered$qend)
    qgaps <- which(unlist(strsplit(as.character(blastn_filtered$qseq[hit[1]]), split = ""))=="-")
    sgaps <- which(unlist(strsplit(as.character(blastn_filtered$sseq[hit[1]]), split = ""))=="-")
    qgaps_inf <- length(which(qgaps<(TAB$new_pos[i]-blastn_filtered$qstart[hit[1]])))
    sgaps_inf <- length(which(sgaps<(TAB$new_pos[i]-blastn_filtered$qstart[hit[1]]+qgaps_inf)))
    if (length(hit)==0){TAB[i,2] <- -1}
    else if (blastn_filtered$sstart[hit[1]]<blastn_filtered$send[hit[1]]){
      pos <- TAB$new_pos[i]-blastn_filtered$qstart[hit[1]]+1+qgaps_inf+blastn_filtered$sstart[hit[1]]-1-sgaps_inf
      TAB[i,1] <- as.character(blastn_filtered$sseqid[hit[1]])
      TAB[i,2] <- as.numeric(pos) 
    }
    else if (blastn_filtered$sstart[hit[1]]>blastn_filtered$send[hit[1]]){
      pos <- -(TAB$new_pos[i]-blastn_filtered$qstart[hit[1]]+1+qgaps_inf)+blastn_filtered$sstart[hit[1]]+1+sgaps_inf
      TAB[i,1] <- as.character(blastn_filtered$sseqid[hit[1]])
      TAB[i,2] <- as.numeric(pos)  
    }
    print(i)
  }
  write.table(TAB, file = paste0("./vcf_NS_bis_", as.character(d), ".txt"))
  print(d)
}

## Combining the different tables of matching sites between Norway and WHite spruce
files <- list.files("./", pattern = "vcf_NS_bis")
TAB <- read.table(files[1])[,1:2]
TAB[,1] <- as.character(TAB[,1])
for(i in 2:length(files)){ 
  file <- read.table(files[i])[,1:2]
  file[,1] <- as.character(file[,1])
  TAB[which(TAB[,2]==-1),] <- file[which(TAB[,2]==-1),]
  print(i)
}
TAB[,3] <- paste(TAB[,1], as.character(TAB[,2]), sep = "_")
write.table(TAB, file = "./chrompos_vcf_NS.txt", quote = F, col.names = F, row.names = F, sep="\t")


EOF

## Second step -> replacing the first two column in the vcf by the good correspondances and filtering the numerous NA obtained
zcat ~/mydata/datashare/Spruce/exome_capture/WES_mapping/ANGSD/ref_Pglauca/Full_Sampling_intersect_poly.vcf.gz | grep -v "##" | awk 'NR==1' > Full_Sampling_intersect_poly_filtered_NS.vcf
zcat ~/mydata/datashare/Spruce/exome_capture/WES_mapping/ANGSD/ref_Pglauca/Full_Sampling_intersect_poly.vcf.gz | grep -v "#" | cut -f4- - | paste chrompos_vcf_NS.txt - | grep -v "JZKD" >> Full_Sampling_intersect_poly_filtered_NS.vcf


### Using SnpEff to find the annotations of the resulted vcf from the Norway spruce genome ###

## Building a new genome database
nano ~/TOOLS/snpEff/snpEff.config
## Adding 
    # Picea abies
    Pabies.genome : Picea_abies

mkdir ~/TOOLS/snpEff/data/Pabies
mv .... ~/TOOLS/snpEff/data/Pabies/sequences.fa.gz # genome
mv .... ~/TOOLS/snpEff/data/Pabies/genes.gff.gz # annotations

## We removed the lines containing CDS, kept the lines containing exon
zcat ~/TOOLS/snpEff/data/Pabies/genes.gff.gz | grep -v "CDS" > ~/TOOLS/snpEff/data/Pabies/genes.gff.gz

java -jar ~/TOOLS/snpEff/snpEff.jar build -gff3 -v Pabies

java -Xmx10g -jar ~/TOOLS/snpEff/snpEff.jar dump -v -txt Pabies > Pabies.txt

# Running snpEFF
java -Xmx10G -jar ~/TOOLS/snpEff/snpEff.jar -c ~/TOOLS/snpEff/snpEff.config -v Pabies ./Full_Sampling_vcf_filtered_NS.vcf > Red_Spruce_intersect_poly_snpeff.vcf



#################################
### Processing SNPeff outputs ###

# Number of variants included in each category
cat Red_Spruce_intersect_poly_snpeff.vcf | awk '{print $8}' | awk -F ';' '{print $6}' | awk -F '|' '{print $2}' | sort | uniq --count > SNPeff_SNPcategory.txt

# Number of variants where the reference allele (white spruce) does not match the genome (norway spruce)
cat Red_Spruce_intersect_poly_snpeff.vcf | grep "WARNING_REF_DOES_NOT_MATCH_GENOME" | wc -l # 97260


### Splitting the SFS according to the different groups ###

# Creating a file with variants scaffold:position matches between Norway and White Spruce
zcat ~/mydata/datashare/Spruce/exome_capture/WES_mapping/ANGSD/ref_Pglauca/Full_Sampling_intersect_poly.vcf.gz | grep -v "##" | cut -f -2 | grep -v "#" | paste -d "\t" - chrompos_vcf_NS.txt | grep -v "NA" > MatchingWhiteNorway_variants.txt

# List of silent variants 
cat Red_Spruce_intersect_poly_snpeff.vcf | grep -E "downstream_gene_variant|intergenic_region|intron_variant|synonymous_variant|upstream_gene_variant" | awk '{print $1,$2}' > list_silent.txt

cat Red_Spruce_intersect_poly_snpeff.vcf | awk '{print $8}' | awk -F ';' '{print $6}' | awk -F '|' '{print $2}' | grep -nrE "downstream_gene_variant|intergenic_region|intron_variant|synonymous_variant|upstream_gene_variant" | awk '{print $1,$2}' > list_silent.temp
cat list_silent.temp | awk -F ':' '{print $3}' > list_silent.txt
cat Red_Spruce_intersect_poly_snpeff.vcf | awk '{print $8}' | awk -F ';' '{print $6}' | awk -F '|' '{print $2}' | grep -nrE "missense_variant|splice_acceptor_variant|splice_donor_variant|splice_region_variant|start_lost|start_retained_variant|stop_gained|stop_lost|stop_retained_variant" | awk '{print $1,$2}' > list_nonsyn.temp
cat list_nonsyn.temp | awk -F ':' '{print $3}' > list_nonsyn.txt
cat Red_Spruce_intersect_poly_snpeff.vcf | awk '{print $8}' | awk -F ';' '{print $6}' | awk -F '|' '{print $2}' | grep -nrE "intron_variant|synonymous_variant" | awk '{print $1,$2}' > list_synonymous_intronic.temp
cat list_synonymous_intronic.temp | awk -F ':' '{print $3}' > list_synonymous_intronic.txt
cat Red_Spruce_intersect_poly_snpeff.vcf | awk '{print $8}' | awk -F ';' '{print $6}' | awk -F '|' '{print $2}' | grep -nrE "synonymous_variant" | awk '{print $1,$2}' > list_synonymous.temp
cat list_synonymous.temp | awk -F ':' '{print $3}' > list_synonymous.txt
cat Red_Spruce_intersect_poly_snpeff.vcf | awk '{print $8}' | awk -F ';' '{print $6}' | awk -F '|' '{print $2}' | grep -nrE "intergenic" | awk '{print $1,$2}' > list_intergenic.temp
cat list_intergenic.temp | awk -F ':' '{print $3}' > list_intergenic.txt

rm list_silent.temp
rm list_nonsyn.temp
rm list_synonymous.temp
rm list_synonymous_intronic.temp
rm list_intergenic.temp

### SFS with only silent or non-synonymous sites ###

# Initial SFS
~/mydata/datashare/Spruce/exome_capture/WES_mapping/ANGSD/ref_Pglauca/REGIONS/CORE_MARGIN_EDGE_intersect.sfs

R --vanilla <<EOF
	
	SFS <- read.table("~/mydata/datashare/Spruce/exome_capture/WES_mapping/ANGSD/ref_Pglauca/REGIONS/CORE_MARGIN_EDGE_intersect.sfs")
	
	match <- read.table("MatchingWhiteNorway_variants.txt")[,-5]
	colnames(match) <- c("Scaffold_WS","Position_WS","Scaffold_NS","Position_NS")
	silent <- read.table("./list_silent.txt")
	nonsyn <- read.table("./list_nonsyn.txt")
	synonymous <- read.table("./list_synonymous.txt")
	synonymous_intronic <- read.table("./list_synonymous_intronic.txt")
	intergenic <- read.table("./list_intergenic.txt")
	
	chrom_pos_silent <- match[which(paste0(match$Scaffold_NS, "_", match$Position_NS)%in%paste0(silent$V1, "_", silent$V2)),]
	SFS_silent <- SFS[which(paste0(SFS[,3], "_", SFS[,4])%in%paste0(chrom_pos_silent$Scaffold_WS, "_", chrom_pos_silent$Position_WS)),]
	write.table(SFS_silent, file = "./SFS_silent.sfs", quote = F, col.names = F, row.names = F, sep="\t")
	
	chrom_pos_nonsyn <- match[which(paste0(match$Scaffold_NS, "_", match$Position_NS)%in%paste0(nonsyn$V1, "_", nonsyn$V2)),]
	SFS_nonsyn <- SFS[which(paste0(SFS[,3], "_", SFS[,4])%in%paste0(chrom_pos_nonsyn$Scaffold_WS, "_", chrom_pos_nonsyn$Position_WS)),]
	write.table(SFS_nonsyn, file = "./SFS_nonsyn.sfs", quote = F, col.names = F, row.names = F, sep="\t")
	
	chrom_pos_synonymous_intronic <- match[which(paste0(match$Scaffold_NS, "_", match$Position_NS)%in%paste0(synonymous_intronic$V1, "_", synonymous_intronic$V2)),]
	SFS_synonymous_intronic <- SFS[which(paste0(SFS[,3], "_", SFS[,4])%in%paste0(chrom_pos_synonymous_intronic$Scaffold_WS, "_", chrom_pos_synonymous_intronic$Position_WS)),]
	write.table(SFS_synonymous_intronic, file = "./SFS_synonymous_intronic.sfs", quote = F, col.names = F, row.names = F, sep="\t")
		
	chrom_pos_synonymous <- match[which(paste0(match$Scaffold_NS, "_", match$Position_NS)%in%paste0(synonymous$V1, "_", synonymous$V2)),]
	SFS_synonymous <- SFS[which(paste0(SFS[,3], "_", SFS[,4])%in%paste0(chrom_pos_synonymous$Scaffold_WS, "_", chrom_pos_synonymous$Position_WS)),]
	write.table(SFS_synonymous, file = "./SFS_synonymous.sfs", quote = F, col.names = F, row.names = F, sep="\t")

	chrom_pos_intergenic <- match[which(paste0(match$Scaffold_NS, "_", match$Position_NS)%in%paste0(intergenic$V1, "_", intergenic$V2)),]
	SFS_intergenic <- SFS[which(paste0(SFS[,3], "_", SFS[,4])%in%paste0(chrom_pos_intergenic$Scaffold_WS, "_", chrom_pos_intergenic$Position_WS)),]
	write.table(SFS_intergenic, file = "./SFS_intergenic.sfs", quote = F, col.names = F, row.names = F, sep="\t")
EOF

# Transform for dadi
~/mydata/Thibaut/RedSpruce_demography/DADI/ref_Pglauca/realsfs2dadi.pl SFS_intergenic.sfs 178 51 110 > RedSpruce_intergenic.sfs

