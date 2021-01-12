##################################
####### DOWNSTREAM ANALYSES ######
##################################

#!/bin/bash

##############################
##### Thetas, Tajima's D #####

# Estimating thetas from bam files using the optimized SFS as a prior and keeping only the intersecting sites
N=3
for POP in CORE EDGE MARGIN 
do
    (
	~/TOOLS/angsd/angsd -bam ./REGIONS/${POP}_bam.list -out ./REGIONS/${POP}/${POP}_intersect -doThetas 1 -doSaf 1 -pest ./REGIONS/${POP}/${POP}_intersect.sfs -anc ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -ref ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -GL 1 -nThreads 6 -sites ./REGIONS/intersect.txt -rf ./REGIONS/intersect.chrs -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20
	~/TOOLS/angsd/misc/thetaStat do_stat ./REGIONS/${POP}/${POP}_intersect.thetas.idx
    sleep $(( (RANDOM % 3) + 1))
    ) &
    # allow only to execute $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
        # wait only for first job
        wait -n
    fi
done


##################################
##### Linkage disequilibrium #####

### LD per region

# Genotype likelihoods for intersected and polymorphic sites (beagle) 
N=3
for POP in MARGIN CORE EDGE 
do
    (
    ~/TOOLS/angsd/angsd -b ./REGIONS/${POP}_bam.list -GL 1 -out ./REGIONS/${POP}/${POP}_intersect_poly -ref ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -anc ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -sites ./REGIONS/intersect.txt -rf ./REGIONS/intersect.chrs -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -skipTriallelic 1 -SNP_pval 1e-6 -doMaf 1 -doMajorMinor 1 -doGlf 2
	# Estimating LD
	zcat REGIONS/${POP}/${POP}_intersect_poly.mafs.gz | awk '{print $1, $2}' > ./REGIONS/${POP}/pos_${POP}_intersect_poly.txt
	count=`cat ./REGIONS/${POP}_bam.list | wc -l`
    N_SITES=$((`zcat ./REGIONS/${POP}/${POP}_intersect_poly.mafs.gz | wc -l`-1))
	~/TOOLS/ngsLD/ngsLD --geno ./REGIONS/${POP}/${POP}_intersect_poly.beagle.gz --probs --n_ind ${count} --n_sites ${N_SITES} --pos ./REGIONS/${POP}/pos_${POP}_intersect_poly.txt --out ./REGIONS/${POP}/${POP}_intersect_maf0.5.LD --min_maf 0.05
    sleep $(( (RANDOM % 3) + 1))
    ) &
    # allow only to execute $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
        # wait only for first job
        wait -n
    fi
done


R --vanilla <<EOF

	# Margin data
	data_margin <- read.table("./REGIONS/MARGIN/MARGIN_intersect_maf0.5.LD", nrows = 1000000)
	data_margin <- data_margin[-which(as.character(data_margin$V1)=="chromo"),]

	# Edge data
	data_edge <- read.table("./REGIONS/EDGE/EDGE_intersect_maf0.5.LD", nrows = 1000000)
	data_edge <- data_edge[-which(as.character(data_edge$V1)=="chromo"),]

	# Core data
	data_core <- read.table("./REGIONS/CORE/CORE_intersect_maf0.5.LD", nrows = 1000000)
	data_core <- data_core[-which(data_core$V1=="chromo"),]

	####################################
	### Fitting an exponential decay ###
	
	# Function of LD decay (from https://www.pnas.org/content/98/20/11479)

	### CORE ###
	# Finding the initial parameters a and b
	Cstart <- c(C=0.1)
	n = 178
	modelC <- nls(LD ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))), data = data.frame(dist=data_core[,5],LD=data_core[,6]), start = Cstart, control=nls.control(maxiter=100)) 
	# extract rho, the recombination parameter, 4Nr
	new.rho <- summary(modelC)$parameters[1]
	# feed in the new value of rho to obtain LD values adjusted for their distances along the chromosome/genome
	fpoints<-((10+new.rho*data_core[,5])/((2+new.rho*data_core[,5])*(11+new.rho*data_core[,5])))*(1+((3+new.rho*data_core[,5])*(12+12*new.rho*data_core[,5]+(new.rho*data_core[,5])^2))/(n*(2+new.rho*data_core[,5])*(11+new.rho*data_core[,5])))
	# final table
	ld.df <- data.frame(distance=data_core[,5],fpoints,rep("Core",length(data_core[,5])))
	colnames(ld.df) <- c("distance","fpoints","Region")
	ld.df_core<-ld.df[order(ld.df$distance),]

	### MARGIN ###
	Cstart <- c(C=0.1)
	n = 51
	modelC <- nls(LD ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))), data = data.frame(dist=data_margin[,5],LD=data_margin[,6]), start = Cstart, control=nls.control(maxiter=100)) 
	new.rho <- summary(modelC)$parameters[1]
	fpoints<-((10+new.rho*data_margin[,5])/((2+new.rho*data_margin[,5])*(11+new.rho*data_margin[,5])))*(1+((3+new.rho*data_margin[,5])*(12+12*new.rho*data_margin[,5]+(new.rho*data_margin[,5])^2))/(n*(2+new.rho*data_margin[,5])*(11+new.rho*data_margin[,5])))
	ld.df <- data.frame(distance=data_margin[,5],fpoints,rep("Margin",length(data_margin[,5])))
	colnames(ld.df) <- c("distance","fpoints","Region")
	ld.df_margin<-ld.df[order(ld.df$distance),]

	### EDGE ###
	Cstart <- c(C=0.1)
	n = 110
	modelC <- nls(LD ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))), data = data.frame(dist=data_edge[,5],LD=data_edge[,6]), start = Cstart, control=nls.control(maxiter=100)) 
	new.rho <- summary(modelC)$parameters[1]
	fpoints<-((10+new.rho*data_edge[,5])/((2+new.rho*data_edge[,5])*(11+new.rho*data_edge[,5])))*(1+((3+new.rho*data_edge[,5])*(12+12*new.rho*data_edge[,5]+(new.rho*data_edge[,5])^2))/(n*(2+new.rho*data_edge[,5])*(11+new.rho*data_edge[,5])))
	ld.df <- data.frame(distance=data_edge[,5],fpoints,rep("Edge",length(data_edge[,5])))
	colnames(ld.df) <- c("distance","fpoints","Region")
	ld.df_edge<-ld.df[order(ld.df$distance),]

	##############
	#### Plot ####

	cols <- c('Core'="#FEDF00", 'Margin'="#377EB8", 'Edge'="#4DAF4A")

	pdf("LD_decay_intersect_maf0.5.perRegion.pdf")
	plot(ld.df_core$distance,ld.df_core$fpoints, xlim=c(0,500), type = "n", xlab = "Distance (bp)", ylab = "LD (r2)")
	lines(ld.df_core$distance,ld.df_core$fpoints, lwd=2,col=cols[1], xlim=c(0,500))
	lines(ld.df_margin$distance,ld.df_margin$fpoints, lwd=2,col=cols[2], xlim=c(0,500))
	lines(ld.df_edge$distance,ld.df_edge$fpoints, lwd=2,col=cols[3], xlim=c(0,500))
	legend('topright','groups',c("Core","Margin","Edge"), lwd = 2, col=cols) #, ncol=5,nrow=2,bty ="n")
	dev.off()

EOF


####################################
##### PCA on the full sampling #####

### Pruning the full sampling dataset

# List of .bam files used
ls ~/mydata/Thibaut/WES_mapping/Mapping_WES/ref_Pglauca/BWA/*_rmdup.sorted.bam > ./all_bam.list

# Estimating genotype likelihood for the polymorphic sites common in the three regional populations
~/TOOLS/angsd/angsd -b ./all_bam.list -ref ~/mydata/datashare/Spruce/exome_capture/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -anc ~/mydata/datashare/Spruce/exome_capture/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -fai ~/mydata/datashare/Spruce/exome_capture/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna.fai -out ./Full_Sampling_intersect_poly -sites ./REGIONS/intersect.txt -rf ./REGIONS/intersect.chrs -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -skipTriallelic 1 -SNP_pval 1e-6 -minMaf 0.003 -GL 1 -doGlf 1 -doCounts 1 -doPost 1 -doMajorMinor 1 -dovcf 1 -doHWE 1 -doGeno 2 -doMaf 1 -postCutoff 0.8

# 1,372,627 polymorphic sites

# LD estimation
zcat ./Full_Sampling_intersect_poly.mafs.gz | awk '{print $1, $2}' > ./pos_Full_Sampling_intersect_poly.txt
count=`cat ./all_bam.list | wc -l`
N_SITES=$((`zcat ./Full_Sampling_intersect_poly.mafs.gz | wc -l`-1))
~/TOOLS/ngsLD/ngsLD --geno ./Full_Sampling_intersect_poly.beagle.gz --probs --n_ind ${count} --n_sites ${N_SITES} --pos ./pos_Full_Sampling_intersect_poly.txt --out ./Full_Sampling_intersect_poly_maf0.5.LD --min_maf 0.05

# Selecting only one polymorphic site every 500bp
R --vanilla <<EOF

	sites <- read.table("./pos_Full_Sampling_intersect_poly.txt", header=T)
	list_sites <- split(sites, as.factor(sites[,1]))
	pruned_sites <- list()
	for(i in 1:length(list_sites)){
		vec <- as.integer(list_sites[[i]][,2])
		res <- NULL
		j <- min(vec)
		while(j < max(vec)){
			sup <- vec[vec>j][1]
			res <- c(res,sup)
			j <- sup+500
		}
		scaffold <- rep(as.character(list_sites[[i]][1,1]), length(res))
		pruned_sites[[i]] <- cbind(scaffold,res)
	}
	list_pruned_sites <- do.call(rbind,pruned_sites)
	colnames(list_pruned_sites) <- colnames(sites)
	write.table(list_pruned_sites, "./Full_Sampling_pruned_sites.txt", quote = FALSE, sep = "\t", row.names=FALSE, col.names=FALSE)
	write.table(list_pruned_sites[,1], "./Full_Sampling_pruned_sites.chr", quote = FALSE, sep = "\t", row.names=FALSE, col.names=FALSE)
	
EOF
# Indexing the list of sites
~/TOOLS/angsd/angsd sites index ./Full_Sampling_pruned_sites.txt

# 114,699 polymorphic sites after pruning

### PCA on the full sampling pruned dataset

# Estimating genotype likelihood and covariance matrix for the prunned dataset
~/TOOLS/angsd/angsd -b ./all_bam.list -ref ~/mydata/datashare/Spruce/exome_capture/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -anc ~/mydata/datashare/Spruce/exome_capture/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -fai ~/mydata/datashare/Spruce/exome_capture/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna.fai -out ./Full_Sampling_intersect_prunned -sites Full_Sampling_pruned_sites.txt -rf Full_Sampling_pruned_sites.chrs -nThreads 2 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -skipTriallelic 1 -SNP_pval 1e-6 -minMaf 0.003 -GL 1 -doHWE 1 -doGeno 2 -doPost 1 -doCounts 1 -doMaf 1 -doMajorMinor 1 -doGlf 1 -dovcf 1 -nInd ${count} -postCutoff 0.8 -doIBS 1 -doCov 1 -makeMatrix 1

# Anylysing the covariance matrix and plotting the results in R, see PCA_redspruce.R script


######################################
### Fst among regional populations ###

### Pruning the dataset

# Re-launching the genotype likelihoods estimation by keeping only intersecting and prunned sites
N=3
for POP in CORE EDGE MARGIN 
do
    (
	~/TOOLS/angsd/angsd -b ./REGIONS/${POP}_bam.list -GL 1 -out ./REGIONS/${POP}/${POP}_intersect_prunned -ref ~/mydata/datashare/Spruce/exome_capture/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -anc ~/mydata/datashare/Spruce/exome_capture/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -sites Full_Sampling_pruned_sites.txt -rf Full_Sampling_pruned_sites.chrs -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -skipTriallelic 1 -doHWE 1 -doGeno 32 -doPost 1 -doCounts 1 -doMaf 1 -doMajorMinor 1 -doGlf 1 -doSaf 1 -fold 0 
    #EM optimization of the sfs
	~/TOOLS/angsd/misc/realSFS ./REGIONS/${POP}/${POP}_intersect_prunned.saf.idx -maxIter 50000 -tole 1e-6 -P 4 > ./REGIONS/${POP}/${POP}_intersect_prunned.sfs
    echo "starting task $POP.."
    sleep $(( (RANDOM % 3) + 1))
    ) &
    # allow only to execute $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
        # wait only for first job
        wait -n
    fi
done


# Estimating the 2d SFS with the pruned datasets
~/TOOLS/angsd/misc/realSFS ./REGIONS/CORE/CORE_intersect_prunned.saf.idx ./REGIONS/EDGE/EDGE_intersect_prunned.saf.idx -P 4 > ./REGIONS/CORE.EDGE_intersect_prunned.ml
~/TOOLS/angsd/misc/realSFS ./REGIONS/CORE/CORE_intersect_prunned.saf.idx ./REGIONS/MARGIN/MARGIN_intersect_prunned.saf.idx -P 4 > ./REGIONS/CORE.MARGIN_intersect_prunned.ml
~/TOOLS/angsd/misc/realSFS ./REGIONS/MARGIN/MARGIN_intersect_prunned.saf.idx ./REGIONS/EDGE/EDGE_intersect_prunned.saf.idx -P 4 > ./REGIONS/MARGIN.EDGE_intersect_prunned.ml

# Prepare the fst for each pair of region
~/TOOLS/angsd/misc/realSFS fst index ./REGIONS/CORE/CORE_intersect_prunned.saf.idx ./REGIONS/EDGE/EDGE_intersect_prunned.saf.idx -sfs ./REGIONS/CORE.EDGE_intersect_prunned.ml -fstout ./REGIONS/CORE.EDGE_intersect_prunned_Fst
~/TOOLS/angsd/misc/realSFS fst index ./REGIONS/CORE/CORE_intersect_prunned.saf.idx ./REGIONS/MARGIN/MARGIN_intersect_prunned.saf.idx -sfs ./REGIONS/CORE.MARGIN_intersect_prunned.ml -fstout ./REGIONS/CORE.MARGIN_intersect_prunned_Fst
~/TOOLS/angsd/misc/realSFS fst index ./REGIONS/MARGIN/MARGIN_intersect_prunned.saf.idx ./REGIONS/EDGE/EDGE_intersect_prunned.saf.idx -sfs ./REGIONS/MARGIN.EDGE_intersect_prunned.ml -fstout ./REGIONS/MARGIN.EDGE_intersect_prunned_Fst

# Get the global estimate (Unweighted and weighted Fst)
~/TOOLS/angsd/misc/realSFS fst stats ./REGIONS/CORE.EDGE_intersect_prunned_Fst.fst.idx # 0.017217 0.022925
~/TOOLS/angsd/misc/realSFS fst stats ./REGIONS/CORE.MARGIN_intersect_prunned_Fst.fst.idx # 0.023190 0.027451
~/TOOLS/angsd/misc/realSFS fst stats ./REGIONS/MARGIN.EDGE_intersect_prunned_Fst.fst.idx # 0.031392 0.027267


#######################################################
##### Preparing a 3D SFS for demographic analysis #####

## 3DSFS estimate for demographic analysis
# Intersected sites
~/TOOLS/angsd/misc/realSFS dadi -P 8 ./REGIONS/CORE/CORE_intersect.saf.idx ./REGIONS/MARGIN/MARGIN_intersect.saf.idx ./REGIONS/EDGE/EDGE_intersect.saf.idx -sfs ./REGIONS/CORE/CORE_intersect.sfs -sfs ./REGIONS/MARGIN/MARGIN_intersect.sfs -sfs ./REGIONS/EDGE/EDGE_intersect.sfs -ref ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -anc ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -nSites 10000000 > ./REGIONS/CORE_MARGIN_EDGE_intersect.sfs
# Pruned sites
~/TOOLS/angsd/misc/realSFS dadi -P 8 ./REGIONS/CORE/CORE_intersect_prunned.saf.idx ./REGIONS/MARGIN/MARGIN_intersect_prunned.saf.idx ./REGIONS/EDGE/EDGE_intersect_prunned.saf.idx -sfs ./REGIONS/CORE/CORE_intersect_prunned.sfs -sfs ./REGIONS/MARGIN/MARGIN_intersect_prunned.sfs -sfs ./REGIONS/EDGE/EDGE_intersect_prunned.sfs -ref ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -anc ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -nSites 10000000 > ./REGIONS/CORE_MARGIN_EDGE_intersect_prunned.sfs

