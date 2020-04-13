###################################
####### GENOTYPE LIKELIHOODS ######
###################################

#!/bin/bash

##########################################
##### Genotype likelihood per REGION #####

mkdir ./REGIONS
mkdir ./REGIONS/CORE
mkdir ./REGIONS/EDGE
mkdir ./REGIONS/MARGIN

# Have to create a list of .bam files for each region
# ./REGIONS/CORE_bam.list
# ./REGIONS/MARGIN_bam.list
# ./REGIONS/EDGE_bam.list

# Estimating the genotype likelihoods and SFS for each regional population separately with all the sites (both mono- and poly-morphic sites by removing the SNP_val option)
N=3
for POP in CORE MARGIN EDGE
do
    (
	count=`cat ./REGIONS/${POP}_bam.list | wc -l`
	MAX=$(($count+$count+$count+$count+$count))
	~/TOOLS/angsd/angsd -b ./REGIONS/${POP}_bam.list -ref ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -anc ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -out ./REGIONS/${POP}/${POP} -nThreads 4 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 2 -setMinDepthInd 2 -setMaxDepthInd 17 -setMinDepth 15 -setMaxDepth ${MAX} -skipTriallelic 1 -GL 1 -doHWE 1 -doGeno 32 -doPost 1 -doCounts 1 -doMaf 1 -doMajorMinor 1 -doGlf 1
    echo "starting task $POP.."
    sleep $(( (RANDOM % 3) + 1))
    ) &
    # allow only to execute $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
        # wait only for first job
        wait -n
    fi
done

# Finding loci intersection among regions
zcat ./REGIONS/CORE/CORE.mafs.gz | awk 'BEGIN { OFS = ":" }{ print $1,$2 }' | sed '1d' | sort > ./REGIONS/CORE/CORE_sites.txt
zcat ./REGIONS/EDGE/EDGE.mafs.gz | awk 'BEGIN { OFS = ":" }{ print $1,$2 }' | sed '1d' | sort > ./REGIONS/EDGE/EDGE_sites.txt
zcat ./REGIONS/MARGIN/MARGIN.mafs.gz | awk 'BEGIN { OFS = ":" }{ print $1,$2 }' | sed '1d' | sort > ./REGIONS/MARGIN/MARGIN_sites.txt

comm -12 ./REGIONS/CORE/CORE_sites.txt ./REGIONS/EDGE/EDGE_sites.txt > ./REGIONS/CORE_EDGE_sites.txt
comm -12 ./REGIONS/CORE_EDGE_sites.txt ./REGIONS/MARGIN/MARGIN_sites.txt > ./REGIONS/CORE_EDGE_MARGIN_sites.txt

sed 's/:/\t/' ./REGIONS/CORE_EDGE_MARGIN_sites.txt | sort -b -k1,1 > ./REGIONS/intersect.txt
cut -f1 ./REGIONS/intersect.txt | uniq | sort > ./REGIONS/intersect.chrs
~/TOOLS/angsd/angsd sites index ./REGIONS/intersect.txt

# Number of total sites with monomorphic and plolymorphic sites
cat ./REGIONS/intersect.txt | wc -l # 42,328,740 sites

# Re-launching the genotype likelihoods estimation by keeping only intersecting sites (monomorphic and polymorphic) and making the optimization step on the SFS
N=3
for POP in CORE EDGE MARGIN 
do
    (
	~/TOOLS/angsd/angsd -b ./REGIONS/${POP}_bam.list -GL 1 -out ./REGIONS/${POP}/${POP}_intersect -ref ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -anc ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/WES_RedSpurce_Pglauca_reduced_reference.fna -sites ./REGIONS/intersect.txt -rf ./REGIONS/intersect.chrs -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -skipTriallelic 1 -doHWE 1 -doGeno 32 -doPost 1 -doCounts 1 -doMaf 1 -doMajorMinor 1 -doGlf 1 -doSaf 1 -fold 0 
    #EM optimization of the sfs
	~/TOOLS/angsd/misc/realSFS ./REGIONS/${POP}/${POP}_intersect.saf.idx -maxIter 50000 -tole 1e-6 -P 4 > ./REGIONS/${POP}/${POP}_intersect.sfs
    echo "starting task $POP.."
    sleep $(( (RANDOM % 3) + 1))
    ) &
    # allow only to execute $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
        # wait only for first job
        wait -n
    fi
done
