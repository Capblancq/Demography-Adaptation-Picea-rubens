### Whole exome sequencing reveals a long-term decline in effective population size of red spruce (Picea rubens)
---------------------

This repository gathers the scripts used to conduct the analyses and produce the figures of the study described in Capblancq et al. 2020, Evolutionary Applications (DOI: ??):

1. ReferenceProduction.sh describes the way we created the reduced reference used for mapping the reads

2. Mapping_reads.sh describes how we mapped the reads and produced the quality statistics

3. Genotype_likelihoods_ANGSD.sh describes how we used ANGSD (http://www.popgen.dk/angsd/index.php/ANGSD#Overview) to produce genotype likelihoods from the .bam files

4. Downstream_anayses_ANGSD.sh describes how we used the genotype likelihoods to produce different estimates of genetic diversity (thetas), genetic structure (PCA and Fst) and linkage disiquilibrium (LD) along the genome.

5. Annotations.sh describes how we annotated the resulting vcf file using the norway spruce reference genome and annotations.
