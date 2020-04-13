################################
##### Reference production #####

#!/bin/bash

### Creating the reference for WES data mapping

# Listing the scaffolds of the Picea glauca genome containing at least one of the probes used for WES 
cat UVM_131901_RG_0806_probes_ref_WS77111.bed | awk '{print $1}' | sort | uniq > list.txt

# Checking the number of lines
wc -l list.txt

# Listing the sequences ID depending on the scaffolds selected
grep -f list.txt WS77111_V1_genomic.fna | awk '{print $1}' | cut -c 2- > ID_list.txt

# Checking if the number of lines is ok
wc -l ID_list.txt

# Extracting the selected scaffolds within genome
samtools faidx WS77111_V1_genomic.fna `cat ID_list.txt` > WES_RedSpurce_reference.fna

# Check if the number of sequences fit
grep -c "^>" WES_RedSpurce_reference.fna

# or maybe
xargs samtools faidx WS77111_V1_genomic.fna.gz < ID_list.txt
