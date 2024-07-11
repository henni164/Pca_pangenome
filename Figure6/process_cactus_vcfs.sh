LIST=($(cut -f 1 haps.txt))

for f in ${LIST[@]}; do sed 's/#CHROM/CHROM/g' Aus_pangenome_cactus.${f}.full.vcf | grep -v "#" | cut -f 3,10-40 > ${f}_setup.txt; done

for f in ${LIST[@]}; do grep -v "#" Aus_pangenome_cactus.${f}.full.vcf | awk -v OFS="\t" -v HAP=$f '{print HAP,$1,$2,$3}' - >> all_varpos_index.txt; done
