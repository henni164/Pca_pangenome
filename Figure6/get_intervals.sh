module load bedtools

LIST2=($(cut -f 1 haplist.txt))

for f in ${LIST2[@]}; do
	LIST=($(cut -f 4 ${f}_intervals.txt | sort | uniq))

	for g in ${LIST[@]}; do grep "${g}	" ${f}_intervals.txt | sed $"s/\r//g" > ${f}_${g}_intervals.txt; done

	for g in ${LIST[@]}; do sort -k1,1 -k2,2n ${f}_${g}_intervals.txt > ${f}_${g}_sorted.txt; done

	for g in ${LIST[@]}; do bedtools merge -s -d 1 -c 4,5,6 -o distinct,mean,distinct -i ${f}_${g}_sorted.txt > ${f}_${g}_merged.txt; done

	for g in ${LIST[@]}; do grep "+" ${f}_${g}_merged.txt |  awk -v OFS="\t" '{if($3 - $2 >= 50000) print $0}' | sed "s/\t+/\t${f}/g" > ${f}_${g}_shared.txt; done

done

LIST=($(cut -f 1 haplist_extend.txt))
index=29
for g in ${LIST[@]}; do
	echo ${index}
	if [ ${g} == "hap3" ] || [ ${g} == "hap4" ]; then
		:
	else
		cp hap3_${g}_shared.txt all_int_${g}.txt
		echo "doing regions from haplotypes ${LIST2[@]:index} on ${g} right now"
		for ((i=${#LIST2[@]}-1; i>=$index; i--))
		do
			f="${LIST2[$i]}"
			if [ ${f} == ${g} ]; then
				:
			else
				bedtools subtract -a ${f}_${g}_shared.txt -b all_int_${g}.txt > ${f}_${g}_subtracts.txt
				cat all_int_${g}.txt ${f}_${g}_subtracts.txt > tmp && mv tmp all_int_${g}.txt
			fi
		done
		let index=${index}-1
	fi

done

cat all_int_* > all_intervals.txt

awk -v OFS="\t" '{if($3 - $2 >= 50000) print $0}' all_intervals.txt > tmp && mv tmp all_intervals.txt
