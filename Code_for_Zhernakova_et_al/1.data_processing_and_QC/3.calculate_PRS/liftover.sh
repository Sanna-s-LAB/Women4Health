plink_f=$1
to=$2
chain_f=$3

awk 'BEGIN {FS=" "}; {OFS="\t"}; {print $1, $4-1, $4, $2}' \
${plink_f}.bim > ${plink_f}.tmp.positions.bed

CrossMap bed $chain_f ${plink_f}.tmp.positions.bed ${plink_f}.tmp.positions.${to}.bed

cut -f4 ${plink_f}.tmp.positions.${to}.bed > ${plink_f}.tmp.successful_snps.txt
awk '$1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {print $4, $1, $3}' \
  ${plink_f}.tmp.positions.${to}.bed > ${plink_f}.tmp.pos_update.txt

mv ${plink_f}.bim ${plink_f}.original.bim

# Update positions in bim file using awk
awk 'BEGIN {OFS="\t"} 
     NR==FNR {pos[$1]=$3; chr[$1]=$2; next} 
     $2 in pos {$1=chr[$2]; $4=pos[$2]} 
     {print}' ${plink_f}.tmp.pos_update.txt ${plink_f}.original.bim > ${plink_f}.bim

plink2 --bfile ${plink_f} --make-pgen --sort-vars --out ${plink_f}.tmp.${to}
plink2 --pfile ${plink_f}.tmp.${to} --make-bed --out ${plink_f}.${to} --extract ${plink_f}.tmp.successful_snps.txt

mv ${plink_f}.original.bim ${plink_f}.bim
rm ${plink_f}.tmp*

