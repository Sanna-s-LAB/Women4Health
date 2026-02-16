
##
### target data QC
##
# 
# cd /mnt/sannaLAB-Temp/dasha/PRS/data/microarray
# 
# plink2 \
#     --bfile W4H_forMergingWith16S_geno01_hwe5e-6 \
#     --chr 1-22,X,Y \
#     --maf 0.01 \
#     --hwe 1e-6 \
#     --geno 0.01 \
#     --mind 0.01 \
#     --make-bed \
#     --out microarray.QC 
#     #--set-all-var-ids @:# 
# 
# sed -i "s:GSA-::g" microarray.QC.bim
# 
# sh ../../scripts/liftover.sh microarray.QC hg19 ../ref_data/hg38ToHg19.over.chain
# 
# plink2 \
# --bfile microarray.QC.hg19 \
# --rm-dup force-first list \
# --make-bed \
# --out microarray.QC.hg19.rmdup
# 
# plink2 \
# --bfile microarray.QC.hg19.rmdup \
# --freq '\
# --out microarray.QC.hg19.rmdup



##
### Base data QC
##
cd /mnt/sannaLAB-Temp/dasha/PRS/data/

protein=$1
prot_id=`awk -F'\t' -v prot="$protein" '$3 == prot {print $1}' /mnt/sannaLAB-Temp/dasha/PRS/data/olink_protein_names_UKB.txt`

echo "Processing $protein, $prot_id"

base_data=pQTLs/combined_${prot_id}_UKBB_proteomics_female_only.tsv.gz 
target_data=microarray/microarray.batch12.QC.hg19
ld_ref=ref_data/1kg_mrcieu/EUR
out_prefix=../results_b12/Plink/${protein}


# subset the base data to keep only SNPs present in the target data and remove duplicates (keep the one that has the same allele as in the target data)
Rscript ../scripts/filter_summary_stats.R $protein $base_data ${target_data}.bim rsid 

plink2 \
    --bfile ${ld_ref} \
    --clump-p1 1e-4 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump ${base_data%gz}in_microarray.rmdup.txt \
    --clump-snp-field rsid \
    --clump-p-field p_value \
    --out ${out_prefix}.tmp_clumped

awk 'NR!=1{print $3}' ${out_prefix}.tmp_clumped.clumps >  ${out_prefix}.tmp_clumped.valid.snp

### Plink ###

echo "5e-08 0 5e-08" > ../results_b12/Plink/range_list 
#echo "0.00001 0 0.00001" >> ../results_b12/Plink/range_list
#echo "0.0001 0 0.0001" >> ../results_b12/Plink/range_list
#echo "0.001 0 0.001" >> ../results_b12/Plink/range_list

awk '{print $10,$8}' ${base_data%gz}in_microarray.rmdup.txt > ${base_data%gz}in_microarray.pvalues


plink2 \
    --bfile $target_data \
    --score ${base_data%gz}in_microarray.rmdup.txt 10 3 5 header \
    --q-score-range ../results_b12/Plink/range_list ${base_data%gz}in_microarray.pvalues \
    --extract ${out_prefix}.tmp_clumped.valid.snp \
    --read-freq ${target_data}.afreq \
    --out ${out_prefix}

rm ${out_prefix}.tmp_clumped.*
rm ${base_data%gz}in_microarray.pvalues
rm ${base_data%gz}in_microarray.rmdup.txt



