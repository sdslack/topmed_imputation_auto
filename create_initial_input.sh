#!/bin/bash

#Set arguments
if [ "$#" -eq  "0" ]
then
   echo "Usage: ${0##*/} <plink_prefix> <out_dir> <code_dir>"
   echo "Script removes SNPs with duplicate positions, removes"
   echo "strand ambiguous SNPs, and performs pre-imputation QC."
   echo "Current script runs PLINK1.9 with --keep-allele-order."
   echo "Crossover is from hg19 to hg38."
   exit
fi

plink_prefix=$1
out_dir=$2
code_dir=$3

#Remove SNPs with duplicate positions
plink --bfile $plink_prefix \
   --list-duplicate-vars suppress-first \
   --keep-allele-order \
   --out ${out_dir}/tmp_dupl_check
   
cat ${out_dir}/tmp_dupl_check.dupvar | sed -e '1d' | \
   cut -f4 > ${out_dir}/tmp_dupl_snpids.txt
plink --bfile $plink_prefix \
   --exclude ${out_dir}/tmp_dupl_snpids.txt \
   --keep-allele-order \
   --make-bed --out ${out_dir}/tmp_no_dupl


#Create bed file to crossover from hg19 to hg38 
cat ${out_dir}/tmp_no_dupl.bim | cut -f1 | sed 's/^/chr/' > ${out_dir}/tmp_c1.txt
cat ${out_dir}/tmp_no_dupl.bim | cut -f4 > ${out_dir}/tmp_c2.txt
cat ${out_dir}/tmp_no_dupl.bim | cut -f4 > ${out_dir}/tmp_c3.txt
cat ${out_dir}/tmp_no_dupl.bim | cut -f2 > ${out_dir}/tmp_c4.txt
paste  ${out_dir}/tmp_c1.txt \
       ${out_dir}/tmp_c2.txt \
       ${out_dir}/tmp_c3.txt \
       ${out_dir}/tmp_c4.txt \
       >  ${out_dir}/tmp_in.bed

#Do crossover
CrossMap.py bed ${code_dir}/hg19ToHg38.over.chain \
   ${out_dir}/tmp_in.bed  \
   ${out_dir}/tmp_out.bed

#Extract only those SNPs that were successfully cross-overed
cut -f4 ${out_dir}/tmp_out.bed > ${out_dir}/tmp_snp_keep.txt
plink --bfile ${out_dir}/tmp_no_dupl \
   --extract ${out_dir}/tmp_snp_keep.txt \
   --keep-allele-order \
   --make-bed --out ${out_dir}/tmp_gwas

#Update bim file positions
Rscript --vanilla ${code_dir}/update_pos.R \
  ${out_dir}/tmp_out.bed ${out_dir}/tmp_gwas.bim

#Remove strand ambiguous SNPs
Rscript --vanilla ${code_dir}/get_strand_amb_SNPs.R ${out_dir}/tmp_no_dupl.bim
plink --bfile ${out_dir}/tmp_gwas \
   --exclude ${out_dir}/tmp_strand_remove_snps.txt \
   --keep-allele-order \
   --make-bed --out ${out_dir}/tmp_gwas_no_AT_CG

#Perform pre-imputation QC - remove monomorphic SNPs, SNPs with high
#missingness, SNPs not in HWE
plink --bfile ${out_dir}/tmp_gwas_no_AT_CG \
   --maf 0.000001 --geno 0.05 --hwe 0.000001 \
   --keep-allele-order \
   --make-bed --out ${out_dir}/pre_qc

#Create vcf files for uploading to imputation server for QC
#Note that the encoding for chromosome is e.g. chr22, not 22
for ((chr=1; chr<=22; chr++)); do
    plink --bfile ${out_dir}/pre_qc \
       --chr $chr --keep-allele-order \
       --recode vcf --out ${out_dir}/tmp_chr${chr}
    vcf-sort ${out_dir}/tmp_chr${chr}.vcf | \
       sed -E 's/^([[:digit:]]+)/chr\1/' | \
       bgzip -c > ${out_dir}/chr${chr}_pre_qc.vcf.gz
done

#Report SNP counts
orig_snp_nr=`wc -l ${plink_prefix}.bim`
crossover_snp_nr=`wc -l ${out_dir}/tmp_gwas.bim`
nonamb_snp_nr=`wc -l ${out_dir}/tmp_gwas_no_AT_CG.bim`
qc_snp_nr=`wc -l ${out_dir}/pre_qc.bim`
echo "Original SNP nr: $orig_snp_nr"
echo "Crossovered SNP nr: $crossover_snp_nr"
echo "Non-ambiguous SNP nr: $nonamb_snp_nr"
echo "Final SNP nr after QC: $qc_snp_nr"

#Save report to text file
echo "Original SNP nr: $orig_snp_nr" > ${out_dir}/create_initial_input_log.txt
echo "Crossovered SNP nr: $crossover_snp_nr" >> ${out_dir}/create_initial_input_log.txt
echo "Non-ambiguous SNP nr: $nonamb_snp_nr" >> ${out_dir}/create_initial_input_log.txt
echo "Final SNP nr after QC: $qc_snp_nr" >> ${out_dir}/create_initial_input_log.txt

#Cleanup
rm ${out_dir}/tmp_*
