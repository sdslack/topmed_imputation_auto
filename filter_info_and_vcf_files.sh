#!/bin/bash

#Set arguments
if [ "$#" -eq  "0" ]
then
   echo "Usage: ${0##*/} <chr> <rsq> <maf> <in_dir> <out_dir>"
   echo "Script filters new TOPMed INFO files formatted as VCFs."
   echo "Keeps TYPED or IMPUTED with Rsq less than given threshold,"
   echo "and filters MAF to given threshold."
   echo "Script expects file formats: chr#.dose.vcf.gz & chr#.info.gz"
   exit
fi

chr=$1
rsq=$2
maf=$3
in_dir=$4
out_dir=$5

# Set filter
to_filt="((INFO/TYPED = 1 | (INFO/IMPUTED = 1 & INFO/R2 > ${rsq})) & INFO/MAF > ${maf})"

# Filter INFO file, for smaller output with kept variables
    # ((typed OR (imputed & R2)) & MAF)
bcftools filter -i \
    "$to_filt" \
    "${in_dir}/chr${chr}.info.gz" -o "${out_dir}/chr${chr}_clean.info"

# Write out list of RSIDs want to keep
bcftools query -f '%ID\n' \
    "${out_dir}/chr${chr}_clean.info" > "${out_dir}/chr${chr}_maf${maf}_rsq${rsq}_snps.txt"

# Filter VCF to these IDs using PLINK
plink2 --vcf "${in_dir}/chr${chr}.dose.vcf.gz" \
  --export vcf 'bgz' \
  --extract "${out_dir}/chr${chr}_maf${maf}_rsq${rsq}_snps.txt" \
  --out "${out_dir}/tmp_chr${chr}_clean"

# Finally, clean up RSIDs that may have appear more than once
bcftools filter -i \
    "$to_filt" \
    "${out_dir}/tmp_chr${chr}_clean.vcf.gz" -o "${out_dir}/chr${chr}_clean.vcf.gz"
tabix "${out_dir}/chr${chr}_clean.vcf.gz"

# Clean up
rm ${out_dir}/tmp_*