#!/bin/bash

#Set arguments
if [ "$#" -eq  "0" ]
then
   echo "Usage: ${0##*/} <chr> <rsq> <maf> <dir>"
   echo "Script filters new TOPMed INFO files formatted as VCFs."
   echo "Keeps TYPED or IMPUTED with Rsq less than given threshold,"
   echo "and filters MAF to given threshold."
   echo "Script expects file formats: chr#.dose.vcf.gz & chr#.info.gz"
   exit
fi

chr=$1
rsq=$2
maf=$3
dir=$4

# Set filter
to_filt="((INFO/TYPED = 1 | (INFO/IMPUTED = 1 & INFO/R2 > ${rsq})) & INFO/MAF > ${maf})"

# Filter INFO file, for smaller output with kept variables
    # ((typed OR (imputed & R2)) & MAF)
bcftools filter -i \
    "$to_filt" \
    "${dir}/chr${chr}.info.gz" -o "${dir}/${chr}_clean.info"

# Filter actual VCF
    # To note: can't just do IDs because some appear more than once above & below
    # filtering criteria
bcftools filter -i \
    "$to_filt" \
    "${dir}/chr${chr}.dose.vcf.gz" -o "${dir}/${chr}_clean.vcf.gz"
tabix "${dir}/${chr}_clean.vcf.gz"
