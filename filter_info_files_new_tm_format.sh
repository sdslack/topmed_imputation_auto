#!/bin/bash

#Set arguments
if [ "$#" -eq  "0" ]
then
   echo "Usage: ${0##*/} <info.gz> <rsq> <maf>"
   echo "Script filters new TOPMed INFO files formatted as VCFs."
   echo "Keeps TYPED or IMPUTED with Rsq less than given threshold,"
   echo "and filters MAF to given threshold."
   exit
fi

info=$1
rsq=$2
maf=$3

info_name=$(basename "$info")
echo $info_name
info_path=$(dirname "$info")

# Filter INFO file
bcftools filter -i \
    "((INFO/TYPED = 1 | (INFO/IMPUTED = 1 & INFO/R2 > ${rsq})) & INFO/MAF > ${maf})" \
    "$info" -o "${info_path}/clean_${info_name}"

# Write out list of IDs want to keep
bcftools query -f '%ID\n' \
    "${info_path}/clean_${info_name}" > "${info}_maf${maf}_rsq${rsq}_snps.txt"
