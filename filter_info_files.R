#!/usr/bin/env Rscript

# To note: This script was modified from the "New Info Files Filter" tool on
# Seven Bridges, which was edited by Erika and Randi from a script by Meher.

args <- commandArgs(trailingOnly = TRUE)
info_file <- args[1]
rsq <- args[2]
maf <- args[3]

print("First trailing arg should be input INFO file, second should be Rsq")
print("filter, third should be MAF filter.")
print("Assumes INFO file has columns 'Rsq', 'MAF', & 'Genotyped'.")

info_file_name <- basename(info_file)
info_file_dir <- dirname(info_file)

snp_out_file <- paste0(info_file_dir, "/", info_file_name,
                       "_maf", maf, "_rsq", rsq, "_snps.txt")
info_out_file <- paste0(info_file_dir, "/clean_", info_file_name)

info <- read.delim(gzfile(info_file), stringsAsFactors=F)

info$Rsq_num <- as.numeric(info$Rsq)
print("Check that NA Rsq values are for Typed Only variants:")
all(info$Genotyped[is.na(info$Rsq_num)] == "Typed_Only")

# Select rows from info that are either genotyped or above the Rsq threshold
# and also are above the MAF threshold
info_clean <- info[
  which(((info$Genotyped == "Typed_Only") |
           (!is.na(info$Rsq_num) & (info$Rsq_num >= rsq)))
        & (info$MAF >= maf)),]
extracted_rows <- (
  (info$Genotyped == "Typed_Only") |
    (!is.na(info$Rsq_num) & (info$Rsq_num >= rsq))) &
  (info$MAF >= maf)
snps <- info$SNP[extracted_rows]

# Write out cleaned files
write.table(info_clean, info_out_file, sep="\t", quote=F, row.names=F)
write.table(snps, snp_out_file,  sep="\t", quote=F, row.names=F, col.names=F)
