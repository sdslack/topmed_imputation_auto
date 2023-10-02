# Edited Script sent by Meher


# Read in tool inputs
source("cwl_inputs.R")


args <- commandArgs(trailingOnly = TRUE)

# The directory where the files will go
snp.out.file <- paste0(file_name, "_maf",maf, "_rsq", rsq, "_snps.txt")
info.out.file <- paste0("clean_", file_name )



info <- read.delim(gzfile(info.file), stringsAsFactors=F)

info$Rsq_num <- as.numeric(info$Rsq)
info_clean <- info[which(((info$Genotyped == "Typed_Only") | (!is.na(info$Rsq_num) & (info$Rsq_num >= rsq)) ) & (info$MAF >= maf)),]
extracted.rows <- ((info$Genotyped == "Typed_Only") | (!is.na(info$Rsq_num) & (info$Rsq_num >= rsq)) ) & (info$MAF >= maf )
snps <- info$SNP[extracted.rows]


# Write out cleaned files
write.table(info_clean, info.out.file, sep="\t", quote=F, row.names=F)
write.table(snps, snp.out.file,  sep="\t", quote=F, row.names=F, col.names=F)

