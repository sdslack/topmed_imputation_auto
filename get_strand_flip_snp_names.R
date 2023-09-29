#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

print("First trailing arg should be pre_qc directory path")
print("Second trailing arg should be post_qc directory path")
print("This _updated version fixes both strand flips and")
print("strand flip and allele switch. But, in the case of")
print("duplicate variants, or duplicate chr:pos locations")
print("where one has an allele mismatch, both variants are")
print("removed due to inability to identify which of the")
print("pair is correct.")

get_strand_flip_snp_names <- function(pre_qc_dir, post_qc_dir) {
  # Load snp-excluded.txt
  snp.excl <- read.delim(paste0(pre_qc_dir, "/snps-excluded.txt"),
                         stringsAsFactors = F)
  
  # Add column with chr:pos for all variants
  snp.excl$chr <- as.numeric(
    gsub("chr", "", unlist(strsplit(snp.excl$X.Position, split=":"))[seq(1,dim(snp.excl)[1]*4,4)]))
  snp.excl$pos <- as.numeric(
    unlist(strsplit(snp.excl$X.Position, split=":"))[seq(2,dim(snp.excl)[1]*4,4)])
  snp.excl$chr.pos <- paste0(snp.excl$chr, ":", snp.excl$pos)
  
  # Load bim so can get actual varID names
  bim <- read.table(paste0(pre_qc_dir, "/pre_qc.bim"),
                    stringsAsFactors = F)
  bim$chr.pos <- paste0(bim$V1, ":", bim$V4)
  
  # Add varID names to snp.excl by chr:pos
  merge <- merge(snp.excl, bim, by = c("chr.pos"))
  
  # Get strand flips & strand flip and allele switch
  snp.frame <- merge[grep("Strand flip", merge$FilterType),]
  
  # Remove strand flip & strand flip and allele switch from from snp.excl
  merge.to.excl <- merge[!merge$FilterType %in% snp.frame$FilterType,]
  
  # Write out SNPs to be excluded using PLINK --exclude
  merge.to.excl <- merge.to.excl[,c("V2")]
  write.table(merge.to.excl, paste0(post_qc_dir, "/tmp_snps-excl.txt"), 
              sep="\t", quote=F, row.names=F, col.names=F)
  
  # For both strand flip & strand flip and allele switch, flip strand
  # Only do if strand flips were found, otherwise create empty list
  if (nrow(snp.frame) > 0) {
    # Get varIDs from .bim by chr:pos for strand flip & strand flip and
    # allele switch
    snps <- bim$V2[bim$chr.pos %in% snp.frame$chr.pos]
    
    # Get ref alleles for strand flip and allele switch
    snp.frame.both <- snp.frame[grep("Strand flip and Allele switch", snp.frame$FilterType),]
    snp.frame.both$ref <- gsub("/[[:alpha:]]", "",
                               gsub("[[:alpha:]]+:", "", snp.frame.both$Info))
    
    # Get varID name with reference a2 allele for PLINK --a2-allele
    snp.frame.both <- snp.frame.both[,c("V2", "ref")]
    
  } else {
    snps <- character()
    snp.frame.both <- character()
    print("No strand flip SNPs were found in the excluded SNPs.")
  }
  
  # Write out SNPs for input into PLINK --flip, will flip both strand
  # flip & strand flip and allele switch
  write.table(snps, paste0(post_qc_dir, "/tmp_flip.txt"), 
              sep="\t", quote=F, row.names=F, col.names=F)
  
  # Write out SNPs for input into PLINK --a2-allele, will allele switch
  # strand flip and allele switch alleles after --flip
  write.table(snp.frame.both, paste0(post_qc_dir, "/tmp_a2-allele.txt"),
              sep="\t", quote = F, row.names = F, col.names = F)
}

get_strand_flip_snp_names(args[1], args[2])
