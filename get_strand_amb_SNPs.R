#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

print("Script is called by create_initial_input.sh. First trailing")
print("arg should be file path to .bim after crossover (if done) or")
print("after removal of SNPs with duplicate positions. Writes new")
print("file tmp_strand_remove_snps.txt.")

get_strand_amb_SNPs <- function(bim_file) {
  bim <- read.table(bim_file, stringsAsFactors=F)
  snps <- bim$V2[((bim$V5 == "A") & (bim$V6 == "T")) |
                   ((bim$V5 == "T") & (bim$V6 == "A")) |
                   ((bim$V5 == "C") & (bim$V6 == "G")) |
                   ((bim$V5 == "G") & (bim$V6 == "C"))]
  path <- dirname(bim_file)
  write.table(x = snps, file = paste0(path, "/tmp_strand_remove_snps.txt"),
              sep="\t", quote=F, row.names=F, col.names=F)
}

get_strand_amb_SNPs(args[1])
