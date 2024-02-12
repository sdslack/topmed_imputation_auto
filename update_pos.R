#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

print("Script is called by create_initial_input.sh. First trailing ")
print("arg should be file path to tmp_out.bed created by CrossMap. ")
print("Second trailing arg should be file path to tmp_gwas.bim with ")
print("only SNPs successfully crossovered. Writes new tmp_gwas.bim, ")
print("or equivalent file name that was passed in.")

update_pos <- function(in_bed, in_bim) {
  in.bed <- read.table(in_bed)[,c(2,4)]
  names(in.bed)[1] <- "NEW.POS"
  in.bim <- read.table(in_bim)[,c(1,2,5,6)]
  in.bim$ORDER <- seq(1, length(in.bim$V1))
  merged <- merge(in.bim, in.bed, by.x="V2", by.y="V4")
  merged <- merged[order(merged$ORDER),]
  out.bim <- data.frame(V1=merged$V1,
                        V2=merged$V2,
                        V3=rep(0,length(in.bim$V1)),
                        V4=merged$NEW.POS,
                        V5=merged$V5,
                        V6=merged$V6)
  write.table(out.bim, in_bim,
              sep="\t", quote=F, row.names=F, col.names=F)
}

update_pos(args[1], args[2])
