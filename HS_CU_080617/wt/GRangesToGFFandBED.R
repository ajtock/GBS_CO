#!/applications/R/R-3.3.2/bin/Rscript

# Convert GRanges objects to GFF files

# Usage:
# ./GRangesToGFFandBED.R GBS_crossovers

library(GenomicRanges)

args <- commandArgs(trailingOnly = T)
libName <- args[1]

load("/projects/ajt200/GBS_CO/HS_CU_080617/wt/COsGRcoords.RData")

COs <- data.frame(COsGRcoords)
COsGRcoords <- NULL

COsGFF <- cbind(COs[,1], rep(".", length(COs[,1])), rep(paste0(libName, "")), COs[,2], COs[,3], rep(".", length(COs[,1])), rep(".", length(COs[,1])), rep(".", length(COs[,1])), rep(".", length(COs[,1])))
colnames(COsGFF) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
write.table(COsGFF, file = paste0(libName, ".gff"), row.names = F, col.names = F, quote = F, sep = "\t")

COsBED <- data.frame(chr = as.character(COs$seqnames),
                     start = as.integer(COs$start-1),
                     end = as.integer(COs$end),
                     stringsAsFactors = F)
write.table(COsBED,
            file = paste0(libName, ".bed"),
            sep = "\t", quote = F,
            row.names = F, col.names = F)
