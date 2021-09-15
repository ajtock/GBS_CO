#!/applications/R/R-3.3.2/bin/Rscript

# Convert GRanges objects into GFF3 files

# Usage:
# ./GRangesTOgff.R CU_cmt3_ColLer_F2_2016

args <- commandArgs(trailingOnly = T)
popName <- args[1]

library(GenomicRanges)

print(popName)
load(paste0("./", popName, "_pooled_GRanges.RData"))
if( exists("popCOsGR") == TRUE ) {
  COsGR <- popCOsGR
}
if( exists("popCOsGR_pooled") == TRUE ) {
  COsGR <- popCOsGR_pooled
}

# Convert into GFF3 format
COs <- data.frame(COsGR)
COsGFF <- data.frame(chr = as.character(COs$seqnames),
                     source = as.character(rep(".")),
                     feature = as.character(rep(paste0(popName, "_CO"))),
                     start = as.integer(COs$start),
                     end = as.integer(COs$end),
                     score = as.character(rep(".")),
                     strand = as.character(rep(".")),
                     phase = as.character(rep(".")),
                     midpoint = as.integer(COs$midpoint))
write.table(COsGFF,
            file = paste0(popName,
                          "_pooled.gff"),
            row.names = F, col.names = F, quote = F, sep = "\t")
