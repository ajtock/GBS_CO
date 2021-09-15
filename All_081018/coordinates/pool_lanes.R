#!/applications/R/R-3.3.2/bin/Rscript

# Usage:
# ./pool_lanes.R "_wt_ColLer_F2_2016" 

args <- commandArgs(trailingOnly = T)
#popName <- "_wt_ColLer_F2_2016"
popName <- args[1]

library(GenomicRanges)

inDir <- "/projects/ajt200/GBS_CO/All_081018/coordinates/"
outDir <- paste0(inDir, "pooled_lanes/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))

# Pool lanes
filePath <- system(paste0("ls ", inDir, "*", popName, "_lane*_GRanges.RData"),
                   intern = T)

popCOsGR_pooled <- GRanges()
for(i in 1:length(filePath)) {
  load(filePath[i])
  popCOsGR$lib <- paste0(i, ".", popCOsGR$lib)
  popCOsGR_pooled <- append(popCOsGR_pooled, popCOsGR)
}
popCOsGR_pooled <- sortSeqlevels(popCOsGR_pooled)
popCOsGR_pooled <- sort(popCOsGR_pooled, ignore.strand = T)
print(popName)
print("Individuals with crossovers:")
print(length(unique(popCOsGR_pooled$lib)))
save(popCOsGR_pooled,
     file = paste0(outDir, popName, "_pooled_GRanges.RData"))

