#!/applications/R/R-3.5.0/bin/Rscript

# Usage from within directory containing crossover csv files:
# ./1_pooled_lane_COs_to_GRanges_Emma.R wt.cos_bothlanes.csv EL_wt_ColBur_F2_2019 

args <- commandArgs(trailingOnly = TRUE)
popCSV <- args[1]
popName <- args[2]

library(GenomicRanges)

# Load crossover intervals
popCOs <- read.csv(popCSV)

# Convert to GRanges object
popCOsGR_pooled <- GRanges(seqnames = paste0("Chr", popCOs$chrs),
                           ranges   = IRanges(start = popCOs$start,
                                              end   = popCOs$stop),
                           strand   = "*",
                           midpoint = popCOs$cos,
                           lib      = paste0(popCOs$librarynum, ".", popCOs$lib))

# Sort crossovers GRanges object by seqnames (chromosome), start and end coordinates
popCOsGR_pooled <- sortSeqlevels(popCOsGR_pooled)
popCOsGR_pooled <- sort(popCOsGR_pooled, ignore.strand = T)

# Report crossover intervals
print(paste0(popName, " crossover intervals:"))
print(popCOsGR_pooled)

# Report number of individuals in population
print("Individuals with crossovers:")
print(length(unique(popCOsGR_pooled$lib)))

# Write crossovers intervals in GRanges format to current working directory
save(popCOsGR_pooled,
     file = paste0("./", popName, "_pooled_GRanges.RData"))
print(paste0(popName, "_pooled_GRanges.RData crossover intervals file written to current working directory"))
