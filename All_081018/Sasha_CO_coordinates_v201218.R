#!/applications/R/R-3.4.0/bin/Rscript

# Convert GBS-derived haplotype files into crossover interval files

# Usage:
# ./Sasha_CO_coordinates_v201218.R "SB_wt_ColCvi_F2_2018_pooled" "v201218"

#prefix <- "SB_wt_ColCvi_F2_2018_pooled"
#suffix <- "v201218"

args <- commandArgs(trailingOnly = T)
prefix <- args[1]
suffix <- args[2]

library(GenomicRanges)
library(regioneR)

inDir <- "/projects/ajt200/GBS_CO/SB_060718/SB_msh2_wt/"
outDir <- "/projects/ajt200/GBS_CO/All_081018/coordinates/pooled_lanes/"
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))

# Convert dataframe into GRanges object
popCOs <- read.csv(paste0(inDir, prefix, "_",  suffix, ".csv"))
colnames(popCOs) <- c("X", "lane", "lib", "chr", "start", "end", "midpoint", "width")
popCOs$lib <- sprintf("%02d", popCOs$lib)
popCOsGR <- GRanges(seqnames = paste0("Chr", popCOs$chr),
                    ranges = IRanges(start = popCOs$start,
                                     end = popCOs$end),
                    strand = "*",
                    midpoint = popCOs$midpoint,
                    lib = paste0(popCOs$lane, ".", popCOs$lib))
popCOsGR <- sortSeqlevels(popCOsGR)
popCOsGR <- sort(popCOsGR, ignore.strand = T)
save(popCOsGR,
     file = paste0(outDir, prefix, "_GRanges.RData"))
maskGR <- toGRanges(data.frame(chr = c("Chr2", "Chr3"),
                               start = c(18286717, 1),
                               end = c(18957092, 13587785)))
maskCOsOverlaps <- findOverlaps(maskGR, popCOsGR,
                                ignore.strand = T,
                                select = "all")
print("Crossovers located within masked regions:")
print(popCOsGR[subjectHits(maskCOsOverlaps)])

