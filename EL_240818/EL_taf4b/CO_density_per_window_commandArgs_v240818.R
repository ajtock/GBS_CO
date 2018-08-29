#!/applications/R/R-3.3.2/bin/Rscript

# Calculate crossover density per window tiled along each chromosome

# Usage:
# ./CO_density_per_window_commandArgs_v240818.R wt.cos_bind_chrs.txt wt.cos_bind_coords.txt wt 100000 100kb

library(segmentSeq)
library(doParallel)
registerDoParallel(cores = 5)
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

inDir <- "/projects/ajt200/GBS_CO/EL_240818/EL_wt/"
chrFile <- "wt.cos_bind_chrs.txt"
coordsFile <- "wt.cos_bind_coords.txt"
genotype <- "wt"
winSize <- 100000
winName <- "100kb"

args <- commandArgs(trailingOnly = TRUE)
inDir <- "./"
chrFile <- args[1]
coordsFile <- args[2]
genotype <- as.character(args[3])
winSize <- as.numeric(args[4])
winName <- as.character(args[5])

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

# Make cumulative genomes
sumchr <- cumsum(c(0, chrLens))
print(sumchr)
#[1] 0  30427671  50125960  73585790  92170846 119146348
sumchr_tot <- sumchr[length(sumchr)]
print(sumchr_tot)

centromeres <- sapply(seq_along(centromeres), function(x) {
  centromeres[x] + sumchr[x]
})
print(centromeres)
pericenStart <- sapply(seq_along(pericenStart), function(x) {
  pericenStart[x] + sumchr[x]
})
print(pericenStart)
pericenEnd <- sapply(seq_along(pericenEnd), function(x) {
  pericenEnd[x] + sumchr[x]
})
print(pericenEnd)

# Load CO coordinate files
CO_midpoints <- data.frame(cbind(read.table(paste0(inDir, chrFile))$V1,
                                 read.table(paste0(inDir, coordsFile))$V1))
# Convert into GRanges object
COsGR <- GRanges(seqnames = paste0("Chr", CO_midpoints[,1]),
                 ranges = IRanges(start = CO_midpoints[,2],
                                  end = CO_midpoints[,2]),
                 strand = "*")

cumCOsChrWinDat <- NULL
for(i in 1:5) {
  # Define windows as GRanges object
  print(i)
  seqWindows <- seq(1, chrLens[i], by = winSize)
  cumWindows <- seqWindows + sumchr[i]
  windowsIRanges <- IRanges(start = seqWindows,
                            width = winSize)
  windowsIRanges <- windowsIRanges[-length(windowsIRanges)]
  windowsIRanges <- append(windowsIRanges,
                           IRanges(start = seqWindows[length(seqWindows)],
                                   end = chrLens[i]))
  windowsGRanges <- GRanges(seqnames = chrs[i], strand = "*",
                            ranges = windowsIRanges)
  print(windowsGRanges)  

  # Define and count COs within windows
  COsChr <- COsGR[seqnames(COsGR) == chrs[i]]
  COsChrWin <- countOverlaps(windowsGRanges, COsChr)
  COsChrWinDat <- cbind(cumWindows, COsChrWin)
  write.table(COsChrWinDat,
              file = paste0(inDir, genotype, "_CO_density_chr", i, "_", winName, ".txt"))
  cumCOsChrWinDat <- rbind(cumCOsChrWinDat, COsChrWinDat)
}
write.table(cumCOsChrWinDat,
            file = paste0(inDir, genotype, "_CO_density_genome_", winName, ".txt"))

