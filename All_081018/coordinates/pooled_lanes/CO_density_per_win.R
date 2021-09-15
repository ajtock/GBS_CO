#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./CO_density_per_win.R "HS_wt_ColLer_F2_2018" 500000 500kb

args <- commandArgs(trailingOnly = T)
popName <- args[1]
winSize <- as.numeric(args[2])
winName <- as.character(args[3])

library(GenomicRanges)

inDir <- "/projects/ajt200/GBS_CO/All_081018/coordinates/pooled_lanes/"
genomeProfilesDir <- paste0(inDir, "genomeProfiles/")
system(paste0("[ -d ", genomeProfilesDir, " ] || mkdir ", genomeProfilesDir))

print(popName)
load(paste0(inDir, popName, "_pooled_GRanges.RData"))
if( exists("popCOsGR") == TRUE ) {
popCOsGR_pooled <- popCOsGR
}
print("Individuals with crossovers:")
print(length(unique(popCOsGR_pooled$lib)))

# Calculate crossover density per window
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
  COsChr <- popCOsGR_pooled[seqnames(popCOsGR_pooled) == chrs[i]]
  COsChrWin <- countOverlaps(windowsGRanges, COsChr)
  COsChrWinDat <- cbind(cumWindows, COsChrWin)
  write.table(COsChrWinDat,
              file = paste0(genomeProfilesDir, popName, "_pooled_CO_density_chr", i, "_", winName, ".txt"))
  cumCOsChrWinDat <- rbind(cumCOsChrWinDat, COsChrWinDat)
}
write.table(cumCOsChrWinDat,
            file = paste0(genomeProfilesDir, popName, "_pooled_CO_density_genome_", winName, ".txt"))
