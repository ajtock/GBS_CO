#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./CO_count_per_scaled_win_TelCen.R "HS_wt_ColLer_F2_2018" 10 10ths

args <- commandArgs(trailingOnly = T)
popName <- args[1]
prop <- as.numeric(args[2])
propName <- args[3]

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

for(i in 1:5) {
  # Define windows as GRanges object
  print(i)
  seqWindowStartsL <- as.integer(seq(from = 1, to = centromeres[i], by = centromeres[i]/prop))
  seqWindowEndsL <- c(seqWindowStartsL[2:length(seqWindowStartsL)]-1,
                      centromeres[i]-1)
  windowsGRangesL <- GRanges(seqnames = chrs[i],
                             ranges = IRanges(start = seqWindowStartsL,
                                              end = seqWindowEndsL),
                             strand = "*")
  print(windowsGRangesL)
  seqWindowStartsR <- as.integer(seq(from = centromeres[i], to = chrLens[i], by = (chrLens[i]-centromeres[i]+1)/prop))
  seqWindowEndsR <- c(seqWindowStartsR[2:length(seqWindowStartsR)]-1,
                      chrLens[i])
  windowsGRangesR <- rev(GRanges(seqnames = chrs[i],
                                 ranges = IRanges(start = seqWindowStartsR,
                                                 end = seqWindowEndsR),
                                 strand = "*"))
  print(windowsGRangesR)

  # Define and count COs within windows
  COsChr <- popCOsGR_pooled[seqnames(popCOsGR_pooled) == chrs[i]]
  COsChrWinL <- countOverlaps(windowsGRangesL, COsChr)
  write.table(COsChrWinL,
              file = paste0(genomeProfilesDir, popName, "_pooled_CO_density_chr", i, "_leftArm_", propName, ".txt"))
  COsChrWinR <- countOverlaps(windowsGRangesR, COsChr)
  write.table(COsChrWinR,
              file = paste0(genomeProfilesDir, popName, "_pooled_CO_density_chr", i, "_rightArm_", propName, ".txt"))
}
