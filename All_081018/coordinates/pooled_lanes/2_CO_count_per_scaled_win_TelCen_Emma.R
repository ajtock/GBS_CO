#!/applications/R/R-3.5.0/bin/Rscript

# Count crossover intervals per proportionally scaled window
# Argument 2 specifies the number of proportionally scaled windows
# to divide each chromosome arm into

# Usage from within directory containing crossover intervals as GRanges (*.RData) file:
# ./2_CO_count_per_scaled_win_TelCen_Emma.R EL_wt_ColBur_F2_2019 10

args <- commandArgs(trailingOnly = T)
popName <- args[1]
prop <- as.numeric(args[2])
propName <- paste0(as.character(prop), "ths")

library(GenomicRanges)

print(popName)
load(paste0("./", popName, "_pooled_GRanges.RData"))
if( exists("popCOsGR") == TRUE ) {
  popCOsGR_pooled <- popCOsGR
}
print("Individuals with crossovers:")
print(length(unique(popCOsGR_pooled$lib)))

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

# Count crossover intervals in proportionally scaled windows
# along left and right chromosome arms
TelCenMatrix <- NULL
for(i in seq_along(chrs)) {
  print(chrs[i])
  # Define GRanges objects containing proportionally scaled windows
  # left and right of centromere
  # Left arm
  seqWindowStartsL <- as.integer(seq(from = 1,
                                     to = centromeres[i],
                                     by = centromeres[i]/prop))
  # Check that last window start coordinate is < centromeres[i]-1
  stopifnot(seqWindowStartsL[length(seqWindowStartsL)] < centromeres[i]-1)
  seqWindowEndsL <- c(seqWindowStartsL[2:length(seqWindowStartsL)]-1,
                      centromeres[i]-1)
  windowsGRangesL <- GRanges(seqnames = chrs[i],
                             ranges = IRanges(start = seqWindowStartsL,
                                              end = seqWindowEndsL),
                             strand = "*")
  print(windowsGRangesL)
  # Check that last window end coordinate == centromere[i]-1
  # (i.e., most proximal coordinate)
  stopifnot(end(windowsGRangesL[length(windowsGRangesL)]) == centromeres[i]-1)

  # Right arm
  seqWindowStartsR <- as.integer(seq(from = centromeres[i],
                                     to = chrLens[i],
                                     by = (chrLens[i]-centromeres[i]+1)/prop))
  # Check that last window start coordinate is < chrLens[i]
  stopifnot(seqWindowStartsR[length(seqWindowStartsR)] < chrLens[i])
  seqWindowEndsR <- c(seqWindowStartsR[2:length(seqWindowStartsR)]-1,
                      chrLens[i])
  # Reverse the order of right-arm windows so that
  # distal windows precede proximal windows
  windowsGRangesR <- rev(GRanges(seqnames = chrs[i],
                                 ranges = IRanges(start = seqWindowStartsR,
                                                 end = seqWindowEndsR),
                                 strand = "*"))
  print(windowsGRangesR)
  # Check that first window end coordinate == chrLens[i]
  # (i.e., most distal coordinate)
  stopifnot(end(windowsGRangesR[1]) == chrLens[i])

  # Count crossover intervals within windows
  # and write to current working directory
  COsChr <- popCOsGR_pooled[seqnames(popCOsGR_pooled) == chrs[i]]
  COsChrWinL <- countOverlaps(windowsGRangesL, COsChr)
  COsChrWinR <- countOverlaps(windowsGRangesR, COsChr)
  COsChrWin <- cbind(COsChrWinL, COsChrWinR)
  TelCenMatrix <- cbind(TelCenMatrix, COsChrWin)
}
write.table(TelCenMatrix,
            file = paste0("./", popName,
                          "_pooled_CO_frequency_",
                          propName, "_TelCenMatrix.txt"))
print(paste0(popName, "_pooled_CO_frequency_", propName, "_TelCenMatrix.txt written to current working directory"))

