#!/applications/R/R-4.0.0/bin/Rscript

# Count crossover intervals per proportionally scaled window
# Argument 2 specifies the number of proportionally scaled windows
# to divide each chromosome arm into
# Apply this script to each population

# Usage from within directory containing crossover intervals as
# space-delimited text (*.txt) file:
# ./F2_CO_count_per_scaled_win_TelCen.R coller.wt 10

#popName <- "coller.wt"
#prop <- 10
#propName <- paste0(as.character(prop), "ths")

args <- commandArgs(trailingOnly = T)
popName <- args[1]
prop <- as.numeric(args[2])
propName <- paste0(as.character(prop), "ths")

library(GenomicRanges)
library(parallel)

outDir <- "genomeProfiles/"
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))

# Load COs and convert into GRanges object
COs <- read.csv(paste0("../COs/", popName, "cos.csv"),
                header = T)

popCOsGR_pooled <- GRanges(seqnames = paste0("Chr",
                                             COs$chr),
                           ranges = IRanges(start = COs$start,
                                            end = COs$end),
                           strand = "*",
                           midpoint = COs$cos,
                           lib = as.character(COs$lib))

libNames <- sort(unique(popCOsGR_pooled$lib))
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

# Create a list of crossover intervals, with each list element
# corresponding to crossover intervals for one F2 individual
indCOsGR_list <- mclapply(seq_along(libNames), function(x) {
  popCOsGR_pooled[popCOsGR_pooled$lib == libNames[x]]
}, mc.cores = detectCores())

# Count crossover intervals in proportionally scaled windows
# along left and right chromosome arms for each F2 individual
indCOs_perWindow_list <- mclapply(seq_along(indCOsGR_list), function(x) {
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
    # for each F2 individual
    COsChr <- indCOsGR_list[[x]][seqnames(indCOsGR_list[[x]]) == chrs[i]]
    COsChrWinL <- countOverlaps(query = windowsGRangesL,
                                subject = COsChr,
                                ignore.strand = TRUE)
    COsChrWinR <- countOverlaps(query = windowsGRangesR,
                                subject = COsChr,
                                ignore.strand = TRUE)
    # Create growing matrix in which each column represents a TEL-CEN
    # chromosome arm and each row represents a proportionally scaled window
    COsChrWin <- cbind(COsChrWinL, COsChrWinR)
    TelCenMatrix <- cbind(TelCenMatrix, COsChrWin)
  }
  TelCenMatrix
}, mc.cores = detectCores())
# Save list of matrices in which each list element is a matrix
# of crossover interval counts within scaled windows for one F2
save(indCOs_perWindow_list,
     file = paste0(outDir, popName,
                   "_pooled_F2_CO_frequency_",
                   propName, "_TelCenMatrix_list.RData"))
print(paste0(popName, "_pooled_F2_CO_frequency_", propName, "_TelCenMatrix_list.RData written to outDir"))

