#!/applications/R/R-3.5.0/bin/Rscript

# Count crossover intervals per proportionally scaled window
# Argument 2 specifies the number of proportionally scaled windows
# to divide each chromosome arm into

# Usage from within directory containing crossover intervals as GRanges (*.RData) file:
# ./2b_F2_CO_count_per_win_Emma.R EL_wt_ColBur_F2_2019 1000000 1Mb

args <- commandArgs(trailingOnly = T)
popName <- args[1]
winSize <- as.numeric(args[2])
winName <- args[3]

library(GenomicRanges)
library(parallel)

print(popName)
load(paste0("./", popName, "_pooled_GRanges.RData"))
if( exists("popCOsGR") == TRUE ) {
  popCOsGR_pooled <- popCOsGR
}
# Convert library (F2 individual) names (e.g., from "1.1" to "1.01") to enable
# sorting individuals by increasing number
popCOsGR_pooled$lib[grep("\\d\\.\\d$",
                         popCOsGR_pooled$lib)] <- paste0(substr(x = popCOsGR_pooled$lib[grep("\\d\\.\\d$",
                                                                                        popCOsGR_pooled$lib)],
                                                                start = 1, stop = 1),
                                                         ".0",
                                                         substr(x = popCOsGR_pooled$lib[grep("\\d\\.\\d$",
                                                                                             popCOsGR_pooled$lib)],
                                                                start = 3, stop = 3))
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

# Genomic definitions with cumulative coordinates
sumchr <- cumsum(c(0, chrLens))
print(sumchr)
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

# Create a list of crossover intervals, with each list element
# corresponding to crossover intervals for one F2 individual
indCOsGR_list <- mclapply(seq_along(libNames), function(x) {
  popCOsGR_pooled[popCOsGR_pooled$lib == libNames[x]]
}, mc.cores = detectCores())

# Count crossover intervals in winName windows along each chromosome
# for each F2 individual
indCOs_perWindow_list <- mclapply(seq_along(indCOsGR_list), function(x) {
  cumWinCOsFreq <- NULL
  for(i in seq_along(chrs)) {
    # Define adjacent windows
    seqWindows <- seq(1, chrLens[i], by = winSize)
    cumWindows <- seqWindows + sumchr[i]
    windowsIRanges <- IRanges(start = seqWindows,
                              width = winSize)
    windowsIRanges <- windowsIRanges[-length(windowsIRanges)]
    windowsIRanges <- append(windowsIRanges,
                             IRanges(start = seqWindows[length(seqWindows)],
                                     end = chrLens[i]))
    windowsGRanges <- GRanges(seqnames = chrs[i],
                              ranges = windowsIRanges,
                              strand = "*")
    print(windowsGRanges)
    
    # Count crossover intervals
    chrCOs <- indCOsGR_list[[x]][seqnames(indCOsGR_list[[x]]) == chrs[i]]
    winCOs <- countOverlaps(query = windowsGRanges,
                            subject = chrCOs,
                            ignore.strand = T)
    COsDF <- data.frame(chr = seqnames(windowsGRanges),
                        chrWindow = as.integer(seqWindows),
                        cumWindow = as.integer(cumWindows),
                        COs = as.numeric(winCOs))
    cumWinCOsFreq <- rbind(cumWinCOsFreq, COsDF)
  }
  cumWinCOsFreq
}, mc.cores = detectCores())
save(indCOs_perWindow_list,
     file = paste0("./", popName,
                   "_pooled_F2_CO_frequency_per_",
                   winName, "_window_list.RData"))
print(paste0(popName, "_pooled_F2_CO_frequency_per_", winName, "_window_list.RData written to current working directory"))
