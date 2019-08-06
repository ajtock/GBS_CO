#!/applications/R/R-3.4.0/bin/Rscript

# Profile mean coverage around COs and random loci

# Usage via Condor submission system on node7:
# csmit -m 20G -c 1 "/applications/R/R-3.4.0/bin/Rscript 5000 5kb 50 ./crossover_SNPfreq_profiles_commandArgs.R coller.filtarb collerF2.complete.tiger.txt"

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc_R3.4.0.R")

library(EnrichedHeatmap)

flankSize <- 5000
flankName <- "5kb"
winSize <- 50
popName <- "coller.filtarb"
SNPsFile <- "collerF2.complete.tiger.txt"

args <- commandArgs(trailingOnly = T)
flankSize <- as.numeric(args[1])
flankName <- as.character(args[2])
winSize <- as.numeric(args[3])
popName <- args[4]
SNPsFile <- args[5]

matDir <- paste0("matrices/")
system(paste0("[ -d ", matDir, " ] || mkdir ", matDir))

# Chromosome definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrStart <- c(rep(1, 5))
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

# Load COs and convert into GRanges object
COs <- read.table(paste0("../../COs/", popName, "cos.txt"),
                  header = T)
COsGR <- GRanges(seqnames = paste0("Chr",
                                   COs$chr),
                 ranges = IRanges(start = COs$start,
                                  end = COs$end),
                 strand = "*",
                 midpoint = COs$cos,
                 lib = paste0(as.character(COs$lane), ".",
                              as.character(COs$lib)))

# Convert library (individual) names (e.g., from "1.1" to "1.01") to enable
# sorting individuals by increasing number
COsGR$lib[grep("\\d\\.\\d$",
               COsGR$lib)] <- paste0(substr(x = COsGR$lib[grep("\\d\\.\\d$",
                                                               COsGR$lib)],
                                            start = 1, stop = 1),
                                     ".0",
                                     substr(x = COsGR$lib[grep("\\d\\.\\d$",
                                                               COsGR$lib)],
                                            start = 3, stop = 3))
libNames <- sort(unique(COsGR$lib))
print("Individuals with crossovers:")
print(length(unique(COsGR$lib)))


# Load SNPs
SNPs <- read.table(paste0("../../SNPs/", SNPsFile),
		   header = T)
SNPsGR <- GRanges(seqnames = paste0("Chr",
				    SNPs$snps.chrs),
		  ranges = IRanges(start = SNPs$X1,
				   end = SNPs$X1),
		  strand = "*",
		  coverage = rep(1, dim(SNPs)[1]))

# For each chromosome, make GRanges object consisting of all SNP-to-SNP ranges
# whose width is <= mean+1 SD of crossover interval widths for that chromosome
interSNPsGR <- GRanges()
for(i in 1:length(chrs)) {
  COsGRchr <- COsGR[seqnames(COsGR) == chrs[i]]
  SNPsGRchr <- SNPsGR[seqnames(SNPsGR) == chrs[i]]
  interSNPsGRchr <- GRanges()
  for(x in (1:(length(SNPsGRchr)-1))) {
    SNPxGRends <- SNPsGRchr[( end(SNPsGRchr) > start(SNPsGRchr[x]) ) &
                            ( (end(SNPsGRchr) - start(SNPsGRchr[x]))+1 ) <=
                            ( mean(width(COsGRchr)) + (1*sd(width(COsGRchr))) )]
    if(length(SNPxGRends) > 0) { 
      SNPxGR <- GRanges(seqnames = seqnames(SNPxGRends),
                        ranges = IRanges(start = start(SNPsGRchr[x]),
                                         end = end(SNPxGRends)),
                        strand = "*")
      interSNPsGRchr <- append(interSNPsGRchr, SNPxGR)
    }
  }
  interSNPsGR <- append(interSNPsGR, interSNPsGRchr)
}


SNPsGRextend <- GRanges(seqnames = seqnames(SNPsGR),
                        ranges = IRanges(start = start(SNPsGR)-mean(width(COsGR)),
                                         end = end(SNPsGR)+mean(width(COsGR))),
                        strand = "*",
                        coverage = rep(1, dim(SNPs)[1]))


# Define function to select randomly positioned loci of the same
# width distribution as COsGR
ranLocStartSelect <- function(n, minStart, maxEnd) {
  ceiling(runif(n = n,
                min = minStart-1,
                max = maxEnd))
  #sample(x = minStart:maxEnd,
  #       size = n,
  #       replace = TRUE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Define seed so that random selections are reproducible
set.seed(93750174)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocGR contains the same number of loci per chromosome as COsGR
chrs <- chrs[grep(genomeName, chrs)]
ranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  COsChrGR <- COsGR[seqnames(COsGR) == chrs[i]]
  regionChrGR <- regionGR[seqnames(regionGR) == chrs[i]]
  # Contract regionChrGR so that random loci do not overlap masked region
  # and do not extend beyond chromosome ends
  end(regionChrGR) <- end(regionChrGR)-max(width(COsChrGR))
  ranLocChrStart <- ranLocStartSelect(coordinates = unlist(as.vector(ranges(regionChrGR))),
                                      n = length(COsChrGR))
  ranLocChrIR <- IRanges(start = ranLocChrStart,
                         width = width(COsChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = ranLocChrIR,
                         strand = "*")
  ranLocGR <- append(ranLocGR, ranLocChrGR)
}


# Define matrix and column mean SNP frequency outfile
outDF <- list(paste0(matDir, popName, "_COs_",
                     "SNP_frequency_feature_smoothed_target_and_",
                     flankName, "_flank_dataframe.txt"),
              paste0(matDir, popName, "_COs_in_",
                     "SNP_frequency_ranLoc_smoothed_target_and_",
                     flankName, "_flank_dataframe.txt"))
outDFcolMeans <- list(paste0(matDir, popName, "_COs_in_",
                             "SNP_frequency_feature_smoothed_target_and_",
                             flankName, "_flank_dataframe_colMeans.txt"),
                      paste0(matDir, popName, "_COs_in_",
                             "SNP_frequency_ranLoc_smoothed_target_and_",
                             flankName, "_flank_dataframe_colMeans.txt"))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
# containing normalised coverage values around target and random loci
covMatrix(signal = SNPsGR,
          feature = COsGR,
          ranLoc = ranLocGR,
          featureSize = mean(width(COsGR)),
          flankSize = flankSize,
          winSize = winSize,
          outDF = outDF,
          outDFcolMeans = outDFcolMeans)
print(paste0(genomeName, "-genome ", region, " ", popName, " COs SNP frequency profile calculation complete"))
