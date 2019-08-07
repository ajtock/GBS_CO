#!/applications/R/R-3.5.0/bin/Rscript

# Profile mean coverage around COs and random loci

# Usage via Condor submission system on node7:
# csmit -m 20G -c 1 "/applications/R/R-3.5.0/bin/Rscript ./crossover_SNPfreq_profiles_commandArgs.R 5000 5kb 50 coller.filtarb collerF2.complete.tiger.txt"

library(GenomicRanges)
library(parallel)

#flankSize <- 5000
#flankName <- "5kb"
#winSize <- 50
#popName <- "coller.filtarb"
#SNPsFile <- "collerF2.complete.tiger.txt"

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
# Correct CO midpoints
COsGR$midpoint <- start(COsGR) +
                    (round(0.5 * ((end(COsGR) - start(COsGR)) + 1)))

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
		  strand = "*")

# For each chromosome, extend SNP coordinates to the
# maximum of crossover interval widths for that chromosome
extendSNPsGR <- GRanges()
for(i in 1:length(chrs)) {
  COsGRchr <- COsGR[seqnames(COsGR) == chrs[i]]
  SNPsGRchr <- SNPsGR[seqnames(SNPsGR) == chrs[i]]
  extendSNPsGRchr <- GRanges(seqnames = seqnames(SNPsGRchr),
                             ranges = IRanges(start = start(SNPsGRchr),
                                              end = start(SNPsGRchr) +
                                                       (max(width(COsGRchr)) - 1)),
#                                                      (round(mean(width(COsGRchr))) - 1) +
#                                                      (round(2*sd(width(COsGRchr))))),
                             strand = "*")
  extendSNPsGR <- append(extendSNPsGR, extendSNPsGRchr)
}

# Count number of SNPs that overlap extended SNP intervals
SNPoverlapsCount <- countOverlaps(query = extendSNPsGR,
                                  subject = SNPsGR,
                                  type = "any")
# Extract qualifying SNPs (those which are within the mean crossover
# interval width of at least one other SNP to their right)
ranLocSNPsGR <- SNPsGR[SNPoverlapsCount > 1]

# Remove SNPs too close to chromosome starts and ends
# [i.e., those within (flankSize - half of the min width of the
# crossover intervals) + (0.5*winSize) of chromosome starts,
# and those within (half of the max width of the crossover intervals)
# + (0.5*winSize) + (whichever is greater of flankSize OR
# half of the max width of the crossover intervals) of chromosome ends]
# (0.5*winSize) represents half of the midpoint-window within which,
# and around which, SNPs will be counted
ranLocSNPsGRrmDistal <- GRanges() 
for(i in 1:length(chrs)) {
  ranLocSNPsGRchr <- ranLocSNPsGR[seqnames(ranLocSNPsGR) == chrs[i]]
  COsGRchr <- COsGR[seqnames(COsGR) == chrs[i]]
  if(flankSize >= (0.5*(max(width(COsGRchr)) - 1))) {
    ranLocSNPsGRchr <- ranLocSNPsGRchr[( start(ranLocSNPsGRchr) >
                                           (flankSize - (0.5*(min(width(COsGRchr)) - 1))) +
                                           (0.5*winSize) ) &
                                       ( start(ranLocSNPsGRchr) <=
                                           chrLens[i] -
                                           (0.5*(max(width(COsGRchr)) - 1)) -
                                           (0.5*winSize) -
                                           flankSize )]
  } else {
    ranLocSNPsGRchr <- ranLocSNPsGRchr[( start(ranLocSNPsGRchr) >
                                           (flankSize - (0.5*(min(width(COsGRchr)) - 1))) +
                                           (0.5*winSize) ) &
                                       ( start(ranLocSNPsGRchr) <=
                                           chrLens[i] -
                                           (0.5*(max(width(COsGRchr)) - 1)) -
                                           (0.5*winSize) -
                                           (0.5*(max(width(COsGRchr)) - 1)) )]
  }                                        
  ranLocSNPsGRrmDistal <- append(ranLocSNPsGRrmDistal, ranLocSNPsGRchr)
}
ranLocSNPsGR <- ranLocSNPsGRrmDistal
rm(ranLocSNPsGRrmDistal)
gc()

# Define function to randomly select qualifying SNPs from ranLocSNPsGR
# as start coordinates and to define end coordinates based on the
# width distribution of COsGR
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = TRUE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Define seed so that random selections are reproducible
set.seed(93750174)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocGR contains the same number of loci per chromosome as COsGR
ranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  ranLocSNPsGRchr <- ranLocSNPsGR[seqnames(ranLocSNPsGR) == chrs[i]]
  COsGRchr <- COsGR[seqnames(COsGR) == chrs[i]]
  ranLocChrStart <- ranLocStartSelect(coordinates = start(ranLocSNPsGRchr),
                                      n = length(COsGRchr))
  ranLocChrIR <- IRanges(start = ranLocChrStart,
                         width = width(COsGRchr))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = ranLocChrIR,
                         strand = "*")
  ranLocGR <- append(ranLocGR, ranLocChrGR)
}
stopifnot(all.equal(width(ranLocGR), width(COsGR)))

# Define ranLocGR midpoints
ranLocGR$midpoint <- start(ranLocGR) +
                       (round(0.5 * ((end(ranLocGR) - start(ranLocGR)) + 1)))

# Create COsGRflank
COsGRflank <- COsGR
start(COsGRflank) <- COsGR$midpoint -
                       ((0.5*winSize)-1) - flankSize
end(COsGRflank) <- COsGR$midpoint +
                     (0.5*winSize) + flankSize
# Create ranLocGRflank
ranLocGRflank <- ranLocGR
start(ranLocGRflank) <- ranLocGR$midpoint -
                          ((0.5*winSize)-1) - flankSize
end(ranLocGRflank) <- ranLocGR$midpoint +
                        (0.5*winSize) + flankSize


## Count SNPs in winSize-bp windows along each locus
## in COsGRflank and ranLocGRflank

# Create GRangesList in which list element is a set of
# GRanges corresponding to winSize-bp windows along
# a given locus
# COs
COsGRflankWinGRL <- list()
for(x in seq_along(COsGRflank)) {
  COsGRflankWinStarts <- seq(from = start(COsGRflank[x]),
                             to = end(COsGRflank[x]),
                             by = winSize)
  COsGRflankWinGR <- GRanges(seqnames = seqnames(COsGRflank[x]),
                             ranges = IRanges(start = COsGRflankWinStarts,
                                              width = winSize),
                             strand = "*")
  stopifnot(start(COsGRflankWinGR)[1] == start(COsGRflank[x]) &
            end(COsGRflankWinGR)[length(COsGRflankWinGR)] == end(COsGRflank[x]) &
            range(width(COsGRflankWinGR)) == winSize &
            length(COsGRflankWinGR) == ((flankSize*2)+winSize)/winSize)
  COsGRflankWinGRL <- c(COsGRflankWinGRL, COsGRflankWinGR)
}

# Count SNPs in each window
COsSNPcountsList <- mclapply(seq_along(COsGRflankWinGRL), function(x) {
  countOverlaps(query = COsGRflankWinGRL[[x]],
                subject = SNPsGR,
                type = "any",
                ignore.strand = TRUE)
}, mc.cores = detectCores())
# Convert into dataframe
COsSNPcountsDF <- data.frame(matrix(data = unlist(COsSNPcountsList),
                                    nrow = length(COsSNPcountsList),
                                    byrow = TRUE))
colnames(COsSNPcountsDF) <- c(
  paste0("u",
         as.character(seq(from = 1,
                          to = ((((flankSize*2)+winSize)/winSize)-1)/2))),
  "t1",
  paste0("d",
         as.character(seq(from = 1,
                          to = ((((flankSize*2)+winSize)/winSize)-1)/2)))
)
write.table(COsSNPcountsDF,
            file = paste0(matDir, popName,
                          "_COs_SNP_frequency_feature_target_and_",
                          flankName, "_flank_dataframe.txt"),
            col.names = T, row.names = F, quote = F, sep = "\t")

# ranLoc
ranLocGRflankWinGRL <- list()
for(x in seq_along(ranLocGRflank)) {
  ranLocGRflankWinStarts <- seq(from = start(ranLocGRflank[x]),
                                to = end(ranLocGRflank[x]),
                                by = winSize)
  ranLocGRflankWinGR <- GRanges(seqnames = seqnames(ranLocGRflank[x]),
                                ranges = IRanges(start = ranLocGRflankWinStarts,
                                                 width = winSize),
                                strand = "*")
  stopifnot(start(ranLocGRflankWinGR)[1] == start(ranLocGRflank[x]) &
            end(ranLocGRflankWinGR)[length(ranLocGRflankWinGR)] == end(ranLocGRflank[x]) &
            range(width(ranLocGRflankWinGR)) == winSize &
            length(ranLocGRflankWinGR) == ((flankSize*2)+winSize)/winSize)
  ranLocGRflankWinGRL <- c(ranLocGRflankWinGRL, ranLocGRflankWinGR)
}

# Count SNPs in each window
ranLocSNPcountsList <- mclapply(seq_along(ranLocGRflankWinGRL), function(x) {
  countOverlaps(query = ranLocGRflankWinGRL[[x]],
                subject = SNPsGR,
                type = "any",
                ignore.strand = TRUE)
}, mc.cores = detectCores())
# Convert into dataframe
ranLocSNPcountsDF <- data.frame(matrix(data = unlist(ranLocSNPcountsList),
                                       nrow = length(ranLocSNPcountsList),
                                       byrow = TRUE))
colnames(ranLocSNPcountsDF) <- c(
  paste0("u",
         as.character(seq(from = 1,
                          to = ((((flankSize*2)+winSize)/winSize)-1)/2))),
  "t1",
  paste0("d",
         as.character(seq(from = 1,
                          to = ((((flankSize*2)+winSize)/winSize)-1)/2)))
)
write.table(ranLocSNPcountsDF,
            file = paste0(matDir, popName,
                          "_COs_SNP_frequency_ranLoc_target_and_",
                          flankName, "_flank_dataframe.txt"),
            col.names = T, row.names = F, quote = F, sep = "\t")


print(paste0(popName, " COs and ranLoc SNP frequency profile calculation complete"))

pdf("test.pdf")
plot(x = 1:201,
     y = colMeans(COsSNPcountsDF),
     ylim = c(min(colMeans(COsSNPcountsDF), colMeans(ranLocSNPcountsDF)),
              max(colMeans(COsSNPcountsDF), colMeans(ranLocSNPcountsDF))),
     type = "l",
     col = "red",
     lwd = 3)
lines(x = 1:201,
      y = colMeans(ranLocSNPcountsDF),
      type = "l",
      col = "grey40",
      lwd = 3)
dev.off()
