#!/applications/R/R-3.5.0/bin/Rscript

# Calculate windowed SNP frequencies around the
# midpoints of CO intervals and random loci

# Usage via Condor submission system on node7:
# csmit -m 20G -c 24 "/applications/R/R-3.5.0/bin/Rscript ./crossover_SNPfreq_profiles_commandArgs_mid2SNPs_minWidth.R 15000 15kb 100 100bp coller.filtarb collerF2.complete.tiger.txt"

#flankSize <- 15000
#flankName <- "15kb"
#winSize <- 100
#winName <- "100bp"
#popName <- "coller.filtarb"
#SNPsFile <- "collerF2.complete.tiger.txt"

args <- commandArgs(trailingOnly = T)
flankSize <- as.numeric(args[1])
flankName <- as.character(args[2])
winSize <- as.numeric(args[3])
winName <- as.character(args[4])
popName <- args[5]
SNPsFile <- args[6]

library(GenomicRanges)
library(parallel)

matDir <- paste0("matrices/")
plotDir <- paste0("plots/")
system(paste0("[ -d ", matDir, " ] || mkdir ", matDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Chromosome definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrStart <- c(rep(1, 5))
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
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
                    (round(0.5 * (end(COsGR) - start(COsGR))))

# Convert library (individual) names (e.g., from "1.1" to "1.01") to enable
# sorting individuals by increasing number
COsGR$lib[grep("\\d+\\.\\d$",
               COsGR$lib)] <- paste0(sub(pattern = "(\\d+)\\.\\d",
                                         replacement = "\\1",
                                         x = COsGR$lib[grep(pattern = "\\d+\\.\\d$",
                                                            x = COsGR$lib)]),
                                     ".0",
                                     sub(pattern = "\\d+\\.(\\d)",
                                         replacement = "\\1",
                                         x = COsGR$lib[grep(pattern = "\\d+\\.\\d$",
                                                            x = COsGR$lib)]))
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
SNPsGR <- sortSeqlevels(SNPsGR)
SNPsGR <- sort(SNPsGR)

# For each chromosome, create sequential SNP intervals
# and remove those with widths < the 50th percentile (median) width and
# > the max width of crossover intervals for that chromosome
interSNPsGR <- GRanges()
for(i in 1:length(chrs)) {
  COsGRchr <- COsGR[seqnames(COsGR) == chrs[i]]
  SNPsGRchr <- SNPsGR[seqnames(SNPsGR) == chrs[i]]
  interSNPsGRchr <- GRanges(seqnames = chrs[i],
                            ranges = IRanges(start = start(SNPsGRchr[1:(length(SNPsGRchr)-1)]),
                                             end = start(SNPsGRchr[2:(length(SNPsGRchr))])),
                            strand = "*")
  interSNPsGRchr <- interSNPsGRchr[width(interSNPsGRchr) >= quantile(width(COsGRchr), 0.50)[[1]] &
                                   width(interSNPsGRchr) <= max(width(COsGRchr))]
  interSNPsGRchr$midpoint <- start(interSNPsGRchr) +
                               (round(0.5 * ((end(interSNPsGRchr) - start(interSNPsGRchr)) + 1)))
  # Remove SNP intervals too close to chromosome starts and ends
  # [i.e., those whose midpoints are within flankSize + (0.5*winSize)
  # of chromosome starts and ends]
  # (0.5*winSize) represents half of the midpoint-window within which,
  # and around which, SNPs will be counted
  interSNPsGRchr <- interSNPsGRchr[( interSNPsGRchr$midpoint >
                                       flankSize +
                                       (0.5*winSize) ) &
                                   ( interSNPsGRchr$midpoint <
                                       chrLens[i] -
                                       flankSize -
                                       (0.5*winSize) )]
  interSNPsGR <- append(interSNPsGR, interSNPsGRchr)
}

# Find and remove inter-SNP intervals that overlap crossover intervals
COs_interSNPs_overlap <- findOverlaps(query = COsGR,
                                      subject = interSNPsGR,
                                      type = "any",
                                      select = "all",
                                      ignore.strand = TRUE)
if(length(COs_interSNPs_overlap) > 0) {
  interSNPsGR <- interSNPsGR[-subjectHits(COs_interSNPs_overlap)]
}
## OR
#interSNPsGR2 <- subsetByOverlaps(x = interSNPsGR,
#                                 ranges = COsGR,
#                                 type = "any",
#                                 invert = TRUE)

# Define function to randomly select qualifying SNP intervals from interSNPsGR,
# with the same number per chromosome as COsGR
ranLocSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = TRUE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Define seed so that random selections are reproducible
set.seed(93750174)

# Apply ranLocSelect() on a per-chromosome basis so that
# ranLocGR contains the same number of loci per chromosome as COsGR
ranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  interSNPsGRchr <- interSNPsGR[seqnames(interSNPsGR) == chrs[i]]
  COsGRchr <- COsGR[seqnames(COsGR) == chrs[i]]
  ranLocChrGR <- ranLocSelect(coordinates = interSNPsGRchr,
                              n = length(COsGRchr))
  ranLocGR <- append(ranLocGR, ranLocChrGR)
}
print(range(width(COsGR)))
print(range(width(ranLocGR)))

# Plot crossover and random loci width distributions
pdf(paste0(plotDir, popName, "_COs_and_ranLoc_width_distributions.pdf"),
    height = 6, width = 10)
par(mfrow = c(1, 2),
    mar = c(6, 6, 2, 2),
    mgp = c(4, 1.5, 0))
hist(width(COsGR),
     breaks = 200,
     col = "grey50",
     border = NA,
     lwd = 2,
     xlim = c(min(c(width(COsGR), width(ranLocGR))),
              max(c(width(COsGR), width(ranLocGR)))),
     ylim = c(0, 250),
     xlab = "Width (bp)",
     ylab = "Intervals",
     main = "Crossovers",
     cex.lab = 2, cex.axis = 2)
abline(v = mean(width(COsGR)),
       col = "red", lty = 2, lwd = 2)
hist(width(ranLocGR),
     breaks = 200,
     col = "grey50",
     border = NA,
     lwd = 2,
     xlim = c(min(c(width(COsGR), width(ranLocGR))),
              max(c(width(COsGR), width(ranLocGR)))),
     ylim = c(0, 250),
     xlab = "Width (bp)",
     ylab = "Intervals",
     main = "Random SNP intervals",
     cex.lab = 2, cex.axis = 2)
abline(v = mean(width(ranLocGR)),
       col = "red", lty = 2, lwd = 2)
dev.off()

KSwidthDist <- ks.test(x = width(COsGR), y = width(ranLocGR))
print(paste0("K-S test P-value = ", KSwidthDist$p.value))
if(KSwidthDist$p.value < 0.05) {
  print("K-S test P-value is < 0.05; crossover and random loci width distributions differ")
}


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
                          flankName, "_flank_", winName, "_win_dataframe.txt"),
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
                          flankName, "_flank_", winName, "_win_dataframe.txt"),
            col.names = T, row.names = F, quote = F, sep = "\t")

print(paste0(popName, " COs and ranLoc SNP frequency profile calculation complete"))

# Function for plotting average SNP frequency profiles
plotAvgSNPfreq <- function(dat1, dat2,
                           col1, col2,
                           mainTitle,
                           flankSize, winSize,
                           flankLabL1, flankLabR1,
                           flankLabL2, flankLabR2,
                           legendLoc, legendLabs) {
  plot(x = 1:(((flankSize*2)+winSize)/winSize),
       y = dat1, col = col1,
       type = "l", lwd = 3, ann = F,
       ylim = c(min(dat1, dat2),
                max(dat1, dat2)),
       xaxt = "n", yaxt = "n")
  lines(x = 1:(((flankSize*2)+winSize)/winSize),
        y = dat2, col = col2,
        type = "l", lwd = 3)
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle)
  axis(side = 2, at = pretty(c(dat1, dat2)), lwd = 2, cex.axis = 0.7)
  mtext(side = 2, line = 2.0, cex = 0.8,
        text = bquote("SNP frequency" ~
                      .(as.character(winSize)) ~ "bp"^-1))
  axis(side = 1, lwd = 2,
       at = c(1,
              ((flankSize/winSize)/2),
              ((flankSize+winSize)/winSize),
              (((flankSize*2)+winSize)/winSize)-((flankSize/winSize)/2),
              (((flankSize*2)+winSize)/winSize)),
       labels = c("", "", "", "", ""))
  mtext(side = 1, line = 0.5, cex = 0.7,
        at = c(1,
               ((flankSize/winSize)/2),
               ((flankSize+winSize)/winSize),
               (((flankSize*2)+winSize)/winSize)-((flankSize/winSize)/2),
               (((flankSize*2)+winSize)/winSize)),
        text = c(flankLabL1, flankLabL2, "Midpoint", flankLabR2, flankLabR1))  
  abline(v = (flankSize+winSize)/winSize, lty = 3, lwd = 2)
  box(lwd = 2)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col1, col2),
         text.col = c(col1, col2),
         text.font = c(1, 1),
         ncol = 1, cex = 0.7, lwd = 2, bty = "n")
}

# Plot
pdf(paste0(plotDir, popName, "_COs_and_ranLoc_SNP_frequency_",
           flankName, "_flank_", winName, "_win.pdf"),
    height = 4, width = 4.5)
par(mar = c(2.1, 3.2, 2.1, 2.1))
par(mgp = c(2.25, 0.75, 0))
plotAvgSNPfreq(dat1 = colMeans(COsSNPcountsDF),
               dat2 = colMeans(ranLocSNPcountsDF),
               col1 = "forestgreen" ,
               col2 = "grey50",
               mainTitle = bquote(.(popName) ~ "(" *
                                  .(prettyNum(dim(COsSNPcountsDF)[1],
                                              big.mark = ",", trim = T)) *
                                  " crossovers)"),
               flankSize = flankSize,
               winSize = winSize,
               flankLabL1 = paste0("-", as.character(flankSize/1000), " kb"),
               flankLabR1 = paste0(as.character(flankSize/1000), " kb"),
               flankLabL2 = paste0("-", as.character((flankSize/1000)/2), " kb"),
               flankLabR2 = paste0(as.character((flankSize/1000)/2), " kb"),
               legendLoc = "bottomleft",
               legendLabs = c("Crossovers", "Random"))
dev.off()
