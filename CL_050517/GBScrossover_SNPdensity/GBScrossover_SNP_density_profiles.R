#!/applications/R/R-3.5.0/bin/Rscript

library(GenomicRanges)

flankSize <- 5000
flankName <- "5kb"

SNPs <- read.table("/home/cal66/GBS/GBS_WT_asy1homo_Col_Ws_summary/SNPs_Col_Ws_190125.txt")
SNPs <- SNPs[,-4]
SNPsGR <- makeGRangesFromDataFrame(SNPs,
                                   ignore.strand = T)

COs <- read.table("/home/cal66/GBS/GBS_WT_asy1homo_Col_Ws_summary/CO_fine_map/asy1homo_Col_Ws_COs_midpoints_GR_width_LESS_1kb_190226.txt")
stopifnot(identical(COs$start, COs$end))
COsGR <- GRanges(seqnames = COs$seqnames,
                 ranges = IRanges(start = COs$start-flankSize,
                                  end = COs$end+flankSize),
                 strand = "*")

SNPcountsMatrix <- NULL
for(x in seq_along(COsGR)) {
  # Create GRanges object containing one range for each
  # position around a crossover midpoint
  posGR <- GRanges(seqnames = seqnames(COsGR[x]),
                   ranges = IRanges(start = start(COsGR[x]):end(COsGR[x]),
                                    end = start(COsGR[x]):end(COsGR[x])),
                   strand = "*")
  # Count number of SNPs overlapping each position
  overlaps <- countOverlaps(query = posGR,
                            subject = SNPsGR)
  # Add row of SNP counts for this crossover to the matrix
  SNPcountsMatrix <- rbind(SNPcountsMatrix, overlaps)
}

# Calculate average SNP frequency profile across all crossovers
avgSNPfreq <- colMeans(SNPcountsMatrix)

# Calculate moving average of current window,
# (N/2)-0.5 previous windows (where N is odd),
# and
# (N/2)-0.5 subsequent windows (where N is odd)
# (the higher N is, the greater the smoothing)
N <- 101
stopifnot(N %% 2 != 0)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

filt_avgSNPfreq <- stats::filter(x = avgSNPfreq,
                                 filter = f,
                                 sides = 2)
filt_avgSNPfreq[1:flank] <- filt_avgSNPfreq[flank+1]
filt_avgSNPfreq[(length(filt_avgSNPfreq)-flank+1):length(filt_avgSNPfreq)] <- filt_avgSNPfreq[(length(filt_avgSNPfreq)-flank)]

pdf("./test.pdf")
plot(x = seq_along(1:width(COsGR[1])),
     y = filt_avgSNPfreq,
     type = "l",
     lwd = 3,
     col = "red")
abline(v = flankSize+1, lty = 3, lwd = 3)
dev.off()

