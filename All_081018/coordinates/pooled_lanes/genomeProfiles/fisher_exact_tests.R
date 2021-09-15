#!/applications/R/R-3.5.0/bin/Rscript

# Perform a Fisher's exact test on crossover counts in each genomic window
# for two genotypes

# Usage from within directory containing windowed CO counts files:
# ./fisher_exact_tests.R HS_wt_ColLer_F2_2018 SB_msh2_ColLer_F2_2018 /projects/ajt200/GBS_CO/SB_060718/SB_msh2_wt/regioneR/SB_ColLer_F2_2018_maskGR.RData 500kb 0.1

library(GenomicRanges)

#pop1Name <- "HS_wt_ColLer_F2_2018"
#pop2Name <- "SB_msh2_ColLer_F2_2018"
#maskFile <- "/projects/ajt200/GBS_CO/SB_060718/SB_msh2_wt/regioneR/SB_ColLer_F2_2018_maskGR.RData"
#winName <- "500kb"
#FDR <- 0.1
#FDRname <- paste0("FDR", as.character(FDR))

args <- commandArgs(trailingOnly = TRUE)
pop1Name <- args[1]
pop2Name <- args[2]
maskFile <- args[3]
winName <- as.character(args[4])
FDR <- as.numeric(args[5])
FDRname <- paste0("FDR", as.character(FDR))

CODir <- "/projects/ajt200/GBS_CO/All_081018/coordinates/pooled_lanes/"
load(paste0(CODir, pop1Name, "_pooled_GRanges.RData"))
if( exists("popCOsGR") == TRUE ) {
popCOsGR_pooled <- popCOsGR
}
print(pop1Name)
print("Individuals with crossovers:")
print(length(unique(popCOsGR_pooled$lib)))
pop1Ind <- length(unique(popCOsGR_pooled$lib))
popCOsGR_pooled <- NULL
popCOsGR <- NULL
load(paste0(CODir, pop2Name, "_pooled_GRanges.RData"))
if( length(popCOsGR) != 0 ) {
popCOsGR_pooled <- popCOsGR
}
print(pop2Name)
print("Individuals with crossovers:")
print(length(unique(popCOsGR_pooled$lib)))
pop2Ind <- length(unique(popCOsGR_pooled$lib))
popCOsGR_pooled <- NULL

inDir <- "./"
outDir <- paste0(inDir, "fisher_exact_tests/")
plotDir <- paste0(outDir, FDRname, "/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))
pop1COdensity <- read.table(paste0(inDir, pop1Name,
                                   "_pooled_CO_density_genome_", winName, ".txt"))
pop2COdensity <- read.table(paste0(inDir, pop2Name,
                                   "_pooled_CO_density_genome_", winName, ".txt"))

# Function to create a 2x2 contingency table where pop1InWin and pop1OutWin are the
# number of crossovers in pop1 located within and outside of a given genomic window
contingencyTable <- function(pop1InWin, pop1OutWin, pop2InWin, pop2OutWin) {
  conTab <- matrix(c(pop1InWin, pop1OutWin, pop2InWin, pop2OutWin), ncol = 2)
  rownames(conTab) <- c("InWin", "OutWin")
  colnames(conTab) <- c(pop1Name, pop2Name)
  conTab
}

# Perform Fisher's exact test on each contingency table and extract P-values
fisherPvals <- sapply(1:dim(pop1COdensity)[1], function(x) {
  fisher.test(contingencyTable(pop1InWin = pop1COdensity$COsChrWin[x],
                               pop1OutWin = sum(pop1COdensity$COsChrWin)-pop1COdensity$COsChrWin[x],
                               pop2InWin = pop2COdensity$COsChrWin[x],
                               pop2OutWin = sum(pop2COdensity$COsChrWin)-pop2COdensity$COsChrWin[x]),
              alternative = "two.sided")$p.value
})

# Alternative way to make contingency table and apply test
## Function to create a 2x2 contingency table where pop1InWin and pop1OutWin are the
## number of crossovers in pop1 located within and outside of a given genomic window
#contingencyTable <- function(pop1InWin, pop1OutWin, pop2InWin, pop2OutWin) {
#  conTab <- matrix(c(pop1InWin, pop2InWin, pop1OutWin, pop2OutWin),
#                   nrow = 2,
#                   dimnames = list(c(pop1Name, pop2Name),
#                                   c("InWin", "OutWin")))
#  conTab
#}
#
## Perform Fisher's exact test on each contingency table and extract P-values
#fisherPvals <- lapply(1:dim(pop1COdensity)[1], function(x) {
#  fisher.test(contingencyTable(pop1InWin = pop1COdensity$COsChrWin[x],
#                               pop2InWin = pop2COdensity$COsChrWin[x],
#                               pop1OutWin = sum(pop1COdensity$COsChrWin)-pop1COdensity$COsChrWin[x],
#                               pop2OutWin = sum(pop2COdensity$COsChrWin)-pop2COdensity$COsChrWin[x]),
#              alternative = "two.sided")$p.value
#})

# Adjust P-values by applying the Benjamini-Hochberg multiple testing correction
fisherAdjPvals <- p.adjust(p = fisherPvals, method = "BH")

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

load(maskFile)
maskStartsCum <- NULL
maskEndsCum <- NULL
for(x in 1:length(chrs)) {
  maskGRchr <- maskGR[seqnames(maskGR) == chrs[x]]
  maskGRchrStarts <- start(maskGRchr)
  maskGRchrStartsCum <- maskGRchrStarts+sumchr[x]
  maskStartsCum <- c(maskStartsCum, maskGRchrStartsCum)
  maskGRchrEnds <- end(maskGRchr)
  maskGRchrEndsCum <- maskGRchrEnds+sumchr[x]
  maskEndsCum <- c(maskEndsCum, maskGRchrEndsCum)
}

# Function to plot genome-scale profile of pop1 vs pop2 (one Y-axis) 
pop1Vpop2GenomePlot <- function(xplot,
                                pop1, pop1Col,
                                pop2, pop2Col,
                                Ylabel,
                                legendLoc,
                                legendLabs) {
  plot(xplot, pop1, type = "l", lwd = 1.5, col = pop1Col,
       ylim = c(min(c(pop1, pop2)),
                max(c(pop1, pop2))),
       xlab = "",
       ylab = "", yaxt = "n",
       main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(pop1, pop2, method = "spearman", use = "pairwise.complete.obs"), digits = 2))))
  lines(xplot, pop2, type = "l", lwd = 1.5, col = pop2Col)
  axis(side = 1, cex.axis = 1, lwd.tick = 1.5)
  axis(side = 2, at = pretty(c(pop1, pop2)), cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, text = Ylabel, col = "black")
  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75, col = "black")
  abline(v = c(maskStartsCum, maskEndsCum), lty = 1, lwd = 0.75, col = "orange")
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(pop1Col, pop2Col),
         text.col = c(pop1Col, pop2Col),
         text.font = c(1, 1),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

# Function to plot -log10-transformed P-values and adjusted P-values derived from Fisher's exact tests
minusLog10PvalPlot <- function(xplot,
                               Pvals, PvalsCol,
                               AdjPvals, AdjPvalsCol) {
  plot(x = xplot, y = -log10(Pvals), col = PvalsCol, type = "l", lwd = 1.5,
       ylim = c(min(-log10(Pvals)),
                max(-log10(Pvals))),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = "")
  mtext(side = 2, line = 2.25, cex = 1, col = PvalsCol,
        text = bquote("-Log"[10]*"("*italic("P")*"-value) "*.(winName)^-1))
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  abline(h = -log10(0.05), lty = 5, lwd = 1, col = PvalsCol)

  par(new = T)
  plot(x = xplot, y = -log10(AdjPvals), col = AdjPvalsCol, type = "l", lwd = 1.5,
       ylim = c(min(-log10(AdjPvals)),
                max(-log10(AdjPvals))),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(pop1Name)*" versus "*.(pop2Name)*" crossover counts per "*.(winName)*" window"),
       cex.main = 1)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1, adj = c(0.5, -2.5), xpd = NA, srt = -90, col = AdjPvalsCol,
       labels = bquote("-Log"[10]*"(BH-adjusted "*italic("P")*"-value) "*.(winName)^-1))
  axis(side = 4, cex.axis = 1, lwd.tick = 1.5)
  axis(side = 1, cex.axis = 1, lwd.tick = 1.5)
  abline(h = -log10(FDR), lty = 5, lwd = 1, col = AdjPvalsCol)

  abline(v = sumchr, lty = 1, lwd = 0.75)
  abline(v = centromeres, lty = 5, lwd = 0.75, col = "black")
  abline(v = c(maskStartsCum, maskEndsCum), lty = 1, lwd = 0.75, col = "orange")
  box(lwd = 1.5)
}

pdf(file = paste0(plotDir, "Fisher_exact_test_Pvals_genomeProfile_", pop1Name, "_vs_", pop2Name, "_", winName, "_", FDRname, ".pdf"),
    height = 6, width = 9)
par(mfrow = c(2, 1))
par(mar = c(2.1, 4.1, 2.1, 4.1))
par(mgp = c(3, 1, 0))

minusLog10PvalPlot(xplot = pop1COdensity$cumWindows,
                   Pvals = fisherPvals,
                   AdjPvals = fisherAdjPvals,
                   PvalsCol = "dodgerblue3",
                   AdjPvalsCol = "red")
pop1Vpop2GenomePlot(xplot = pop1COdensity$cumWindows,
                    pop1 = pop1COdensity$COsChrWin/pop1Ind,
                    pop2 = pop2COdensity$COsChrWin/pop2Ind,
                    Ylabel = bquote("Normalised crossover density "*.(winName)^-1),
                    legendLoc = "topleft",
                    legendLabs = c(pop1Name, pop2Name),
                    pop1Col = "grey50",
                    pop2Col = "forestgreen")
dev.off()

