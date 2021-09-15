#!/applications/R/R-3.3.2/bin/Rscript

#######################################################################################
# Base relative frequency around target and random loci                               #
#######################################################################################

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/baseRelativeFrequency.R")
source("/projects/ajt200/Rfunctions/baseRelativeFrequencyPlot.R")
library(Biostrings)
library(segmentSeq)
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(regioneR)

# Chromosome definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
mask <- toGRanges(data.frame(rep(chrs, 2),
                             c(rep(1, 5), chrLens-5000),
                             c(rep(5000, 5), chrLens)))

outDir <- "./"
plotDir <- "./plots/"

# Import GBS COs as GRanges object
load("/projects/ajt200/GBS_CO/HS_CU_080617/wt/COsGRcoords.RData")
targetsGR <- COsGRcoords
print(length(targetsGR))

# Define target loci midpoints around which flanking sequence base relative frequencies will be calculated
targetMidpointsGR <- GRanges(seqnames = seqnames(targetsGR),
                             ranges = IRanges(start = round((start(targetsGR))+((end(targetsGR)-start(targetsGR))/2)),
                                              end = round((start(targetsGR))+((end(targetsGR)-start(targetsGR))/2))),
                             strand = "*")

# Chromosome sequence definitions
chr1 <- Athaliana$Chr1
chr2 <- Athaliana$Chr2
chr3 <- Athaliana$Chr3
chr4 <- Athaliana$Chr4
chr5 <- Athaliana$Chr5

# Use baseRelativeFreq() function to calculate base relative frequency around target loci and random loci
# For stand-aware analysis (e.g., around TSS or TTS), must separate into plus- and minus-strand loci
baseRelativeFreq(targets = targetMidpointsGR, genome = genome, mask = mask,
                 flankSize = 5000, locSize = 10001,
                 outDir = outDir, targetName = "GBS_crossovers")

# Load base relative frequency files
target.baseFreq <- list(read.table(file = paste0(outDir, "GBS_crossovers_target.a.txt"))[,1],
                        read.table(file = paste0(outDir, "GBS_crossovers_target.t.txt"))[,1],
                        read.table(file = paste0(outDir, "GBS_crossovers_target.g.txt"))[,1],
                        read.table(file = paste0(outDir, "GBS_crossovers_target.c.txt"))[,1])
ranLoc.baseFreq <- list(read.table(file = paste0(outDir, "GBS_crossovers_ranLoc.a.txt"))[,1],
                        read.table(file = paste0(outDir, "GBS_crossovers_ranLoc.t.txt"))[,1],
                        read.table(file = paste0(outDir, "GBS_crossovers_ranLoc.g.txt"))[,1],
                        read.table(file = paste0(outDir, "GBS_crossovers_ranLoc.c.txt"))[,1])

# Target and random locus summed A+T and G+C relative frequencies
target.ATfreq <- target.baseFreq[[1]]+target.baseFreq[[2]]
target.GCfreq <- target.baseFreq[[3]]+target.baseFreq[[4]]
ranLoc.ATfreq <- ranLoc.baseFreq[[1]]+ranLoc.baseFreq[[2]]
ranLoc.GCfreq <- ranLoc.baseFreq[[3]]+ranLoc.baseFreq[[4]]

ma <- function(x, n = 101) {
  filter(x, filter = rep(1/n, n), sides = 2, circular = T)
}

target.baseFreq <- lapply(seq_along(target.baseFreq), function(x) {
  ma(target.baseFreq[[x]])
})
ranLoc.baseFreq <- lapply(seq_along(ranLoc.baseFreq), function(x) {
  ma(ranLoc.baseFreq[[x]])
})

target.ATfreq <- ma(target.ATfreq)
target.GCfreq <- ma(target.GCfreq)
ranLoc.ATfreq <- ma(ranLoc.ATfreq)
ranLoc.GCfreq <- ma(ranLoc.GCfreq)

pdf(paste0(plotDir, "GBS_crossovers_base_relative_frequency_diffYaxes_v250619.pdf"), height = 5, width = 6)
par(mfrow = c(2, 2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))
flankSize <- 5000
xplot <- seq(-flankSize, flankSize, by = 1)
mycols <- c("magenta2", "springgreen2")
mergeBaseFreqPlotDiffY(at.coords = target.ATfreq, gc.coords = target.GCfreq,
                       at.ran.coords = ranLoc.ATfreq, gc.ran.coords = ranLoc.GCfreq,
                       flankSize = 5000, flankLabL = "-5 kb", flankLabR = "+5 kb",
                       midpointLab1 = "Midpoint", midpointLab2 = "Midpoint", mycols = mycols, xplot = xplot,
                       mainTitle1 = bquote("Crossovers ("*italic("n")*" = "*.(prettyNum(length(targetsGR),
                                                                                        big.mark = ",", trim = T))*")"),
                       mainTitle2 = bquote("Random loci ("*italic("n")*" = "*.(prettyNum(length(targetsGR),
                                                                                         big.mark = ",", trim = T))*")"))
mycols2 <- c("springgreen2", "springgreen4", "magenta2", "magenta4")
baseFreqPlotDiffY(coords = target.baseFreq,
                  ran.coords = ranLoc.baseFreq,
                  flankSize = 5000, flankLabL = "-5 kb", flankLabR = "+5 kb",
                  midpointLab1 = "Midpoint", midpointLab2 = "Midpoint", mycols = mycols2, xplot = xplot,
                  mainTitle1 = bquote("Crossovers ("*italic("n")*" = "*.(prettyNum(length(targetsGR),
                                                                                   big.mark = ",", trim = T))*")"),
                  mainTitle2 = bquote("Random loci ("*italic("n")*" = "*.(prettyNum(length(targetsGR),
                                                                                    big.mark = ",", trim = T))*")"))
dev.off()

pdf(paste0(plotDir, "GBS_crossovers_peak_base_relative_frequency_v250619.pdf"), height = 5, width = 6)
par(mfrow = c(2, 2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))
flankSize <- 5000
xplot <- seq(-flankSize, flankSize, by = 1)
mycols <- c("magenta2", "springgreen2")
mergeBaseFreqPlot(at.coords = target.ATfreq, gc.coords = target.GCfreq,
                  at.ran.coords = ranLoc.ATfreq, gc.ran.coords = ranLoc.GCfreq,
                  flankSize = 5000, flankLabL = "-5 kb", flankLabR = "+5 kb",
                  midpointLab1 = "Midpoint", midpointLab2 = "Midpoint", mycols = mycols, xplot = xplot,
                  mainTitle1 = bquote("Crossovers ("*italic("n")*" = "*.(prettyNum(length(targetsGR),
                                                                                   big.mark = ",", trim = T))*")"),
                  mainTitle2 = bquote("Random loci ("*italic("n")*" = "*.(prettyNum(length(targetsGR),
                                                                                    big.mark = ",", trim = T))*")"))
mycols2 <- c("springgreen2", "springgreen4", "magenta2", "magenta4")
baseFreqPlot(coords = target.baseFreq,
             ran.coords = ranLoc.baseFreq,
             flankSize = 5000, flankLabL = "-5 kb", flankLabR = "+5 kb",
             midpointLab1 = "Midpoint", midpointLab2 = "Midpoint", mycols = mycols2, xplot = xplot,
             mainTitle1 = bquote("Crossovers ("*italic("n")*" = "*.(prettyNum(length(targetsGR),
                                                                              big.mark = ",", trim = T))*")"),
             mainTitle2 = bquote("Random loci ("*italic("n")*" = "*.(prettyNum(length(targetsGR),
                                                                               big.mark = ",", trim = T))*")"))
dev.off()

pdf(paste0(plotDir, "GBS_crossovers_peak_base_relative_frequency_AT_and_GC_diffYlims_v250619.pdf"), height = 10, width = 6)
par(mfrow = c(4, 2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))
flankSize <- 5000
xplot <- seq(-flankSize, flankSize, by = 1)
mycols <- c("magenta2", "springgreen2")
ATorGCmergeBaseFreqPlot(coords = target.ATfreq,
                        ran.coords = ranLoc.ATfreq,
                        flankSize = 5000, flankLabL = "-5 kb", flankLabR = "+5 kb",
                        midpointLab1 = "Midpoint", midpointLab2 = "Midpoint", mycols = mycols[2], xplot = xplot,
                        mainTitle1 = bquote("Crossovers ("*italic("n")*" = "*.(prettyNum(length(targetsGR),
                                                                                         big.mark = ",", trim = T))*")"),
                        mainTitle2 = bquote("Random loci ("*italic("n")*" = "*.(prettyNum(length(targetsGR),
                                                                                          big.mark = ",", trim = T))*")"),
                        legendLab = "A+T")
ATorGCmergeBaseFreqPlot(coords = target.GCfreq,
                        ran.coords = ranLoc.GCfreq,
                        flankSize = 5000, flankLabL = "-5 kb", flankLabR = "+5 kb",
                        midpointLab1 = "Midpoint", midpointLab2 = "Midpoint", mycols = mycols[1], xplot = xplot,
                        mainTitle1 = bquote("Crossovers ("*italic("n")*" = "*.(prettyNum(length(targetsGR),
                                                                                         big.mark = ",", trim = T))*")"),
                        mainTitle2 = bquote("Random loci ("*italic("n")*" = "*.(prettyNum(length(targetsGR),
                                                                                          big.mark = ",", trim = T)
)*")"),
                        legendLab = "G+C")
mycols2 <- c("springgreen2", "springgreen4", "magenta2", "magenta4")
ATorGCbaseFreqPlot(coords = list(target.baseFreq[[1]], target.baseFreq[[2]]),
                   ran.coords = list(ranLoc.baseFreq[[1]], ranLoc.baseFreq[[2]]),
                   flankSize = 5000, flankLabL = "-5 kb", flankLabR = "+5 kb",
                   midpointLab1 = "Midpoint", midpointLab2 = "Midpoint", mycols = mycols2[1:2], xplot = xplot,
                   mainTitle1 = bquote("Crossovers ("*italic("n")*" = "*.(prettyNum(length(targetsGR),
                                                                                    big.mark = ",", trim = T))*")"),
                   mainTitle2 = bquote("Random loci ("*italic("n")*" = "*.(prettyNum(length(targetsGR),
                                                                                     big.mark = ",", trim = T)
)*")"),
                   legendLab = c("A", "T"))
ATorGCbaseFreqPlot(coords = list(target.baseFreq[[3]], target.baseFreq[[4]]),
                   ran.coords = list(ranLoc.baseFreq[[3]], ranLoc.baseFreq[[4]]),
                   flankSize = 5000, flankLabL = "-5 kb", flankLabR = "+5 kb",
                   midpointLab1 = "Midpoint", midpointLab2 = "Midpoint", mycols = mycols2[3:4], xplot = xplot,
                   mainTitle1 = bquote("Crossovers ("*italic("n")*" = "*.(prettyNum(length(targetsGR),
                                                                                    big.mark = ",", trim = T))*")"),
                   mainTitle2 = bquote("Random loci ("*italic("n")*" = "*.(prettyNum(length(targetsGR),
                                                                                     big.mark = ",", trim = T)
)*")"),
                   legendLab = c("G", "C"))
dev.off()
