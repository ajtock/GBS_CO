#!/applications/R/R-3.5.0/bin/Rscript

# Plot average coverage around target and random loci

# Usage:
# ./crossover_ProfilesPlot.R SPO11_1_oligos_RPI1 'SPO11-1' MNase 'Nucleosomes' red darkcyan 5000 5kb '5 kb' 20

#libName1 <- "SPO11_1_oligos_RPI1"
#libName1Plot <- "SPO11-1"
#libName2 <- "MNase"
#libName2Plot <- "Nucleosomes"
#colour1 <- "red"
#colour2 <- "darkcyan"
#flankSize <- 5000
#flankName <- "5kb"
#flankNamePlot <- "5 kb"
#winSize <- 20

args <- commandArgs(trailingOnly = T)
libName1 <- args[1]
libName1Plot <- args[2]
libName2 <- args[3]
libName2Plot <- args[4]
colour1 <- args[5]
colour2 <- args[6]
flankSize <- as.numeric(args[7])
flankName <- as.character(args[8])
flankNamePlot <- as.character(args[9])
winSize <- as.numeric(args[10])

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/plotAvgCov_plotAvgMeth_target_ranLoc.R")

library(parallel)

matDir <- "matrices/"
plotDir <- "plots/"

libNames <- c(libName1, libName2)
libNamesPlot <- c(libName1Plot, libName2Plot)

# Obtain feature names
featureNames <- c(
                  "wt",
                  "hcr1"
                 )

# Define column mean coverage outfiles (profiles)
outDFcolMeansList <- lapply(seq_along(featureNames), function(x) {
  lapply(seq_along(libNames), function(y) {
    list(paste0(matDir, libNames[[y]],
                "_norm_cov_", featureNames[[x]], "_crossovers_smoothed_target_and_", flankName, "_flank_dataframe_colMeans.txt"),
         paste0(matDir, libNames[[y]],
                "_norm_cov_", featureNames[[x]], "_ranLoc_smoothed_target_and_", flankName, "_flank_dataframe_colMeans.txt"))
  })
})

# Read in target and ranLoc mean coverage profiles
targetCovDat <- mclapply(seq_along(featureNames), function(x) {
  lapply(seq_along(libNames), function(y) {
    read.table(file = outDFcolMeansList[[x]][[y]][[1]],
               header = T)[,1]
  })
}, mc.cores = length(featureNames))

ranLocCovDat <- mclapply(seq_along(featureNames), function(x) {
  lapply(seq_along(libNames), function(y) {
    read.table(file = outDFcolMeansList[[x]][[y]][[2]],
               header = T)[,1]
  })
}, mc.cores = length(featureNames))

# Plot mean libName1 vs libNlibName1 coverage profiles around target and random loci
Ymin <- NULL
Ymax <- NULL
for(x in 1:length(targetCovDat)) {
  YminFeature <- min(c(targetCovDat[[x]][[1]],
                       targetCovDat[[x]][[2]],
                       ranLocCovDat[[x]][[1]],
                       ranLocCovDat[[x]][[2]]))
  Ymin <- min(c(Ymin, YminFeature))
  YmaxFeature <- max(c(targetCovDat[[x]][[1]],
                       targetCovDat[[x]][[2]],
                       ranLocCovDat[[x]][[1]],
                       ranLocCovDat[[x]][[2]]))
  Ymax <- max(c(Ymax, YmaxFeature))
}

pdf(paste0(plotDir, "hcr1_paper_GBScrossover_profiles_",
           libName1, "_and_", libName2, "_v060121_Pohang.pdf"),
    height = 2.5*length(featureNames), width = 6)
par(mfrow = c(length(featureNames), 2))
par(mar = c(2.1, 3.4, 2.1, 3.4))
par(mgp = c(2.25, 1, 0))

for(x in 1:length(targetCovDat)) {
  plotAvgCovMidpoint_1v2_oneY_Ylim(featureTitle = bquote(.(featureNames[x])~"crossovers"),
                                   xplot = seq_along(targetCovDat[[x]][[1]]),
                                   dat1 = targetCovDat[[x]][[1]],
                                   dat2 = targetCovDat[[x]][[2]],
                                   ranDat1 = ranLocCovDat[[x]][[1]],
                                   ranDat2 = ranLocCovDat[[x]][[2]],
                                   col1 = colour1, col2 = colour2,
                                   Ylim = c(Ymin, Ymax),
                                   Ylabel = bquote("Log"[2]*"(signal/input)"),
                                   winSize = winSize,
                                   flankSize = flankSize,
                                   flankLabL = paste0("-", flankNamePlot),
                                   flankLabR = paste0("+", flankNamePlot),
                                   legendLabs = libNamesPlot,
                                   targetLegendLoc = "left",
                                   ranLocLegendLoc = "topleft")
}
dev.off()

