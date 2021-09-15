#!/applications/R/R-3.5.0/bin/Rscript

# Plot average coverage around target and random loci

# Usage:
# ./crossover_ProfilesPlot.R SPO11_1_oligos_RPI1 SPO11-1-oligos MNase Nucleosomes red darkcyan 5000 5kb '5 kb' 20

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

matDir <- "./matrices/"
plotDir <- "./plots/"

libNames <- c(libName1, libName2)
libNamesPlot <- c(libName1Plot, libName2Plot)

# Define column mean coverage outfile (mean profiles)
outDFCM <- lapply(seq_along(libNames), function(x)
  list(paste0(matDir, libNames[[x]],
              "_norm_cov_COs_mat1_smoothed_midpoint_and_",
              flankName, "flank_dataframe_colMeans.txt"),
       paste0(matDir, libNames[[x]],
              "_norm_cov_ranLoc_mat2_smoothed_midpoint_and_",
              flankName, "flank_dataframe_colMeans.txt")))

# Read in target and ranLoc mean coverage profiles
target_covDat <- lapply(seq_along(libNames), function(x) {
  read.table(file = outDFCM[[x]][[1]])
})
ranLoc_covDat <- lapply(seq_along(libNames), function(x) {
  read.table(file = outDFCM[[x]][[2]])
})


# Plot mean REC8 vs other coverage profiles around target and random loci
pdf(paste0(plotDir, "BethRowan_crossover_profiles_",
           libName1, "_and_", libName2, ".pdf"),
    height = 2.5, width = 6)
par(mfrow = c(1, 2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))

plotAvgCovMidpoint_1v2_oneY(xplot = seq(from = 1, to = length(target_covDat[[1]][,1]), by = 1),
                            dat1 = target_covDat[[1]][,1],
                            dat2 = target_covDat[[2]][,1],
                            ranDat1 = ranLoc_covDat[[1]][,1],
                            ranDat2 = ranLoc_covDat[[2]][,1],
                            col1 = colour1, col2 = colour2,
                            Ylabel = bquote("Z-standardized log"[2]*"(ChIP/control)"),
                            winSize = winSize,
                            flankSize = flankSize,
                            flankLabL = paste0("-", flankNamePlot),
                            flankLabR = paste0("+", flankNamePlot),
                            legendLabs = c(libNamesPlot),
                            legendLoc = "left")
dev.off()

