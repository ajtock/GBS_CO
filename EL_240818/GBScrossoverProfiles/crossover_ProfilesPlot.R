#!/applications/R/R-3.5.0/bin/Rscript

# Plot average coverage around target and random loci

# Usage:
# ./crossover_ProfilesPlot.R Col_RNAseq_meiocyte_Rep1 'Col-0' taf4b_RNAseq_meiocyte_Rep1 'taf4b-1' red dodgerblue3 5000 5kb '5 kb' 20 Rep1

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
RepNo <- args[11]

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/plotAvgCov_plotAvgMeth_target_ranLoc.R")

matDir <- "./matrices/"
plotDir <- "./plots/"

libNames <- c(libName1, libName2)
libNamesPlot <- c(libName1Plot, libName2Plot)

# wt crossovers
# Define column mean coverage outfile (mean profiles)
wt_outDFCM <- lapply(seq_along(libNames), function(x)
  list(paste0(matDir, libNames[[x]],
              "_norm_cov_wt_COs_mat1_smoothed_midpoint_and_",
              flankName, "flank_dataframe_colMeans.txt"),
       paste0(matDir, libNames[[x]],
              "_norm_cov_wt_ranLoc_mat2_smoothed_midpoint_and_",
              flankName, "flank_dataframe_colMeans.txt")))

# Read in target and ranLoc mean coverage profiles
wt_target_covDat <- lapply(seq_along(libNames), function(x) {
  read.table(file = wt_outDFCM[[x]][[1]])
})
wt_ranLoc_covDat <- lapply(seq_along(libNames), function(x) {
  read.table(file = wt_outDFCM[[x]][[2]])
})

# taf4b crossovers
# Define column mean coverage outfile (mean profiles)
taf4b_outDFCM <- lapply(seq_along(libNames), function(x)
  list(paste0(matDir, libNames[[x]],
              "_norm_cov_taf4b_COs_mat1_smoothed_midpoint_and_",
              flankName, "flank_dataframe_colMeans.txt"),
       paste0(matDir, libNames[[x]],
              "_norm_cov_taf4b_ranLoc_mat2_smoothed_midpoint_and_",
              flankName, "flank_dataframe_colMeans.txt")))

# Read in target and ranLoc mean coverage profiles
taf4b_target_covDat <- lapply(seq_along(libNames), function(x) {
  read.table(file = taf4b_outDFCM[[x]][[1]])
})
taf4b_ranLoc_covDat <- lapply(seq_along(libNames), function(x) {
  read.table(file = taf4b_outDFCM[[x]][[2]])
})

# Plot mean REC8 vs other coverage profiles around target and random loci
Ylim_inclRan <- c(min(c(wt_target_covDat[[1]][,1],
                        wt_target_covDat[[2]][,1],
                        wt_ranLoc_covDat[[1]][,1],
                        wt_ranLoc_covDat[[2]][,1],
                        taf4b_target_covDat[[1]][,1],
                        taf4b_target_covDat[[2]][,1],
                        taf4b_ranLoc_covDat[[1]][,1],
                        taf4b_ranLoc_covDat[[2]][,1]),
                      na.rm = T),
                  max(c(wt_target_covDat[[1]][,1],
                        wt_target_covDat[[2]][,1],
                        wt_ranLoc_covDat[[1]][,1],
                        wt_ranLoc_covDat[[2]][,1],
                        taf4b_target_covDat[[1]][,1],
                        taf4b_target_covDat[[2]][,1],
                        taf4b_ranLoc_covDat[[1]][,1],
                        taf4b_ranLoc_covDat[[2]][,1]),
                      na.rm = T))
Ylim_exclRan <- c(min(c(wt_target_covDat[[1]][,1],
                        wt_target_covDat[[2]][,1],
                        taf4b_target_covDat[[1]][,1],
                        taf4b_target_covDat[[2]][,1]),
                      na.rm = T),
                  max(c(wt_target_covDat[[1]][,1],
                        wt_target_covDat[[2]][,1],
                        taf4b_target_covDat[[1]][,1],
                        taf4b_target_covDat[[2]][,1]),
                      na.rm = T))

pdf(paste0(plotDir, "Col_and_taf4b_GBScrossover_profiles_",
           libName1, "_and_", libName2, ".pdf"),
    height = 7.5, width = 6)
par(mfrow = c(3, 2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))

plotAvgCovMidpoint_1v2_oneY_Ylim(featureTitle = bquote(.(libNamesPlot[1])~"crossovers"),
                                 xplot = seq(from = 1, to = length(wt_target_covDat[[1]][,1]), by = 1),
                                 dat1 = wt_target_covDat[[1]][,1],
                                 dat2 = wt_target_covDat[[2]][,1],
                                 ranDat1 = wt_ranLoc_covDat[[1]][,1],
                                 ranDat2 = wt_ranLoc_covDat[[2]][,1],
                                 col1 = colour1, col2 = colour2,
                                 Ylim = Ylim_inclRan,
                                 Ylabel = paste0("RNA-seq ", RepNo, " (CPM)"),
                                 winSize = winSize,
                                 flankSize = flankSize,
                                 flankLabL = paste0("-", flankNamePlot),
                                 flankLabR = paste0("+", flankNamePlot),
                                 legendLabs = libNamesPlot,
                                 legendLoc = "topleft")
plotAvgCovMidpoint_1v2_oneY_Ylim(featureTitle = bquote(italic(.(libNamesPlot[2]))~"crossovers"),
                                 xplot = seq(from = 1, to = length(taf4b_target_covDat[[1]][,1]), by = 1),
                                 dat1 = taf4b_target_covDat[[1]][,1],
                                 dat2 = taf4b_target_covDat[[2]][,1],
                                 ranDat1 = taf4b_ranLoc_covDat[[1]][,1],
                                 ranDat2 = taf4b_ranLoc_covDat[[2]][,1],
                                 col1 = colour1, col2 = colour2,
                                 Ylim = Ylim_inclRan,
                                 Ylabel = paste0("RNA-seq ", RepNo, " (CPM)"),
                                 winSize = winSize,
                                 flankSize = flankSize,
                                 flankLabL = paste0("-", flankNamePlot),
                                 flankLabR = paste0("+", flankNamePlot),
                                 legendLabs = libNamesPlot,
                                 legendLoc = "topleft")
plotAvgCovMidpoint_1v2_oneY_Ylim_noRan(featureTitle = bquote(.(libNamesPlot[1])~"crossovers"),
                                       xplot = seq(from = 1, to = length(wt_target_covDat[[1]][,1]), by = 1),
                                       dat1 = wt_target_covDat[[1]][,1],
                                       dat2 = wt_target_covDat[[2]][,1],
                                       col1 = colour1, col2 = colour2,
                                       Ylim = Ylim_exclRan,
                                       Ylabel = paste0("RNA-seq ", RepNo, " (CPM)"),
                                       winSize = winSize,
                                       flankSize = flankSize,
                                       flankLabL = paste0("-", flankNamePlot),
                                       flankLabR = paste0("+", flankNamePlot),
                                       legendLabs = libNamesPlot,
                                       legendLoc = "topleft")
plotAvgCovMidpoint_1v2_oneY_Ylim_noRan(featureTitle = bquote(italic(.(libNamesPlot[2]))~"crossovers"),
                                       xplot = seq(from = 1, to = length(taf4b_target_covDat[[1]][,1]), by = 1),
                                       dat1 = taf4b_target_covDat[[1]][,1],
                                       dat2 = taf4b_target_covDat[[2]][,1],
                                       col1 = colour1, col2 = colour2,
                                       Ylim = Ylim_exclRan,
                                       Ylabel = paste0("RNA-seq ", RepNo, " (CPM)"),
                                       winSize = winSize,
                                       flankSize = flankSize,
                                       flankLabL = paste0("-", flankNamePlot),
                                       flankLabR = paste0("+", flankNamePlot),
                                       legendLabs = libNamesPlot,
                                       legendLoc = "topleft")
dev.off()

