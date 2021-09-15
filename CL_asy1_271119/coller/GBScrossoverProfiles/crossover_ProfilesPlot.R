#!/applications/R/R-3.5.0/bin/Rscript

# Plot average coverage around target and random loci

# Usage:
# ./crossover_ProfilesPlot.R 'SPO11_1_oligos_RPI1,MNase' 'SPO11-1-oligos,Nucleosomes' 'red,forestgreen' 'coller,colws,asy1hom,asy1het' 5000 5kb '5 kb' 50

#libNames <- unlist(strsplit("SPO11_1_oligos_RPI1,MNase",
#                            split = ","))
#libNamesPlot <- unlist(strsplit("SPO11-1-oligos,Nucleosomes",
#                                split = ","))
#colours <- unlist(strsplit("red,forestgreen",
#                           split = ","))
#genotypes <- unlist(strsplit("coller,colws,asy1hom,asy1het",
#                             split = ","))
#flankSize <- 5000
#flankName <- "5kb"
#flankNamePlot <- "5 kb"
#winSize <- 50

args <- commandArgs(trailingOnly = T)
libNames <- unlist(strsplit(args[1],
                            split = ","))
libNamesPlot <- unlist(strsplit(args[2],
                                split = ","))
colours <- unlist(strsplit(args[3],
                           split = ","))
genotypes <- unlist(strsplit(args[4],
                             split = ","))
flankSize <- as.numeric(args[5])
flankName <- as.character(args[6])
flankNamePlot <- as.character(args[7])
winSize <- as.numeric(args[8])

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/plotAvgCov_plotAvgMeth_target_ranLoc.R")

inDir <- "/projects/ajt200/GBS_CO/CL_asy1_271119/"
plotDir <- "./plots/"

# Load column mean coverage profiles
# The first level in this list of lists of lists corresponds to libNames,
# the second level in this list of lists of lists corresponds to genotypes,
# and the third level in this list of list of lists corresponds to features and random loci
avgProfilesList <- lapply(seq_along(libNames), function(x) {
  lapply(seq_along(genotypes), function(y) {
    list(
         read.table(paste0(inDir, genotypes[y],
                           "/GBScrossoverProfiles/matrices/",
                           libNames[x], "_norm_cov_feature_smoothed_target_midpoint_and_",
                           flankName, "_flank_dataframe_colMeans.txt"),
                    header = T)[,1],
         read.table(paste0(inDir, genotypes[y],
                           "/GBScrossoverProfiles/matrices/",
                           libNames[x], "_norm_cov_ranLoc_smoothed_target_midpoint_and_",
                           flankName, "_flank_dataframe_colMeans.txt"),
                    header = T)[,1]
    )
  })
})

ymin_list <- c(
  lapply(seq_along(avgProfilesList), function(x) {
    min(sapply(seq_along(avgProfilesList[[x]]), function(y) {
      min(c(avgProfilesList[[x]][[y]][[1]],
            avgProfilesList[[x]][[y]][[2]]))
    }))
  })
)
ymax_list <- c(
  lapply(seq_along(avgProfilesList), function(x) {
    max(sapply(seq_along(avgProfilesList[[x]]), function(y) {
      max(c(avgProfilesList[[x]][[y]][[1]],
            avgProfilesList[[x]][[y]][[2]]))
    }))
  })
)

pdf(paste0(plotDir,
           "GBS_crossover_profiles_",
           paste0(genotypes, collapse = "_"), "_",
           paste0(libNames, collapse = "_"),
           ".pdf"),
    height = 2.5*length(genotypes), width = 6)
par(mfrow = c(length(genotypes), 2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))

for(y in seq_along(genotypes)) {
  plotAvgCovMidpoint_1v2_twoY_Ylims(featureTitle = bquote(.(genotypes[y])~"crossovers"),
                                    xplot = seq(from = 1,
                                                to = length(avgProfilesList[[1]][[y]][[1]]),
                                                by = 1),
                                    dat1 = avgProfilesList[[1]][[y]][[1]],
                                    dat2 = avgProfilesList[[2]][[y]][[1]],
                                    ranDat1 = avgProfilesList[[1]][[y]][[2]],
                                    ranDat2 = avgProfilesList[[2]][[y]][[2]],
                                    col1 = colours[1], col2 = colours[2],
                                    Ylim1 = c(ymin_list[[1]], ymax_list[[1]]),
                                    Ylim2 = c(ymin_list[[2]], ymax_list[[2]]),
                                    Ylabel1 = libNamesPlot[1],
                                    Ylabel2 = libNamesPlot[2],
                                    winSize = winSize,
                                    flankSize = flankSize,
                                    flankLabL = paste0("-", flankNamePlot),
                                    flankLabR = paste0("+", flankNamePlot))
}
dev.off()
