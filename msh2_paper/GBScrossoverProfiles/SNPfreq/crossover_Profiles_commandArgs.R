#!/applications/R/R-3.4.0/bin/Rscript

# Profile mean coverage around crossovers and random loci

# Usage via Condor submission system on node7:
# csmit -m 20G -c 1 "./crossover_Profiles_commandArgs.R 5000 5kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI1"

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc_R3.4.0.R")

library(EnrichedHeatmap)
library(genomation)
library(regioneR)

args <- commandArgs(trailingOnly = T)
flankSize <- as.numeric(args[1])
flankName <- as.character(args[2])
winSize <- as.numeric(args[3])
covDatPath <- as.character(args[4])
libName <- as.character(args[5])

matDir <- "./matrices/"
plotDir <- "./plots/"
system(paste0("[ -d ", matDir, " ] || mkdir ", matDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Chromosome definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrStart <- c(rep(1, 5))
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))

# Import crossovers as GRanges objects
# Obtain feature names
featureNames <- NULL
for(x in 1:length(system("ls *.cut.txt",
                         intern = T))) {
  featureNamesElement <- substr(system("ls *.cut.txt",
                               intern = T)[x],
                               start = 1,
                               stop = nchar(system("ls *.cut.txt",
                                                   intern = T)[x])-8)
  featureNames <- c(featureNames, featureNamesElement)
}

featureListTab <- mclapply(seq_along(featureNames), function(x) {
  read.table(paste0(featureNames[x], ".cut.txt"))
}, mc.cores = length(featureNames))

colnames(featureListTab[[4]]) <- c("lane", "lib",
                                   "chrs", "start", "stop", "cos",
                                   "width")
colnames(featureListTab[[6]]) <- c("chrs", "start", "stop", "cos",
                                   "width")
colnames(featureListTab[[7]]) <- c("lane", "lib",
                                   "chrs", "start", "stop", "cos",
                                   "width", "geno")
colnames(featureListTab[[8]]) <- c("lane", "lib",
                                   "chrs", "start", "stop", "cos",
                                   "width", "geno")
colnames(featureListTab[[9]]) <- c("lane", "lib",
                                   "chrs", "start", "stop", "cos",
                                   "width", "geno")
colnames(featureListTab[[10]]) <- c("chrs", "start", "stop", "cos",
                                    "Selected.420")

featuresListGR <- list()
for(x in 1:length(featureListTab)) {
  COsGR <- GRanges(seqnames = paste0("Chr", featureListTab[[x]]$chrs),
                   ranges = IRanges(start = featureListTab[[x]]$cos-1,
                                    end = featureListTab[[x]]$cos+1),
                   strand = "*")
  featuresListGR <- c(featuresListGR, COsGR)
}

# Specify locations of normalised per base coverage files
libPath <- system(paste0("ls ", covDatPath), intern = T)

# Import coverage files as GRanges objects and assign to library names
covGR <- readGeneric(libPath, meta.col = list(coverage = 4))
assign(paste0(libName), covGR)

mclapply(seq_along(featuresListGR), function(x) {
  # Generate GRanges object containing random loci of same number and
  # size distribution as genesGR
  set.seed(856291)
  ranLocGR <- randomizeRegions(featuresListGR[[x]],
                               genome = genome,
                               per.chromosome = TRUE,
                               allow.overlaps = TRUE)
  
  # Define matrix and column mean coverage outfile (mean profiles)
  outDF <- list(paste0(matDir, libName,
                       "_norm_cov_", featureNames[[x]], "_crossovers_smoothed_target_and_", flankName, "_flank_dataframe.txt"),
                paste0(matDir, libName,
                       "_norm_cov_", featureNames[[x]], "_ranLoc_smoothed_target_and_", flankName, "_flank_dataframe.txt"))
  outDFcolMeans <- list(paste0(matDir, libName,
                               "_norm_cov_", featureNames[[x]], "_crossovers_smoothed_target_and_", flankName, "_flank_dataframe_colMeans.txt"),
                        paste0(matDir, libName,
                               "_norm_cov_", featureNames[[x]], "_ranLoc_smoothed_target_and_", flankName, "_flank_dataframe_colMeans.txt"))
  
  # Run covMatrix() function on each coverage GRanges object to obtain matrices
  ## containing normalised coverage values around target and random loci
  covMatrix(signal = covGR,
            feature = featuresListGR[[x]],
            ranLoc = ranLocGR,
            featureSize = mean(width(featuresListGR[[x]])),
            flankSize = flankSize,
            winSize = winSize,
            outDF = outDF,
            outDFcolMeans = outDFcolMeans)
  print(paste0(libName, " profile calculation complete"))
}, mc.cores = length(featuresListGR))
