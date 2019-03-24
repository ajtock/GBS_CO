#!/applications/R/R-3.3.2/bin/Rscript

# Profile mean coverage around COs and random loci

# Usage via Condor submission system on node7:
# csmit -m 50G -c 1 "/applications/R/R-3.3.2/bin/Rscript crossover_Profiles_commandArgs.R /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed SPO11_1_oligos_RPI1 5000 5kb"

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc.R")

library(EnrichedHeatmap)
library(genomation)
library(regioneR)

args <- commandArgs(trailingOnly = T)
covDatPath <- as.character(args[1])
libName <- as.character(args[2])
flankSize <- as.numeric(args[3])
flankName <- as.character(args[4])

matDir <- "/projects/ajt200/GBS_CO/BethRowan_190219/BR_wt_ColLer_F2/GBScrossoverProfiles/matrices/"

# Chromosome definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))

# Import GBS COs as GRanges object
COs <- read.csv("./CO.all.with.420.flagged.csv")
COsGR <- GRanges(seqnames = paste0("Chr", COs$chr),
                 ranges = IRanges(start = COs$breakpoint.pos-1,
                                  end = COs$breakpoint.pos+1),
                 strand = "*")
print("***********COs***********")
print(length(COsGR))

# Generate GRanges object containing random loci of same number and
# size distribution as COsGR
ranLocGR <- randomizeRegions(COsGR,
                             genome = genome,
                             per.chromosome = TRUE,
                             allow.overlaps = TRUE)

# Specify locations of normalised per base coverage files
libPath <- system(paste0("ls ", covDatPath), intern = T)

# Import coverage files as GRanges objects and assign to library names
gr <- readGeneric(libPath, meta.col = list(coverage = 4))
assign(paste0(libName), gr)

# Define column mean coverage outfile (mean profiles)
outDFCM <- list(list(paste0(matDir, libName,
                            "_norm_cov_COs_mat1_smoothed_midpoint_and_",
                            flankName, "flank_dataframe_colMeans.txt"),
                     paste0(matDir, libName,
                            "_norm_cov_ranLoc_mat2_smoothed_midpoint_and_",
                            flankName, "flank_dataframe_colMeans.txt")))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
covMatrix(signal = gr, target = COsGR, ranLoc = ranLocGR,
          targetSize = mean(width(COsGR)), flankSize = flankSize, winSize = 20,
          outDFCM = outDFCM, x = 1)
print(paste0(libName, " profile calculation complete"))
