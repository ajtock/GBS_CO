#!/applications/R/R-3.3.2/bin/Rscript

# Profile mean coverage around COs and random loci

# Usage via Condor submission system on node7:
# csmit -m 50G -c 1 "/applications/R/R-3.3.2/bin/Rscript crossover_Profiles_commandArgs.R 5000 5kb 20 /projects/ajt200/BAM_masters/RNAseq_meiocyte_Feng/snakemake_RNAseq_STAR/mapped/both/pb/Col_RNAseq_meiocyte_Rep1_MappedOn_TAIR10_chr_all_both_sort_norm.perbase Col_RNAseq_meiocyte_Rep1"

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc.R")

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
 
# Import GBS wt_COs as GRanges object
wt_COs <- load("/projects/ajt200/GBS_CO/All_081018/coordinates/pooled_lanes/EL_wt_ColBur_F2_2019_pooled_GRanges.RData")
wt_COsGR <- popCOsGR_pooled
rm(popCOsGR_pooled)
wt_COsGR <- GRanges(seqnames = seqnames(wt_COsGR),
                    ranges = IRanges(start = wt_COsGR$midpoint-1,
                                     end = wt_COsGR$midpoint+1),
                    strand = "*")
print("***********wt_COs***********")
print(length(wt_COsGR))

# Generate GRanges object containing random loci of same number and
# size distribution as wt_COsGR
set.seed(753492)
wt_ranLocGR <- randomizeRegions(wt_COsGR,
                                genome = genome,
                                per.chromosome = TRUE,
                                allow.overlaps = TRUE)

# Import GBS taf4b_COs as GRanges object
taf4b_COs <- load("/projects/ajt200/GBS_CO/All_081018/coordinates/pooled_lanes/EL_taf4b_ColBur_F2_2019_pooled_GRanges.RData")
taf4b_COsGR <- popCOsGR_pooled
rm(popCOsGR_pooled)
taf4b_COsGR <- GRanges(seqnames = seqnames(taf4b_COsGR),
                       ranges = IRanges(start = taf4b_COsGR$midpoint-1,
                                        end = taf4b_COsGR$midpoint+1),
                       strand = "*")
print("***********taf4b_COs***********")
print(length(taf4b_COsGR))

# Generate GRanges object containing random loci of same number and
# size distribution as taf4b_COsGR
set.seed(639572)
taf4b_ranLocGR <- randomizeRegions(taf4b_COsGR,
                                   genome = genome,
                                   per.chromosome = TRUE,
                                   allow.overlaps = TRUE)

# Specify locations of normalised per base coverage files
libPath <- system(paste0("ls ", covDatPath), intern = T)

# Import coverage files as GRanges objects and assign to library names
cov <- read.table(libPath, header = F)
cov <- cov[cov[,1] != "mitochondria" & cov[,1] != "chloroplast",]
gr <- GRanges(seqnames = paste0("Chr", cov[,1]),
              ranges = IRanges(start = cov[,2],
                               end = cov[,2]),
              strand = "*",
              coverage = cov[,3])
assign(paste0(libName), gr)

# Define column mean coverage outfile (mean profiles)
wt_outDFCM <- list(list(paste0(matDir, libName,
                               "_norm_cov_wt_COs_mat1_smoothed_midpoint_and_",
                               flankName, "flank_dataframe_colMeans.txt"),
                        paste0(matDir, libName,
                               "_norm_cov_wt_ranLoc_mat2_smoothed_midpoint_and_",
                               flankName, "flank_dataframe_colMeans.txt")))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
covMatrix(signal = gr, target = wt_COsGR, ranLoc = wt_ranLocGR,
          targetSize = mean(width(wt_COsGR)), flankSize = flankSize, winSize = winSize,
          outDFCM = wt_outDFCM, x = 1)

# Define column mean coverage outfile (mean profiles)
taf4b_outDFCM <- list(list(paste0(matDir, libName,
                                  "_norm_cov_taf4b_COs_mat1_smoothed_midpoint_and_",
                                  flankName, "flank_dataframe_colMeans.txt"),
                           paste0(matDir, libName,
                                  "_norm_cov_taf4b_ranLoc_mat2_smoothed_midpoint_and_",
                                  flankName, "flank_dataframe_colMeans.txt")))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
covMatrix(signal = gr, target = taf4b_COsGR, ranLoc = taf4b_ranLocGR,
          targetSize = mean(width(taf4b_COsGR)), flankSize = flankSize, winSize = winSize,
          outDFCM = taf4b_outDFCM, x = 1)

print(paste0(libName, " profile calculation complete"))
