#!/applications/R/R-3.3.2/bin/Rscript

# Based on histograms of COs per individual generated with hist_COs_per_ind.R ,
# remove individuals that are obvious outliers
# (e.g., 1 individual with 24 COs in CL_asy1het_ColLer_F2_2018)

# Usage:
# ./drop_outliers.R CL_asy1het_ColLer_F2_2018

args <- commandArgs(trailingOnly = T)
popName <- "CL_asy1het_ColLer_F2_2018"
popName <- args[1]

library(GenomicRanges)

inDir <- "/projects/ajt200/GBS_CO/All_081018/coordinates/pooled_lanes/"

load(paste0(inDir, popName, "_pooled_GRanges.RData"))
COs_per_ind <- NULL
for(j in 1:length(unique(popCOsGR_pooled$lib))) {
   ind_COsGR <- popCOsGR_pooled[popCOsGR_pooled$lib == unique(popCOsGR_pooled$lib)[j]]
   COs_per_ind <- c(COs_per_ind, length(ind_COsGR))
}
quant <- quantile(COs_per_ind, 0.999)[[1]]
outlier_lib <- unique(popCOsGR_pooled$lib)[COs_per_ind > quant]
outlier_COs <- popCOsGR_pooled[popCOsGR_pooled$lib == outlier_lib[1]]
popCOsGR_pooled <- popCOsGR_pooled[popCOsGR_pooled$lib != outlier_lib[1]]
system(paste0("mv ", inDir, popName, "_pooled_GRanges.RData ", inDir, popName, "_pooled_GRanges_incl_outlier_lib2.82.RData"))
save(popCOsGR_pooled,
     file = paste0(inDir, popName, "_pooled_GRanges.RData"))

