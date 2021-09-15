#!/applications/R/R-3.4.0/bin/Rscript

# Plot heatmaps of features sorted by coverage levels between start and end sites

# Usage:
# /applications/R/R-3.4.0/bin/Rscript features_heatmap_sorted.R SPO11_1_oligos_RPI1 COs 2000 2kb '2 kb' 20 20bp bodies

#libName <- "SPO11_1_oligos_RPI1"
#featureName <- "COs"
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 20
#binName <- "20bp"
#region <- "bodies"

args <- commandArgs(trailingOnly = T)
libName <- args[1]
featureName <- args[2]
upstream <- as.numeric(args[3])
downstream <- as.numeric(args[3])
flankName <- args[4]
flankNamePlot <- args[5]
binSize <- as.numeric(args[6])
binName <- args[7]
region <- args[8]

library(EnrichedHeatmap)
library(parallel)
library(circlize)
library(RColorBrewer)
library(tidyr)
library(ggplot2)
library(doParallel)

matDir <- "./matrices/"
plotDir <- "./plots/"
regionPlotDir <- paste0(plotDir, region,
                        "_by_", libName, "/")
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))
system(paste0("[ -d ", regionPlotDir, " ] || mkdir ", regionPlotDir))
plotDir <- regionPlotDir

ChIPseqNames <- c( 
                  "REC8_MYC_Rep1",
                  "REC8_HA_Rep2",
                  "REC8_HA_Rep1",
                  "kss_REC8_HA_Rep1",
                  "kss_REC8_HA_Rep2",
                  "control_MYC_Rep1",
                  "control_HA_Rep1",
                  "MNase",
                  "SPO11_1_ChIP4",
                  "SPO11_1_ChIP13",
                  "H3K9me2",
                  "kss_H3K9me2",
                  "cmt3_H3K9me2",
                  "H3K27me1",
                  "H3K27me3",
                  "H3K27me3_SRR1509478",
                  "H2AW",
                  "H2A",
                  "H2AX",
                  "H2AZ",
                  "H3K4me1",
                  "H3K4me2",
                  "H3K4me3_ChIP14",
                  "H3K4me3_ChIP15",
                  "SPO11_1_oligos_RPI1",
                  "SPO11_1_oligos_RPI3",
                  "SPO11_1_oligos_RPI8",
                  "kss_SPO11_1_oligos_RPI34",
                  "kss_SPO11_1_oligos_RPI35",
                  "PolIV_Rep2",
                  "PolV"
             )
RNAseqNames <- c( 
                 "WT_RNAseq_Chris_Rep1",
                 "WT_RNAseq_Chris_Rep2",
                 "kss_RNAseq_Chris_Rep1",
                 "kss_RNAseq_Chris_Rep2",
                 "WT_RNAseq_Kyuha_Rep1",
                 "WT_RNAseq_Kyuha_Rep2",
                 "WT_RNAseq_Kyuha_Rep3",
                 "WT_RNAseq_meiocyte_Rep1",
                 "WT_RNAseq_meiocyte_Rep2",
                 "WT_RNAseq_meiocyte_Rep3"
                )
methNames <- c(
               "mCG",
               "mCHG",
               "mCHH",
               "kss_mCG",
               "kss_mCHG",
               "kss_mCHH",
               "cmt3_mCG",
               "cmt3_mCHG",
               "cmt3_mCHH"
              )
ChIPseqNamesPlot <- c( 
                      "REC8-Myc Rep1",
                      "REC8-HA Rep2",
                      "REC8-HA Rep1",
                      "kss REC8-HA Rep1",
                      "kss REC8-HA Rep2",
                      "Control-Myc",
                      "Control-HA",
                      "MNase",
                      "SPO11-1 ChIP Rep1",
                      "SPO11-1 ChIP Rep2",
                      "H3K9me2",
                      "kss H3K9me2",
                      "cmt3 H3K9me2",
                      "H3K27me1",
                      "H3K27me3",
                      "H3K27me3 (SRR1509478)",
                      "H2A.W",
                      "H2A",
                      "H2A.X",
                      "H2A.Z",
                      "H3K4me1",
                      "H3K4me2",
                      "H3K4me3 Rep2",
                      "H3K4me3 Rep3",
                      "SPO11-1-oligos Rep1",
                      "SPO11-1-oligos Rep2",
                      "SPO11-1-oligos Rep3",
                      "kss SPO11-1-oligos Rep1",
                      "kss SPO11-1-oligos Rep2",
                      "Pol IV",
                      "Pol V"
                 )
RNAseqNamesPlot <- c(
                     "RNA-seq (floral; Chris) Rep1",
                     "RNA-seq (floral; Chris) Rep2",
                     "kss RNA-seq (floral; Chris) Rep1",
                     "kss RNA-seq (floral; Chris) Rep2",
                     "RNA-seq (floral; Kyuha) Rep1",
                     "RNA-seq (floral; Kyuha) Rep2",
                     "RNA-seq (floral; Kyuha) Rep3",
                     "RNA-seq (meiocyte) Rep1",
                     "RNA-seq (meiocyte) Rep2",
                     "RNA-seq (meiocyte) Rep3"
                    )
methNamesPlot <- c(
                   "mCG",
                   "mCHG",
                   "mCHH",
                   "kss mCG",
                   "kss mCHG",
                   "kss mCHH",
                   "cmt3 mCG",
                   "cmt3 mCHG",
                   "cmt3 mCHH"
                  )

# Load matrix and extract region for sorting of features
mat1 <- as.matrix(read.table(paste0(matDir,
                                    libName,
                                    "_norm_cov_feature_smoothed_target_and_",
                                    flankName, "_flank_dataframe.txt"),
                             header = T))
bodyLength <- (dim(mat1)[2]-((upstream+downstream)/binSize))*binSize
if( region == "promoters" ) {
  mat1Region <- mat1[,(((upstream-500)/binSize)+1):(upstream/binSize)]
} else if ( region == "terminators" ) {
  mat1Region <- mat1[,(((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(500/binSize))]
} else if ( region == "bodies" ) {
  mat1Region <- mat1[,((upstream/binSize)+1):((upstream+bodyLength)/binSize)]
} else {
  print("The region name provided does not match 'promoters', 'terminators', or 'bodies'")
}
mat1RegionRowMeans <- rowMeans(mat1Region, na.rm = T)
mat1RegionRowMeansSorted <- sort.int(mat1RegionRowMeans,
                                     decreasing = T,
                                     index.return = T,
                                     na.last = T)
mat1RegionRowSums <- rowSums(mat1Region, na.rm = T)
mat1RegionRowSumsSorted <- sort.int(mat1RegionRowSums,
                                    decreasing = T,
                                    index.return = T,
                                    na.last = T)
mat1RegionRowMedians <- apply(X = mat1Region,
                              MARGIN = 1,
                              FUN = median)
mat1RegionRowMediansSorted <- sort.int(mat1RegionRowMedians,
                                       decreasing = T,
                                       index.return = T,
                                       na.last = T)

# Load feature coverage matrix for each ChIP-seq dataset
mats <- mclapply(seq_along(ChIPseqNames), function(x) {
  as.matrix(read.table(paste0(matDir,
                              ChIPseqNames[x],                            
                              "_norm_cov_feature_smoothed_target_and_",
                              flankName, "_flank_dataframe.txt"),
                       header = T))
}, mc.cores = length(ChIPseqNames))

# Sort by decreasing mat1RegionRowMeans
matsSorted <- mclapply(seq_along(mats), function(x) {
  mats[[x]][sort.int(mat1RegionRowMeans,
                     decreasing = T,
                     index.return = T,
                     na.last = T)$ix,]
}, mc.cores = length(ChIPseqNames))

# Convert matrices to class "normalizedMatrix" with associated attributes
# for use by EnrichedHeatmap function
for(x in seq_along(matsSorted)) {
  attr(matsSorted[[x]], "upstream_index")         = 1:(upstream/binSize)
  attr(matsSorted[[x]], "target_index")           = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
  attr(matsSorted[[x]], "downstream_index")       = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
  attr(matsSorted[[x]], "extend")                 = c(upstream, downstream)
  attr(matsSorted[[x]], "smooth")                 = TRUE
  attr(matsSorted[[x]], "signal_name")            = ChIPseqNamesPlot[x]
  attr(matsSorted[[x]], "target_name")            = featureName
  attr(matsSorted[[x]], "target_is_single_point") = FALSE
  attr(matsSorted[[x]], "background")             = 0
  attr(matsSorted[[x]], "signal_is_categorical")  = FALSE
  class(matsSorted[[x]]) = c("normalizedMatrix", "matrix")
}

# Load feature coverage matrix for each RNA-seq dataset
RNAseq <- mclapply(seq_along(RNAseqNames), function(x) {
  as.matrix(read.table(paste0(matDir,
                              RNAseqNames[x],                            
                              "_norm_cov_feature_smoothed_target_and_",
                              flankName, "_flank_dataframe.txt"),
                       header = T))
}, mc.cores = length(RNAseqNames))

# Sort by decreasing mat1RegionRowMeans
RNAseqSorted <- mclapply(seq_along(RNAseq), function(x) {
  RNAseq[[x]][sort.int(mat1RegionRowMeans,
                       decreasing = T,
                       index.return = T,
                       na.last = T)$ix,]
}, mc.cores = length(RNAseqNames))

# Convert matrices to class "normalizedMatrix" with associated attributes
# for use by EnrichedHeatmap function
for(x in seq_along(RNAseqSorted)) {
  attr(RNAseqSorted[[x]], "upstream_index")         = 1:(upstream/binSize)
  attr(RNAseqSorted[[x]], "target_index")           = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
  attr(RNAseqSorted[[x]], "downstream_index")       = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
  attr(RNAseqSorted[[x]], "extend")                 = c(upstream, downstream)
  attr(RNAseqSorted[[x]], "smooth")                 = TRUE
  attr(RNAseqSorted[[x]], "signal_name")            = RNAseqNamesPlot[x]
  attr(RNAseqSorted[[x]], "target_name")            = featureName
  attr(RNAseqSorted[[x]], "target_is_single_point") = FALSE
  attr(RNAseqSorted[[x]], "background")             = 0
  attr(RNAseqSorted[[x]], "signal_is_categorical")  = FALSE
  class(RNAseqSorted[[x]]) = c("normalizedMatrix", "matrix")
}

# Load feature DNA methylation matrix for each context
meth <- mclapply(seq_along(methNames), function(x) {
  as.matrix(read.table(paste0(matDir,
                              methNames[x],                            
                              "_norm_cov_feature_smoothed_target_and_",
                              flankName, "_flank_dataframe.txt"),
                       header = T))
}, mc.cores = length(methNames))

# Sort by decreasing mat1RegionRowMeans
methSorted <- mclapply(seq_along(meth), function(x) {
  meth[[x]][sort.int(mat1RegionRowMeans,
                     decreasing = T,
                     index.return = T,
                     na.last = T)$ix,]
}, mc.cores = length(methNames))

# Convert matrices to class "normalizedMatrix" with associated attributes
# for use by EnrichedHeatmap function
for(x in seq_along(methSorted)) {
  attr(methSorted[[x]], "upstream_index")         = 1:(upstream/binSize)
  attr(methSorted[[x]], "target_index")           = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
  attr(methSorted[[x]], "downstream_index")       = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
  attr(methSorted[[x]], "extend")                 = c(upstream, downstream)
  attr(methSorted[[x]], "smooth")                 = TRUE
  attr(methSorted[[x]], "signal_name")            = methNamesPlot[x]
  attr(methSorted[[x]], "target_name")            = featureName
  attr(methSorted[[x]], "target_is_single_point") = FALSE
  attr(methSorted[[x]], "background")             = 0
  attr(methSorted[[x]], "signal_is_categorical")  = FALSE
  class(methSorted[[x]]) = c("normalizedMatrix", "matrix")
}

if(featureName == "genes") {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

# Heatmap plotting function
featureHeatmap <- function(matSorted,
                           col_fun,
                           colour,
                           datName,
                           rowOrder) {
  print(EnrichedHeatmap(mat = matSorted,
                        col = col_fun,
                        top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = colour,
                                                                                              lwd = 3),
                                                                                    yaxis_side = "right",
                                                                                    yaxis_facing = "right",
                                                                                    yaxis_gp = gpar(fontsize = 10),
                                                                                    pos_line_gp = gpar(col = "black",
                                                                                                       lty = 2,
                                                                                                       lwd = 2))),
                        top_annotation_height = unit(2, "cm"),
                        name = "Coverage",
                        row_order = rowOrder,
                        column_title = datName,
                        axis_name = c(paste0("-", flankNamePlot),
                                      featureStartLab, featureEndLab,
                                      paste0("+", flankNamePlot)),
                        axis_name_gp = gpar(fontsize = 12),
                        border = FALSE,
                        pos_line_gp = gpar(col = "black", lty = 2, lwd = 2),
                        use_raster = FALSE))
}

# Define heatmap colours
rich8to6equal <- c("#000041", "#0000CB", "#0081FF", "#FDEE02", "#FFAB00", "#FF3300")

library(doParallel)
registerDoParallel(cores = length(ChIPseqNames))
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

# Plot
# ChIP-seq
foreach(x = 1:length(ChIPseqNames)) %dopar% {
  ChIPseq_col_fun <- colorRamp2(c(-1.0, -0.5, 0.0, 0.5 , 1.0, 1.5),
                                rich8to6equal)
  pdf(paste0(plotDir, ChIPseqNames[x], "_around_", featureName,
             "_heatmap_ordered_by_", libName, "_in_", region, ".pdf"),
      width = 5, height = 10)
  featureHeatmap(matSorted = matsSorted[[x]],
                 col_fun = ChIPseq_col_fun,
                 colour = "red",
                 datName = ChIPseqNamesPlot[x],
                 rowOrder = c(1:dim(matsSorted[[x]])[1]))
  dev.off()
}
# RNA-seq
foreach(x = 1:length(RNAseqNames)) %dopar% {
  RNAseq_col_fun <- colorRamp2(quantile(RNAseqSorted[[x]],
                                        c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
                                        na.rm = T),
                               rich8to6equal)
  pdf(paste0(plotDir, RNAseqNames[x], "_around_", featureName,
             "_heatmap_ordered_by_", libName, "_in_", region, ".pdf"),
      width = 5, height = 10)
  featureHeatmap(matSorted = RNAseqSorted[[x]],
                 col_fun = RNAseq_col_fun,
                 colour = "red",
                 datName = RNAseqNamesPlot[x],
                 rowOrder = c(1:dim(RNAseqSorted[[x]])[1]))
  dev.off()
}
# DNA methylation
foreach(x = 1:length(methNames)) %dopar% {
  meth_col_fun <- colorRamp2(c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                             rich8to6equal)
  pdf(paste0(plotDir, methNames[x], "_around_", featureName,
             "_heatmap_ordered_by_", libName, "_in_", region, ".pdf"),
      width = 5, height = 10)
  featureHeatmap(matSorted = methSorted[[x]],
                 col_fun = meth_col_fun,
                 colour = "red",
                 datName = methNamesPlot[x],
                 rowOrder = c(1:dim(methSorted[[x]])[1]))
  dev.off()
}
