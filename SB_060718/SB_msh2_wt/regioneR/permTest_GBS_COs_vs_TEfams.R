#!/applications/R/R-3.3.2/bin/Rscript

# Use permutation test function in regioneR to determine if crossover intervals overlap features
# of interest (e.g., other peaks, genes, TEs) more or less than expected by chance

# Usage on hydrogen node7:
# csmit -m 50G -c 48 "/applications/R/R-3.3.2/bin/Rscript permTest_GBS_COs_vs_TEfams.R SB_msh2_ColCvi_F2_2018 SB_ColCvi_F2_2018_maskGR.RData"

library(regioneR)

args <- commandArgs(trailingOnly = T)
popName <- args[1]
maskFile <- args[2]

inDir <- "/projects/ajt200/GBS_CO/All_081018/coordinates/pooled_lanes/"
outDir <- "/projects/ajt200/GBS_CO/SB_060718/SB_msh2_wt/regioneR/"
DNAplotDir <- "/projects/ajt200/GBS_CO/SB_060718/SB_msh2_wt/regioneR/plots/hist/TEsDNA/"
RNAplotDir <- "/projects/ajt200/GBS_CO/SB_060718/SB_msh2_wt/regioneR/plots/hist/TEsRNA/"

# Genomic definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
load(maskFile)
mask <- maskGR
print(paste0(maskFile, " masked regions:"))
print(mask)

# Load GBS crossover intervals
load(paste0(inDir, popName, "_pooled_GRanges.RData"))
print(popCOsGR)
popCOsGR_maskOverlaps <- findOverlaps(mask, popCOsGR, ignore.strand = T, select = "all")
print("Crossover intervals within masked regions:")
print(popCOsGR_maskOverlaps)

DNAfamNames <- c("dna", "heli", "ptmari", "mudr", "enspm", "hat", "harbinger")
RNAfamNames <- c("rna", "gypsy", "copia", "linel1", "sine")
DNAdir <- "/projects/ajt200/TAIR10/TE_classes/DNA/"
RNAdir <- "/projects/ajt200/TAIR10/TE_classes/RNA/"

### DNA TEs
TEsDNAGR <- lapply(seq_along(DNAfamNames), function(x) {
  TEsDNA <- read.table(paste0(DNAdir,
                              "TAIR10_Buisine_TEs_strand_tab_ann_",
                              DNAfamNames[x], ".txt"), header = T)
  GRanges(seqnames = TEsDNA$chr,
          ranges = IRanges(start = TEsDNA$start,
                           end = TEsDNA$end),
          strand = "*")
})

# Remove TEs located within pericentromeric regions
TEsDNAGRmasked <- lapply(seq_along(TEsDNAGR), function(x) {
  TEsDNAGR_maskOverlaps <- findOverlaps(mask, TEsDNAGR[[x]], ignore.strand = TRUE, select = "all")
  TEsDNAGR[[x]][-subjectHits(TEsDNAGR_maskOverlaps)]
})
 
# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions in B as in A
set.seed(123)
ptCOsTEsDNAPerChrom <- lapply(seq_along(TEsDNAGRmasked), function(x) {
  permTest(A = popCOsGR, B = TEsDNAGRmasked[[x]], genome = genome, mask = mask,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE, per.chromosome = TRUE,
           evaluate.function = numOverlaps, count.once = TRUE,
           ntimes = 10000, mc.set.seed = FALSE, mc.cores = detectCores())
})

for(i in 1:length(ptCOsTEsDNAPerChrom)) {
  assign(paste0(DNAfamNames[i]), ptCOsTEsDNAPerChrom[[i]])
}
save(ptCOsTEsDNAPerChrom,
     file = paste0(outDir,
                   "permTest_", popName, "_COs_vs_TEsDNA.RData"))

# Summarise results in a table
noOfFeatures <- NULL
expected <- NULL
observed <- NULL
pval <- NULL
zscore <- NULL
for(i in 1:length(ptCOsTEsDNAPerChrom)) {
  noOfFeaturesi <- print(length(TEsDNAGRmasked[[i]]))
  noOfFeatures <- c(noOfFeatures, noOfFeaturesi)
  expectedi <- print(round(mean(ptCOsTEsDNAPerChrom[[i]]$numOverlaps$permuted)))
  expected <- c(expected, expectedi)
  observedi <- print(ptCOsTEsDNAPerChrom[[i]]$numOverlaps$observed)
  observed <- c(observed, observedi)
  pvali <- print(round(ptCOsTEsDNAPerChrom[[i]]$numOverlaps$pval, 4))
  pval <- c(pval, pvali)
  zscorei <- print(round(ptCOsTEsDNAPerChrom[[i]]$numOverlaps$zscore, 4))
  zscore <- c(zscore, zscorei)
}
ptCOsTEsDNAPerChromDataFrame <- cbind(noOfFeatures, expected, observed, pval, zscore)
write.table(ptCOsTEsDNAPerChromDataFrame,
            file = paste0(outDir,
                          "permTest_", popName, "_COs_vs_TEsDNA_DataFrame.txt"),
            sep = "\t", row.names = F)

# plot graphical summaries of results
for(i in 1:length(ptCOsTEsDNAPerChrom)) {
  pdf(file = paste0(DNAplotDir, DNAfamNames[i], "_permTest_nperm10000_", popName, "_COs_perChrom.pdf"), width = 10, height = 7)
  plot(ptCOsTEsDNAPerChrom[[i]], main = paste0(popName, " COs vs ", DNAfamNames[i]), xlab = "Number of overlaps", ylab = "Density")
  dev.off()

  # Using the localZScore() function, evaluate whether the association between COs and TEsDNAfam is highly dependent on their exact position
  lz_1kb <- localZScore(pt = ptCOsTEsDNAPerChrom[[i]], A = popCOsGR, B = TEsDNAGRmasked[[i]],
                        window = 1000, step = 50, count.once = TRUE)
  lz_10kb <- localZScore(pt = ptCOsTEsDNAPerChrom[[i]], A = popCOsGR, B = TEsDNAGRmasked[[i]],
                         window = 10000, step = 500, count.once = TRUE)
  lz_custom <- localZScore(pt = ptCOsTEsDNAPerChrom[[i]], A = popCOsGR, B = TEsDNAGRmasked[[i]],
                           window = 10*mean(width(popCOsGR)), step = mean(width(popCOsGR))/2, count.once = TRUE)
  win <- as.character(round((10*mean(width(popCOsGR)))/1000))
  step <- as.character(round(mean(width(popCOsGR))/2))
  pdf(file = paste0(DNAplotDir, DNAfamNames[i], "_localZscore_permTest_nperm10000_", popName, "_COs_w1kb_s50bp_w10kb_s500bp_w", win ,"kb_s", step, "bp_perChrom.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0(popName, " COs vs ", DNAfamNames[i], " (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(popName, " COs vs ", DNAfamNames[i], " (1-kb shift)"))
  plot(lz_10kb, main = paste0(popName, " COs vs ", DNAfamNames[i], " (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(popName, " COs vs ", DNAfamNames[i], " (10-kb shift)"))
  plot(lz_custom, main = paste0(popName, " COs vs ", DNAfamNames[i], " (~", win, "-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(popName, " COs vs ", DNAfamNames[i], " (~", win, "-kb shift)"))
  dev.off()
}

### RNA TEs
TEsRNAGR <- lapply(seq_along(RNAfamNames), function(x) {
  TEsRNA <- read.table(paste0(RNAdir,
                              "TAIR10_Buisine_TEs_strand_tab_ann_",
                              RNAfamNames[x], ".txt"), header = T)
  GRanges(seqnames = TEsRNA$chr,
          ranges = IRanges(start = TEsRNA$start,
                           end = TEsRNA$end),
          strand = "*")
})

# Remove TEs located within pericentromeric regions
TEsRNAGRmasked <- lapply(seq_along(TEsRNAGR), function(x) {
  TEsRNAGR_maskOverlaps <- findOverlaps(mask, TEsRNAGR[[x]], ignore.strand = TRUE, select = "all")
  TEsRNAGR[[x]][-subjectHits(TEsRNAGR_maskOverlaps)]
})
 
# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions in B as in A
set.seed(123)
ptCOsTEsRNAPerChrom <- lapply(seq_along(TEsRNAGRmasked), function(x) {
  permTest(A = popCOsGR, B = TEsRNAGRmasked[[x]], genome = genome, mask = mask,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE, per.chromosome = TRUE,
           evaluate.function = numOverlaps, count.once = TRUE,
           ntimes = 10000, mc.set.seed = FALSE, mc.cores = detectCores())
})

for(i in 1:length(ptCOsTEsRNAPerChrom)) {
  assign(paste0(RNAfamNames[i]), ptCOsTEsRNAPerChrom[[i]])
}
save(ptCOsTEsRNAPerChrom,
     file = paste0(outDir,
                   "permTest_", popName, "_COs_vs_TEsRNA.RData"))

# Summarise results in a table
noOfFeatures <- NULL
expected <- NULL
observed <- NULL
pval <- NULL
zscore <- NULL
for(i in 1:length(ptCOsTEsRNAPerChrom)) {
  noOfFeaturesi <- print(length(TEsRNAGRmasked[[i]]))
  noOfFeatures <- c(noOfFeatures, noOfFeaturesi)
  expectedi <- print(round(mean(ptCOsTEsRNAPerChrom[[i]]$numOverlaps$permuted)))
  expected <- c(expected, expectedi)
  observedi <- print(ptCOsTEsRNAPerChrom[[i]]$numOverlaps$observed)
  observed <- c(observed, observedi)
  pvali <- print(round(ptCOsTEsRNAPerChrom[[i]]$numOverlaps$pval, 4))
  pval <- c(pval, pvali)
  zscorei <- print(round(ptCOsTEsRNAPerChrom[[i]]$numOverlaps$zscore, 4))
  zscore <- c(zscore, zscorei)
}
ptCOsTEsRNAPerChromDataFrame <- cbind(noOfFeatures, expected, observed, pval, zscore)
write.table(ptCOsTEsRNAPerChromDataFrame,
            file = paste0(outDir,
                          "permTest_", popName, "_COs_vs_TEsRNA_DataFrame.txt"),
            sep = "\t", row.names = F)

# plot graphical summaries of results
for(i in 1:length(ptCOsTEsRNAPerChrom)) {
  pdf(file = paste0(RNAplotDir, RNAfamNames[i], "_permTest_nperm10000_", popName, "_COs_perChrom.pdf"), width = 10, height = 7)
  plot(ptCOsTEsRNAPerChrom[[i]], main = paste0(popName, " COs vs ", RNAfamNames[i]), xlab = "Number of overlaps", ylab = "Density")
  dev.off()

  # Using the localZScore() function, evaluate whether the association between COs and TEsRNAfam is highly dependent on their exact position
  lz_1kb <- localZScore(pt = ptCOsTEsRNAPerChrom[[i]], A = popCOsGR, B = TEsRNAGRmasked[[i]],
                        window = 1000, step = 50, count.once = TRUE)
  lz_10kb <- localZScore(pt = ptCOsTEsRNAPerChrom[[i]], A = popCOsGR, B = TEsRNAGRmasked[[i]],
                         window = 10000, step = 500, count.once = TRUE)
  lz_custom <- localZScore(pt = ptCOsTEsRNAPerChrom[[i]], A = popCOsGR, B = TEsRNAGRmasked[[i]],
                           window = 10*mean(width(popCOsGR)), step = mean(width(popCOsGR))/2, count.once = TRUE)
  win <- as.character(round((10*mean(width(popCOsGR)))/1000))
  step <- as.character(round(mean(width(popCOsGR))/2))
  pdf(file = paste0(RNAplotDir, RNAfamNames[i], "_localZscore_permTest_nperm10000_", popName, "_COs_w1kb_s50bp_w10kb_s500bp_w", win ,"kb_s", step, "bp_perChrom.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0(popName, " COs vs ", RNAfamNames[i], " (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(popName, " COs vs ", RNAfamNames[i], " (1-kb shift)"))
  plot(lz_10kb, main = paste0(popName, " COs vs ", RNAfamNames[i], " (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(popName, " COs vs ", RNAfamNames[i], " (10-kb shift)"))
  plot(lz_custom, main = paste0(popName, " COs vs ", RNAfamNames[i], " (~", win, "-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(popName, " COs vs ", RNAfamNames[i], " (~", win, "-kb shift)"))
  dev.off()
}
