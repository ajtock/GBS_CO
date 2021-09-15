#!/applications/R/R-3.3.2/bin/Rscript

# Use permutation test function in regioneR to determine
# if loci overlap features of interest more or less
# than expected by chance

# Usage on hydrogen node7:
# csmit -m 100G -c 32 "/applications/R/R-3.3.2/bin/Rscript permTest_ovule_hypermeth_vs_others.R ovule_hypermeth_loci 10000"

lociName <- "ovule_hypermeth_loci"
perms <- 10000

args <- commandArgs(trailingOnly = T)
lociName <- args[1]
perms <- as.numeric(args[2])

library(regioneR)

# Genomic definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))

outDir <- "./"

plotDir <- "histograms/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Import loci and convert into a GRanges object
loci <- read.table("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/ovule_hyper_DMR_siRNA_Julie_Law/ovule_hyper_siRNA.bed", header = F, stringsAsFactors = F)
lociGR <- GRanges(seqnames = paste0("Chr", loci$V1),
                  ranges = IRanges(start = loci$V2-1,
                                   end = loci$V3),
                  strand = "*",
                  ID = loci$V4)
print("***********Ovule loci***********")
print(lociGR)

# Others
load("/projects/ajt200/GBS_CO/HS_CU_080617/wt/COsGRcoords.RData")
COsGR <- COsGRcoords
print(length(COsGR))
#[1] 3320

RowanCOs <- read.csv("/projects/ajt200/GBS_CO/BethRowan_190219/BR_wt_ColLer_F2/CO.all.with.420.flagged.csv", header = T)
RowanCOsGR <- GRanges(seqnames = paste0("Chr", RowanCOs$chr),
                      ranges = IRanges(start = RowanCOs$block1.end.pos,
                                       end = RowanCOs$block2.start.pos),
                      strand = "*")
print("***********RowanCOs***********")
print(length(RowanCOsGR))
#[1] 17077

load("/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/armrangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData")
load("/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/perirangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData")
SPO11GR <- sort(c(armrangerPeaksGRmerge, perirangerPeaksGRmerge))
armrangerPeaksGRmerge <- NULL
perirangerPeaksGRmerge <- NULL
strand(SPO11GR) <- "*"
print(length(SPO11GR))
#[1] 5914

load("/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/armPeaksSH99GRmerge.RData")
load("/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/periPeaksSH99GRmerge.RData")
nucleRnucsGR <- sort(c(armPeaksSH99GRmerge, periPeaksSH99GRmerge))
armPeaksSH99GRmerge <- NULL
periPeaksSH99GRmerge <- NULL
strand(nucleRnucsGR) <- "*"
print(length(nucleRnucsGR))
#[1] 57734

otherNames <- c(
                "COsGR",
                "RowanCOsGR",
                "SPO11GR",
                "nucleRnucsGR"
               )
grl <- c(
         "COsGR" = COsGR,
         "RowanCOsGR" = RowanCOsGR,
         "SPO11GR" = SPO11GR,
         "nucleRnucsGR" = nucleRnucsGR
        )

# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions as in A
set.seed(38402)
ptLociOtherPerChrom <- lapply(seq_along(grl), function(x) {
  permTest(A = lociGR,
           B = grl[[x]],
           genome = genome,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE,
           per.chromosome = TRUE,
           evaluate.function = numOverlaps,
           count.once = TRUE,
           ntimes = perms,
           mc.set.seed = FALSE,
           mc.cores = detectCores())
})

for(i in 1:length(ptLociOtherPerChrom)) {
  assign(paste0(otherNames[i]), ptLociOtherPerChrom[[i]])
}
save(ptLociOtherPerChrom,
     file = paste0(outDir,
                   "permTest_", as.character(perms), "perms_", lociName, "_vs_others.RData"))

# Summarise results in a table
featureName <- NULL
noOfFeatures <- NULL
expected <- NULL
observed <- NULL
pval <- NULL
zscore <- NULL
for(i in 1:length(ptLociOtherPerChrom)) {
  featureNamei <- print(otherNames[i])
  featureName <- c(featureName, featureNamei)
  noOfFeaturesi <- print(length(grl[[i]]))
  noOfFeatures <- c(noOfFeatures, noOfFeaturesi)
  expectedi <- print(round(mean(ptLociOtherPerChrom[[i]]$numOverlaps$permuted)))
  expected <- c(expected, expectedi)
  observedi <- print(ptLociOtherPerChrom[[i]]$numOverlaps$observed)
  observed <- c(observed, observedi)
  pvali <- print(round(ptLociOtherPerChrom[[i]]$numOverlaps$pval, 4))
  pval <- c(pval, pvali)
  zscorei <- print(round(ptLociOtherPerChrom[[i]]$numOverlaps$zscore, 4))
  zscore <- c(zscore, zscorei)
}
ptLociOtherPerChromDataFrame <- cbind(featureName, noOfFeatures, expected, observed, pval, zscore)
colnames(ptLociOtherPerChromDataFrame) <- c("feature", "n", "expected", "observed", "pval", "zscore")
write.table(ptLociOtherPerChromDataFrame,
            file = paste0(outDir,
                          "permTest_", as.character(perms), "perms_", lociName, "_vs_others_DataFrame.txt"),
            sep = "\t", quote = F, row.names = F)

# plot graphical summaries of results
for(i in 1:length(ptLociOtherPerChrom)) {
  pdf(file = paste0(plotDir, otherNames[i], "_permTest_", as.character(perms), "perms_", lociName, "_perChrom.pdf"), width = 10, height = 7)
  plot(ptLociOtherPerChrom[[i]], main = paste0(lociName, " vs ", otherNames[i]), xlab = "Number of overlaps", ylab = "Density")
  dev.off()

  # Using the localZScore() function, evaluate whether the association between loci and other is highly dependent on their exact position
  lz_1kb <- localZScore(pt = ptLociOtherPerChrom[[i]], A = lociGR, B = grl[[i]],
                        window = 1000, step = 50, count.once = TRUE)
  lz_10kb <- localZScore(pt = ptLociOtherPerChrom[[i]], A = lociGR, B = grl[[i]],
                         window = 10000, step = 500, count.once = TRUE)
  lz_custom <- localZScore(pt = ptLociOtherPerChrom[[i]], A = lociGR, B = grl[[i]],
                           window = 10*mean(width(lociGR)), step = mean(width(lociGR))/2, count.once = TRUE)
  win <- as.character(round((10*mean(width(lociGR)))/1000))
  step <- as.character(round(mean(width(lociGR))/2))
  pdf(file = paste0(plotDir, otherNames[i], "_localZscore_permTest_", as.character(perms), "perms_", lociName, "_w1kb_s50bp_w10kb_s500bp_w", win ,"kb_s", step, "bp_perChrom.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0(lociName, " vs ", otherNames[i], " (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(lociName, " vs ", otherNames[i], " (1-kb shift)"))
  plot(lz_10kb, main = paste0(lociName, " vs ", otherNames[i], " (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(lociName, " vs ", otherNames[i], " (10-kb shift)"))
  plot(lz_custom, main = paste0(lociName, " vs ", otherNames[i], " (~", win, "-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(lociName, " vs ", otherNames[i], " (~", win, "-kb shift)"))
  dev.off()
}
