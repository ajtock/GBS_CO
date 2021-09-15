#!/applications/R/R-3.3.2/bin/Rscript

# Use permutation test function in regioneR to determine if crossover intervals overlap features
# of interest (e.g., other peaks, genes, TEs) more or less than expected by chance

# Usage on hydrogen node7:
# csmit -m 50G -c 48 "/applications/R/R-3.3.2/bin/Rscript permTest_GBS_COs_vs_others.R SB_msh2_ColCvi_F2_2018 SB_ColCvi_F2_2018_maskGR.RData"

library(regioneR)

args <- commandArgs(trailingOnly = T)
popName <- args[1]
maskFile <- args[2]

inDir <- "/projects/ajt200/GBS_CO/All_081018/coordinates/pooled_lanes/"
outDir <- "/projects/ajt200/GBS_CO/SB_060718/SB_msh2_wt/regioneR/"
plotDir <- "/projects/ajt200/GBS_CO/SB_060718/SB_msh2_wt/regioneR/plots/hist/"

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
                             
# Others
# MSH4_Rep1 peaks
load(paste0("/home/ajt200/analysis/160902_Sasha_ChIP_MSH4_Rep1/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
            "MSH4_Rep1_armrangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
load(paste0("/home/ajt200/analysis/160902_Sasha_ChIP_MSH4_Rep1/peaks/PeakRanger1.18/ranger/p0.001_q0.01/",
            "MSH4_Rep1_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
MSH4_Rep1GR <- sort(c(armrangerPeaksGRmergedOverlaps, perirangerPeaksGRmergedOverlaps))
armrangerPeaksGRmergedOverlaps <- NULL
perirangerPeaksGRmergedOverlaps <- NULL
strand(MSH4_Rep1GR) <- "*"
print("***********peaks***********")
print(MSH4_Rep1GR)
print(length(MSH4_Rep1GR))
#[1] 40556
MSH4_Rep1GR_maskOvelaps <- findOverlaps(mask, MSH4_Rep1GR, ignore.strand = TRUE, select = "all")
MSH4_Rep1GR <- MSH4_Rep1GR[-subjectHits(MSH4_Rep1GR_maskOvelaps)]
print(length(MSH4_Rep1GR))
#[1] 36327

# REC8_HA_Rep1 peaks
load(paste0("/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/",
            "REC8_HA_Rep1_armrangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
load(paste0("/home/ajt200/analysis/REC8_pooled/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.001_q0.01/",
            "REC8_HA_Rep1_perirangerPeaksGRmergedOverlaps_minuslog10_p0.001_q0.01_noMinWidth.RData"))
REC8_HA_Rep1GR <- sort(c(armrangerPeaksGRmergedOverlaps, perirangerPeaksGRmergedOverlaps))
armrangerPeaksGRmergedOverlaps <- NULL
perirangerPeaksGRmergedOverlaps <- NULL
strand(REC8_HA_Rep1GR) <- "*"
print("***********peaks***********")
print(REC8_HA_Rep1GR)
print(length(REC8_HA_Rep1GR))
#[1] 87738
REC8_HA_Rep1GR_maskOvelaps <- findOverlaps(mask, REC8_HA_Rep1GR, ignore.strand = TRUE, select = "all")
REC8_HA_Rep1GR <- REC8_HA_Rep1GR[-subjectHits(REC8_HA_Rep1GR_maskOvelaps)]
print(length(REC8_HA_Rep1GR))
#[1] 79405

# Well-positioned nucleosomes (nucleR)
load("/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/armPeaksSH99GRmerge.RData")
load("/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/periPeaksSH99GRmerge.RData")
nucleRnucsGR <- sort(c(armPeaksSH99GRmerge, periPeaksSH99GRmerge))
armPeaksSH99GRmerge <- NULL
periPeaksSH99GRmerge <- NULL
strand(nucleRnucsGR) <- "*"
print(length(nucleRnucsGR))
#[1] 57734
nucleRnucsGR_maskOvelaps <- findOverlaps(mask, nucleRnucsGR, ignore.strand = TRUE, select = "all")
nucleRnucsGR <- nucleRnucsGR[-subjectHits(nucleRnucsGR_maskOvelaps)]
print(length(nucleRnucsGR))
#[1] 52531

# SPO11-1-oligo hotspots
load("/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/armrangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData")
load("/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/perirangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData")
SPO11GR <- sort(c(armrangerPeaksGRmerge, perirangerPeaksGRmerge))
armrangerPeaksGRmerge <- NULL
perirangerPeaksGRmerge <- NULL
strand(SPO11GR) <- "*"
print(length(SPO11GR))
#[1] 5914
SPO11GR_maskOvelaps <- findOverlaps(mask, SPO11GR, ignore.strand = TRUE, select = "all")
SPO11GR <- SPO11GR[-subjectHits(SPO11GR_maskOvelaps)]
print(length(SPO11GR))
#[1] 5382

# SPO11-1 ChIP4 peaks
load("/projects/ajt200/BAM_masters/SPO11_ChIP/Xiaohui_BAM_and_coverage_files/WT/peaks/PeakRanger1.18/ranger/REC8_MYC_Rep2_input_p0.05_q0.05/armrangerPeaksGRmerge_WT_SPO11_ChIP4_p0.05_q0.05_noMinWidth.RData")
load("/projects/ajt200/BAM_masters/SPO11_ChIP/Xiaohui_BAM_and_coverage_files/WT/peaks/PeakRanger1.18/ranger/REC8_MYC_Rep2_input_p0.05_q0.05/perirangerPeaksGRmerge_WT_SPO11_ChIP4_p0.05_q0.05_noMinWidth.RData")
SPO11_ChIP4GR <- sort(c(armrangerPeaksGRmerge, perirangerPeaksGRmerge))
armrangerPeaksGRmerge <- NULL
perirangerPeaksGRmerge <- NULL
strand(SPO11_ChIP4GR) <- "*"
print(length(SPO11_ChIP4GR))
#[1] 25825
SPO11_ChIP4GR_maskOvelaps <- findOverlaps(mask, SPO11_ChIP4GR, ignore.strand = TRUE, select = "all")
SPO11_ChIP4GR <- SPO11_ChIP4GR[-subjectHits(SPO11_ChIP4GR_maskOvelaps)]
print(length(SPO11_ChIP4GR))
#[1] 23550

# H3K4me3 peaks
load("/projects/ajt200/BAM_masters/H3K4me3/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/armrangerPeaksGRmerge_WT_H3K4me3_ChIP14_WT_H3K4me3_ChIP15_idr0.05_noMinWidth.RData")
load("/projects/ajt200/BAM_masters/H3K4me3/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/perirangerPeaksGRmerge_WT_H3K4me3_ChIP14_WT_H3K4me3_ChIP15_idr0.05_noMinWidth.RData")
H3K4me3GR <- sort(c(armrangerPeaksGRmerge, perirangerPeaksGRmerge))
armrangerPeaksGRmerge <- NULL
perirangerPeaksGRmerge <- NULL
strand(H3K4me3GR) <- "*"
print(length(H3K4me3GR))
#[1] 13951
H3K4me3GR_maskOvelaps <- findOverlaps(mask, H3K4me3GR, ignore.strand = TRUE, select = "all")
H3K4me3GR <- H3K4me3GR[-subjectHits(H3K4me3GR_maskOvelaps)]
print(length(H3K4me3GR))
#[1] 12642

# H3K9me2 peaks (ranger)
load("/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/ranger/p0.05_q0.05/WT_H3K9me2_ChIP_armrangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth.RData")
load("/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/ranger/p0.05_q0.05/WT_H3K9me2_ChIP_perirangerPeaksGRmergedOverlaps_minuslog10_p0.05_q0.05_noMinWidth.RData")
H3K9me2GR <- sort(c(armrangerPeaksGRmergedOverlaps, perirangerPeaksGRmergedOverlaps))
armrangerPeaksGRmergedOverlaps <- NULL
perirangerPeaksGRmergedOverlaps <- NULL 
strand(H3K9me2GR) <- "*"
print(length(H3K9me2GR))
#[1] 20289
H3K9me2GR_maskOvelaps <- findOverlaps(mask, H3K9me2GR, ignore.strand = TRUE, select = "all")
H3K9me2GR <- H3K9me2GR[-subjectHits(H3K9me2GR_maskOvelaps)]
print(length(H3K9me2GR))
#[1] 19443

# Genes
genes <- read.table(file = "/projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt", header = T)
genesGR <- GRanges(seqnames = genes$chr,
                   ranges = IRanges(start = genes$start,
                                    end = genes$end),
                   strand = "*")
seqlevels(genesGR) <- sub("", "Chr", seqlevels(genesGR))
print(length(genesGR))
#[1] 27204
genesGR_maskOvelaps <- findOverlaps(mask, genesGR, ignore.strand = TRUE, select = "all")
genesGR <- genesGR[-subjectHits(genesGR_maskOvelaps)]
print(length(genesGR))
#[1] 24547

# Gene promoters (TSS-500 bp to TSS-1 bp)
genesGRprom <- GRanges(seqnames = genes$chr,
                       ranges = IRanges(start = genes$start,
                                        end = genes$end),
                       strand = genes$strand)
promotersGR <- promoters(genesGRprom, upstream = 500, downstream = 0)
seqlevels(promotersGR) <- sub("", "Chr", seqlevels(promotersGR))
strand(promotersGR) <- "*"
print(length(promotersGR))
#[1] 27204
promotersGR_maskOvelaps <- findOverlaps(mask, promotersGR, ignore.strand = TRUE, select = "all")
promotersGR <- promotersGR[-subjectHits(promotersGR_maskOvelaps)]
print(length(promotersGR))
#[1] 24553

# Gene 5´ ends (TSS to +500 bp)
TSSdownstream500GR <- promoters(genesGRprom, upstream = 0, downstream = 500)
seqlevels(TSSdownstream500GR) <- sub("", "Chr", seqlevels(TSSdownstream500GR))
strand(TSSdownstream500GR) <- "*"
print(length(TSSdownstream500GR))
#[1] 27204
TSSdownstream500GR_maskOvelaps <- findOverlaps(mask, TSSdownstream500GR, ignore.strand = TRUE, select = "all")
TSSdownstream500GR <- TSSdownstream500GR[-subjectHits(TSSdownstream500GR_maskOvelaps)]
print(length(TSSdownstream500GR))
#[1] 24554

# Gene terminators (TTS+1 bp to TSS+500 bp)
source("/projects/ajt200/Rfunctions/TTSplus.R")
genesGRterm <- GRanges(seqnames = genes$chr,
                       ranges = IRanges(start = genes$start,
                                        end = genes$end),
                       strand = genes$strand)
terminatorsGR <- TTSplus(genesGRterm, upstream = -1, downstream = 500)
seqlevels(terminatorsGR) <- sub("", "Chr", seqlevels(terminatorsGR))
strand(terminatorsGR) <- "*"
print(length(terminatorsGR))
#[1] 27204
terminatorsGR_maskOvelaps <- findOverlaps(mask, terminatorsGR, ignore.strand = TRUE, select = "all")
terminatorsGR <- terminatorsGR[-subjectHits(terminatorsGR_maskOvelaps)]
print(length(terminatorsGR))
#[1] 24547

# Gene 3´ ends (TTS-500 bp to TTS)
TTSupstream500GR <- TTSplus(genesGRterm, upstream = 499, downstream = 0)
seqlevels(TTSupstream500GR) <- sub("", "Chr", seqlevels(TTSupstream500GR))
strand(TTSupstream500GR) <- "*"
print(length(TTSupstream500GR))
#[1] 27204
TTSupstream500GR_maskOvelaps <- findOverlaps(mask, TTSupstream500GR, ignore.strand = TRUE, select = "all")
TTSupstream500GR <- TTSupstream500GR[-subjectHits(TTSupstream500GR_maskOvelaps)]
print(length(TTSupstream500GR))
#[1] 24550

# Import exons as GRanges object
genes_exons <- read.table("/projects/ajt200/TAIR10/all_exons.txt",
                          header = T)
exons <- genes_exons[genes_exons$ge.ex == "exon",]
levels(exons$strand) <- c("-", "*", "+")
exonsGR <- GRanges(seqnames = exons$chr,
                   ranges = IRanges(start = exons$start, end = exons$end),
                   strand = "*")
print(length(exonsGR))
#[1] 152155
exonsGR_maskOvelaps <- findOverlaps(mask, exonsGR, ignore.strand = TRUE, select = "all")
exonsGR <- exonsGR[-subjectHits(exonsGR_maskOvelaps)]
print(length(exonsGR))
#[1] 137898

# Import introns tables and convert to GRanges object
intronsPlus <- read.table("/projects/ajt200/TAIR10/all_plus_introns.txt", header = T)
intronsMinus <- read.table("/projects/ajt200/TAIR10/all_minus_introns.txt", header = T)
intronsPlusGR <- GRanges(seqnames = paste0("Chr", intronsPlus$chr),
                         ranges = (IRanges(start = intronsPlus$all.intron.starts+1,
                                           end = intronsPlus$all.intron.stops-1)),
                         strand = "*")
intronsMinusGR <- GRanges(seqnames = paste0("Chr", intronsMinus$chr),
                          ranges = (IRanges(start = intronsMinus$all.intron.stops+1,
                                            end = intronsMinus$all.intron.starts-1)),
                          strand = "*")
intronsGR <- sort(append(intronsPlusGR, intronsMinusGR), by = ~ seqnames + start + end)
print(length(intronsGR))
#[1] 119808
intronsGR_maskOvelaps <- findOverlaps(mask, intronsGR, ignore.strand = TRUE, select = "all")
intronsGR <- intronsGR[-subjectHits(intronsGR_maskOvelaps)]
print(length(intronsGR))
#[1] 108606

TEs <- read.table(file = "/projects/ajt200/TAIR10/TAIR10_Buisine_TEs_strand_tab_ann.txt", header=T)
TEsGR <- GRanges(seqnames = TEs$Chr, ranges = IRanges(start = TEs$start, end = TEs$end), strand = "*")
print(length(TEsGR))
#[1] 31189
TEsGR_maskOvelaps <- findOverlaps(mask, TEsGR, ignore.strand = TRUE, select = "all")
TEsGR <- TEsGR[-subjectHits(TEsGR_maskOvelaps)]
print(length(TEsGR))
#[1] 28562

otherNames <- c("MSH4_Rep1GR", "REC8_HA_Rep1GR", "nucleRnucsGR",
                "SPO11GR", "SPO11_ChIP4GR",
                "H3K4me3GR", "H3K9me2GR",
                "genesGR", "promotersGR", "terminatorsGR",
                "TSSdownstream500GR", "TTSupstream500GR",
                "exonsGR", "intronsGR", "TEsGR")

grl <- GRangesList("MSH4_Rep1GR" = MSH4_Rep1GR, "REC8_HA_Rep1GR" = REC8_HA_Rep1GR, "nucleRnucsGR" = nucleRnucsGR,
                   "SPO11GR" = SPO11GR, "SPO11_ChIP4GR" = SPO11_ChIP4GR,
                   "H3K4me3GR" = H3K4me3GR, "H3K9me2GR" = H3K9me2GR,
                   "genesGR" = genesGR, "promotersGR" = promotersGR, "terminatorsGR" = terminatorsGR,
                   "TSSdownstream500GR" = TSSdownstream500GR, "TTSupstream500GR" = TTSupstream500GR,
                   "exonsGR" = exonsGR, "intronsGR" = intronsGR, "TEsGR" = TEsGR)

# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions in B as in A
set.seed(123)
ptCOsOtherPerChrom <- lapply(seq_along(grl), function(x) {
  permTest(A = popCOsGR, B = grl[[x]], genome = genome, mask = mask,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE, per.chromosome = TRUE,
           evaluate.function = numOverlaps, count.once = TRUE,
           ntimes = 10000, mc.set.seed = FALSE, mc.cores = detectCores())
})

for(i in 1:length(ptCOsOtherPerChrom)) {
  assign(paste0(otherNames[i]), ptCOsOtherPerChrom[[i]])
}
save(ptCOsOtherPerChrom,
     file = paste0(outDir,
                   "permTest_", popName, "_COs_vs_others.RData"))

# Summarise results in a table
noOfFeatures <- NULL
expected <- NULL
observed <- NULL
pval <- NULL
zscore <- NULL
for(i in 1:length(ptCOsOtherPerChrom)) {
  noOfFeaturesi <- print(length(grl[[i]]))
  noOfFeatures <- c(noOfFeatures, noOfFeaturesi)
  expectedi <- print(round(mean(ptCOsOtherPerChrom[[i]]$numOverlaps$permuted)))
  expected <- c(expected, expectedi)
  observedi <- print(ptCOsOtherPerChrom[[i]]$numOverlaps$observed)
  observed <- c(observed, observedi)
  pvali <- print(round(ptCOsOtherPerChrom[[i]]$numOverlaps$pval, 4))
  pval <- c(pval, pvali)
  zscorei <- print(round(ptCOsOtherPerChrom[[i]]$numOverlaps$zscore, 4))
  zscore <- c(zscore, zscorei)
}
ptCOsOtherPerChromDataFrame <- cbind(noOfFeatures, expected, observed, pval, zscore)
write.table(ptCOsOtherPerChromDataFrame,
            file = paste0(outDir,
                          "permTest_", popName, "_COs_vs_others_DataFrame.txt"),
            sep = "\t", row.names = F)

# plot graphical summaries of results
for(i in 1:length(ptCOsOtherPerChrom)) {
  pdf(file = paste0(plotDir, otherNames[i], "_permTest_nperm10000_", popName, "_COs_perChrom.pdf"), width = 10, height = 7)
  plot(ptCOsOtherPerChrom[[i]], main = paste0(popName, " COs vs ", otherNames[i]), xlab = "Number of overlaps", ylab = "Density")
  dev.off()

  # Using the localZScore() function, evaluate whether the association between COs and other is highly dependent on their exact position
  lz_1kb <- localZScore(pt = ptCOsOtherPerChrom[[i]], A = popCOsGR, B = grl[[i]],
                        window = 1000, step = 50, count.once = TRUE)
  lz_10kb <- localZScore(pt = ptCOsOtherPerChrom[[i]], A = popCOsGR, B = grl[[i]],
                         window = 10000, step = 500, count.once = TRUE)
  lz_custom <- localZScore(pt = ptCOsOtherPerChrom[[i]], A = popCOsGR, B = grl[[i]],
                           window = 10*mean(width(popCOsGR)), step = mean(width(popCOsGR))/2, count.once = TRUE)
  win <- as.character(round((10*mean(width(popCOsGR)))/1000))
  step <- as.character(round(mean(width(popCOsGR))/2))
  pdf(file = paste0(plotDir, otherNames[i], "_localZscore_permTest_nperm10000_", popName, "_COs_w1kb_s50bp_w10kb_s500bp_w", win ,"kb_s", step, "bp_perChrom.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0(popName, " COs vs ", otherNames[i], " (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(popName, " COs vs ", otherNames[i], " (1-kb shift)"))
  plot(lz_10kb, main = paste0(popName, " COs vs ", otherNames[i], " (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(popName, " COs vs ", otherNames[i], " (10-kb shift)"))
  plot(lz_custom, main = paste0(popName, " COs vs ", otherNames[i], " (~", win, "-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(popName, " COs vs ", otherNames[i], " (~", win, "-kb shift)"))
  dev.off()
}

