#!/applications/R/R-3.5.0/bin/Rscript

#
# author: Andy Tock (adapted from original script logistic2.R by Tom Hardcastle)
# contact: ajt200@cam.ac.uk
# date: 04.05.2020
# 

# Build GLM with GBS-derived Col/Ler crossovers as the response variable and
# SPO11-1-oligos, MNase-seq-derived nucleosome occupancy, H3K4me3 ChIP-seq,
# BS-seq-derived DNA methylation, Col/Ler SNPs, and interval width as pedictor variables

library(segmentSeq)
library(glm2)
library(MASS)

# Genomic definitions
genomeFAI <- read.table("/projects/ajt200/TAIR10/TAIR10_chr_all.fa.fai",
                        colClasses = c(rep(NA, 2), rep("NULL", 3)), header = F)[1:5,]
chrs <- paste0("Chr", genomeFAI$V1)
chrLens <- genomeFAI$V2

# Load SNPs (481252 SNPs in BC.complete.tiger.txt)
# and create a GRanges object of 481257 intervals within which a crossover can be detected
# There are 5 more intervals than SNPs due to the interval
# between the last SNP in each chromosome and the end of that chromosome
SNPs <- read.table("BC.complete.tiger.txt",
                   header = T)
print(dim(SNPs))
#[1] 481252      6
SNPsGR <- do.call("c", lapply(1:5, function(chr) {
  chrSNPs <- SNPs[SNPs[,1] == chr,]
  GRanges(seqnames = chrs[chr],
          ranges = IRanges(start = c(1, chrSNPs[,2]),
                           end = c(chrSNPs[,2], chrLens[chr])))
}))
print(SNPsGR)

# Load 3320 wild type crossovers
COsBED <- read.table("/projects/ajt200/GBS_CO/HS_CU_080617/wt/GBS_crossovers.bed",
                     header = F)
# Note addition of 1 to start coordinates as these are in BED format
# (i.e., 0-based start coordinates)
# NOTE: Tom merged overlapping COs (i.e., reduce(COsGR), which gives 3115 crossovers) 
COsGR <- GRanges(seqnames = COsBED$V1,
                 ranges = IRanges(start = COsBED$V2+1,
                                  end = COsBED$V3))
COsGR <- reduce(COsGR)

# Identify which SNP intervals overlap at least one crossover interval
# Tom set overlapType to "within", meaning that called overlaps correspond to
# SNP intervals that are completely contained by at least one crossover interval
COevents <- getOverlaps(coordinates = SNPsGR,
                        segments = COsGR,
                        overlapType = "within",
                        whichOverlaps = FALSE,
                        ignoreStrand = TRUE)
# Number of SNP intervals completely contained by at least one crossover interval
print(sum(COevents))
#[1] 4484
# Number of SNP intervals not completely contained by at least one crossover interval
print(sum(!COevents))
#[1] 476773

# Define regionsGR object as SNP intervals that are
# not completely contained by at least one crossover interval
# + all crossover intervals

regionsGR <- c(SNPsGR[!COevents], COsGR)
regionsGR <- sort(regionsGR)
print(length(regionsGR))
#[1] 479888

# Load SPO11-1-oligos, nuclesome occupancy, H3K4me3 ChIP-seq and DNA methylation data
# and calculate sum and mean in each region

## with z-score standardisation
#SPO11 <- read.table(paste0("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/",
#                           "log2wtSPO11oligoMeanAllRepsNakedDNA_norm_allchrs_coverage_coord_tab.bed"),
#                    colClasses = c(NA, rep("NULL", 2), NA),
#                    header = F)
# without z-score standardisation
SPO11 <- read.table(paste0("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/",
                           "log2wtSPO11oligoMeanAllRepsNakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed"),
                    colClasses = c(NA, rep("NULL", 2), NA),
                    header = F)
colnames(SPO11) <- c("chr", "signal")
#SPO11_Rep1 <- read.table(paste0("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/",
#                                "log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed"),
#                         colClasses = c(NA, rep("NULL", 2), NA),
#                         header = F)
#SPO11_Rep2 <- read.table(paste0("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/",
#                                "log2wtSPO11oligoRPI3NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed"),
#                         colClasses = c(NA, rep("NULL", 2), NA),
#                         header = F)
#SPO11_Rep3 <- read.table(paste0("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/",
#                                "log2wtSPO11oligoRPI8NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed"),
#                         colClasses = c(NA, rep("NULL", 2), NA),
#                         header = F)
#SPO11_df <- data.frame(Rep1 = SPO11_Rep1$V4,
#                       Rep2 = SPO11_Rep2$V4,
#                       Rep3 = SPO11_Rep3$V4,
#                       stringsAsFactors = F)
#SPO11 <- data.frame(chr = SPO11_Rep1$V1,
#                    signal = rowMeans(SPO11_df),
#                    stringsAsFactors = F)
splitSPO11 <- list()
for(x in 1:5) {
  chr_SPO11 <- SPO11[SPO11$chr == chrs[x],]$signal
  chr_regionsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  splitSPO11 <- c(splitSPO11,
                  split(x = chr_SPO11,
                        f = c(rep(1:length(chr_regionsGR),
                                  times = width(chr_regionsGR)-1),
                              length(chr_regionsGR))))
}
meanSPO11 <- sapply(splitSPO11, mean)
sumSPO11 <- sapply(splitSPO11, sum)

## with z-score standardisation
#MNase <- read.table(paste0("/projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_trim51_input/log2ChIPinput/",
#                           "log2wtNucNakedDNA_norm_allchrs_coverage_coord_tab.bed"),
#                    colClasses = c(NA, rep("NULL", 2), NA),
#                    header = F)
# without z-score standardisation
MNase <- read.table(paste0("/projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/noZscore/",
                           "log2wtNucNakedDNAuntrimmed_norm_allchrs_coverage_coord_tab_noZscore.bed"),
                    colClasses = c(NA, rep("NULL", 2), NA),
                    header = F)
colnames(MNase) <- c("chr", "signal")
splitMNase <- list()
for(x in 1:5) {
  chr_MNase <- MNase[MNase$chr == chrs[x],]$signal
  chr_regionsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  splitMNase <- c(splitMNase,
                  split(x = chr_MNase,
                        f = c(rep(1:length(chr_regionsGR),
                                  times = width(chr_regionsGR)-1),
                              length(chr_regionsGR))))
}
meanMNase <- sapply(splitMNase, mean)
sumMNase <- sapply(splitMNase, sum)

## with z-score standardisation
#H3K4me3 <- read.table(paste0("/projects/ajt200/BAM_masters/H3K4me3/WT/coverage/log2ChIPinput/",
#                             "log2wtH3K4me3ChIPwtH3K9me2input_norm_allchrs_coverage_coord_tab.bed"),
#                      colClasses = c(NA, rep("NULL", 2), NA),
#                      header = F)
# without z-score standardisation
H3K4me3 <- read.table(paste0("/projects/ajt200/BAM_masters/H3K4me3/WT/coverage/log2ChIPinput/noZscore/",
                             "log2wtH3K4me3ChIPwtH3K9me2input_noZscore_norm_allchrs_coverage_coord_tab.bed"),
                      colClasses = c(NA, rep("NULL", 2), NA),
                      header = F)
colnames(H3K4me3) <- c("chr", "signal")
#H3K4me3_Rep1 <- read.table(paste0("/projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/noZscore/",
#                                  "WT_H3K4me3_ChIP12_log2ChIPinput_norm_allchrs_coverage_coord_tab_noZscore.bed"),
#                           colClasses = c(NA, rep("NULL", 2), NA),
#                           header = F)
#H3K4me3_Rep2 <- read.table(paste0("/projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/noZscore/",
#                                  "WT_H3K4me3_ChIP14_log2ChIPinput_norm_allchrs_coverage_coord_tab_noZscore.bed"),
#                           colClasses = c(NA, rep("NULL", 2), NA),
#                           header = F)
#H3K4me3_Rep3 <- read.table(paste0("/projects/ajt200/BAM_masters/H3K4me3/replicates/coverage/log2ChIPinput/noZscore/",
#                                  "WT_H3K4me3_ChIP15_log2ChIPinput_norm_allchrs_coverage_coord_tab_noZscore.bed"),
#                           colClasses = c(NA, rep("NULL", 2), NA),
#                           header = F)
#H3K4me3_df <- data.frame(
#                         Rep1 = H3K4me3_Rep1$V4,
#                         Rep2 = H3K4me3_Rep2$V4,
#                         Rep3 = H3K4me3_Rep3$V4,
#                         stringsAsFactors = F)
#H3K4me3 <- data.frame(chr = H3K4me3_Rep1$V1,
#                      signal = rowMeans(H3K4me3_df),
#                      stringsAsFactors = F)
splitH3K4me3 <- list()
for(x in 1:5) {
  chr_H3K4me3 <- H3K4me3[H3K4me3$chr == chrs[x],]$signal
  chr_regionsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  splitH3K4me3 <- c(splitH3K4me3,
                  split(x = chr_H3K4me3,
                        f = c(rep(1:length(chr_regionsGR),
                                  times = width(chr_regionsGR)-1),
                              length(chr_regionsGR))))
}
meanH3K4me3 <- sapply(splitH3K4me3, mean)
sumH3K4me3 <- sapply(splitH3K4me3, sum)

# DNA methylation
CG <- read.table(paste0("/home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/",
                        "GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed"),
                 colClasses = c(rep(NA, 2), "NULL", NA),
                 header = F)
CHG <- read.table(paste0("/home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/",
                         "GSM980986_WT_rep2_CHG.wig.bed.gr.tab.bed"),
                  colClasses = c(rep(NA, 2), "NULL", NA),
                  header = F)
CHH <- read.table(paste0("/home/ajt200/BS_Seq/Stroud_2013/WT_rep2/wig/bed/",
                         "GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed"),
                  colClasses = c(rep(NA, 2), "NULL", NA),
                  header = F)
## mC method below would lead to big underestimation of mean methylation
## levels where CG-context and/or CHG-context DNA metyhlation information is
## absent in a given region
#mC <- rbind(CG, CHG, CHH)
#mC <- mC[with(mC, order(V1, V2)),]
#meanmC <- NULL
meanCG <- NULL
meanCHG <- NULL
meanCHH <- NULL
meanDNAmeth <- NULL
for(x in 1:5) {
  windowsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  windowsGR <- GRanges(seqnames = chrs[x],
                       ranges = IRanges(start = start(windowsGR),
                                        end = end(windowsGR)-1))
  end(windowsGR[length(windowsGR)]) <- chrLens[x]
  # Define DNA methylation coordinates as GRanges objects
  # and calculate mean methylation proportions in each window
  ## mC
  #chr_mC <- mC[mC[,1] == paste0("chr", x),]
  #chr_mC_GR <- GRanges(seqnames = chrs[x],
  #                     ranges = IRanges(start = chr_mC[,2],
  #                                      width = 1),
  #                     strand = "*")
  #mCoverlaps <- getOverlaps(coordinates = windowsGR,
  #                          segments = chr_mC_GR,
  #                          overlapType = "overlapping",
  #                          whichOverlaps = TRUE,
  #                          ignoreStrand = TRUE)
  #mCwinVals <- sapply(mCoverlaps, function(x) mean(as.numeric(chr_mC[,3][x])))
  # CG
  chr_CG <- CG[CG[,1] == paste0("chr", x),]
  chr_CG_GR <- GRanges(seqnames = chrs[x],
                       ranges = IRanges(start = chr_CG[,2],
                                        width = 1),
                       strand = "*")
  CGoverlaps <- getOverlaps(coordinates = windowsGR,
                            segments = chr_CG_GR,
                            overlapType = "overlapping",
                            whichOverlaps = TRUE,
                            ignoreStrand = TRUE)
  CGwinVals <- sapply(CGoverlaps, function(x) mean(as.numeric(chr_CG[,3][x])))
  # CHG
  chr_CHG <- CHG[CHG[,1] == paste0("chr", x),]
  chr_CHG_GR <- GRanges(seqnames = chrs[x],
                        ranges = IRanges(start = chr_CHG[,2],
                                         width = 1),
                        strand = "*")
  CHGoverlaps <- getOverlaps(coordinates = windowsGR,
                             segments = chr_CHG_GR,
                             overlapType = "overlapping",
                             whichOverlaps = TRUE,
                             ignoreStrand = TRUE)
  CHGwinVals <- sapply(CHGoverlaps, function(x) mean(as.numeric(chr_CHG[,3][x])))
  # CHH
  chr_CHH <- CHH[CHH[,1] == paste0("chr", x),]
  chr_CHH_GR <- GRanges(seqnames = chrs[x],
                        ranges = IRanges(start = chr_CHH[,2],
                                         width = 1),
                        strand = "*")
  CHHoverlaps <- getOverlaps(coordinates = windowsGR,
                             segments = chr_CHH_GR,
                             overlapType = "overlapping",
                             whichOverlaps = TRUE,
                             ignoreStrand = TRUE)
  CHHwinVals <- sapply(CHHoverlaps, function(x) mean(as.numeric(chr_CHH[,3][x])))
  # Mean of all 3 contexts
  CwinVals <- sapply(seq_along(CGwinVals), function(x) {
                mean(c(CGwinVals[x], CHGwinVals[x], CHHwinVals[x]))
              })
  # Combine values from all chromosomes
  #meanmC <- c(meanmC, mCwinVals)
  meanCG <- c(meanCG, CGwinVals)
  meanCHG <- c(meanCHG, CHGwinVals)
  meanCHH <- c(meanCHH, CHHwinVals)
  meanDNAmeth <- c(meanDNAmeth, CwinVals)
}

# Load larger Col/Ler SNPs data set (519834 SNPs) for calculating SNP frequency
# within each interval in regionsGR
newSNPs <- read.table("/projects/ajt200/GBS_CO/msh2/SNPs/collerF2.complete.tiger.txt",
                      header = T)[,1:2]
colnames(newSNPs) <- c("chr", "pos")
print(dim(newSNPs))
#[1] 519834      6
SNPs <- SNPs[,1:2]
colnames(SNPs) <- c("chr", "pos")
print(dim(SNPs))
#[1] 481252      2
allSNPs <- rbind(SNPs[,1:2], newSNPs[,1:2]) 
print(dim(allSNPs))
#[1] 1001086       2
allSNPs <- unique(allSNPs)
print(dim(allSNPs))
#[1] 534775      2
allSNPsGR <- GRanges(seqnames = paste0("Chr", allSNPs[,1]),
                     ranges = IRanges(start = allSNPs[,2],
                                      width = 1))
# Count SNPs within each interval in which a crossover could potentially have been detected
# and divide by the width of the corresponding interval
meanSNPs <- NULL
sumSNPs <- NULL
for(x in 1:5) {
  windowsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  chr_allSNPsGR <- allSNPsGR[seqnames(allSNPsGR) == chrs[x]]
  winSNPs <- countOverlaps(query = windowsGR,
                           subject = chr_allSNPsGR,
                           type = "any",
                           ignore.strand = TRUE)
  meanSNPs <- c(meanSNPs, winSNPs/width(windowsGR))
  sumSNPs <- c(sumSNPs, winSNPs)
} 

# Count SNPs in windows if width winSize nt, with a step of stepSize nt
winSize <- 1000
stepSize <- 1
SNPsPB <- data.frame()
for(x in 1:5) {
  # Define sliding windows of width winSize nt,
  # with a step of stepSize nt
  ## Note: the commented-out code would create windows of winSize nt only,
  ## whereas the active code creates windows decreasing from winSize nt to stepSize nt
  ## at the right-hand end of each chromosome ( from chrLens[x]-winSize to chrLens[x] ),
  winStarts <- seq(from = 1,
                   to = chrLens[x],
#                   to = chrLens[x]-winSize,
                   by = stepSize)
  stopifnot(winStarts[length(winStarts)] == chrLens[x])
#  if(chrLens[x] - winStarts[length(winStarts)] >= winSize) {
#    winStarts <- c(winStarts,
#                   winStarts[length(winStarts)]+stepSize)
#  }
  winEnds <- seq(from = winStarts[1]+winSize-1,
                 to = chrLens[x],
                 by = stepSize)
  stopifnot(winEnds[length(winEnds)] == chrLens[x])
  winEnds <- c(winEnds,
               rep(chrLens[x], times = length(winStarts)-length(winEnds)))
  stopifnot(length(winStarts) == length(winEnds))

  windowsGR <- GRanges(seqnames = chrs[x],
                       ranges = IRanges(start = winStarts,
                                        end = winEnds),
                       strand = "*")

  # Count SNPs in sliding windows
  chr_allSNPsGR <- allSNPsGR[seqnames(allSNPsGR) == chrs[x]]
  winSNPs <- countOverlaps(query = windowsGR,
                           subject = chr_allSNPsGR,
                           type = "any",
                           ignore.strand = TRUE)
  SNPsPBDF <- data.frame(chr = chrs[x],
                         signal = winSNPs/width(windowsGR)) 
  SNPsPB <- rbind(SNPsPB, SNPsPBDF) 
}

splitSNPsPB <- list()
for(x in 1:5) {
  chr_SNPsPB <- SNPsPB[SNPsPB$chr == chrs[x],]$signal
  chr_regionsGR <- regionsGR[seqnames(regionsGR) == chrs[x]]
  splitSNPsPB <- c(splitSNPsPB,
                   split(x = chr_SNPsPB,
                         f = c(rep(1:length(chr_regionsGR),
                                   times = width(chr_regionsGR)-1),
                               length(chr_regionsGR))))
}
meanSNPsPB <- sapply(splitSNPsPB, mean)
sumSNPsPB <- sapply(splitSNPsPB, sum)

# Create data object for model
dat <- cbind.data.frame(CO = regionsGR %in% COsGR,
                        meanSPO11 = meanSPO11,
                        meanMNase = meanMNase,
                        meanH3K4me3 = meanH3K4me3,
                        meanDNAmeth = meanDNAmeth,
                        meanSNPs = meanSNPs,
                        meanSNPsPB = meanSNPsPB,
                        sumSPO11 = sumSPO11,
                        sumMNase = sumMNase,
                        sumH3K4me3 = sumH3K4me3,
                        sumSNPs = sumSNPs,
                        sumSNPsPB = sumSNPsPB,
                        width = width(regionsGR))
save(dat, file = "df_for_GLM_binomial_logit_1kbWin_1bpStep_SNPsPB.RData")

# Discard rows with missing data
dat <- dat[!is.na(dat$meanDNAmeth),]

# Build binomial GLM with "logit" link function
glmCO <- glm2(formula = CO ~ (meanSPO11 + meanMNase + meanH3K4me3 + meanDNAmeth + meanSNPsPB + width)^2,
              family = binomial(link="logit"),
              control = glm.control(maxit = 100000),
              data = dat)
#Warning message:
#glm.fit2: fitted probabilities numerically 0 or 1 occurred
## Variables as in GLM built by Ian for Choi et al. (2018) Genome Res.:
#glmCO <- glm2(formula = CO ~ meanSPO11 + meanMNase + meanH3K4me3 + meanDNAmeth + width +
#              meanSPO11:meanMNase + meanSPO11:meanH3K4me3 + meanSPO11:width +
#              meanMNase:meanH3K4me3 + meanMNase:meanDNAmeth +
#              meanH3K4me3:meanDNAmeth +
#              meanDNAmeth:width,
#              family = binomial(link="logit"),
#              data = dat)
#Warning message:
#glm.fit2: fitted probabilities numerically 0 or 1 occurred
#glm_select <- glmCO

glm_stepAIC <- stepAIC(object = glmCO, direction = "both")
print("stepAIC-selected model formula:")
print(glm_stepAIC$formula)
glm_select <- glm2(formula = glm_stepAIC$formula, 
                   family = binomial(link="logit"),
                   control = glm.control(maxit = 100000),
                   data = dat)
glm_summary <- summary(glm_select)
glm_coeffs <- glm_summary$coefficients
glm_predict <- predict(glm_select, type = "response")
glm_formula <- glm_select$formula
save(glm_stepAIC, file = "GLM_binomial_logit_1kbWin_1bpStep_SNPsPB_stepAIC.RData")
save(glm_select, file = "GLM_binomial_logit_1kbWin_1bpStep_SNPsPB.RData")
save(glm_summary, file = "GLM_binomial_logit_1kbWin_1bpStep_SNPsPB_summary.RData")
write.csv(glm_coeffs, file = "GLM_binomial_logit_1kbWin_1bpStep_SNPsPB_coeff.csv")
save(glm_predict, file = "GLM_binomial_logit_1kbWin_1bpStep_SNPsPB_predict.RData")
save(glm_formula, file = "GLM_binomial_logit_1kbWin_1bpStep_SNPsPB_stepAIC_selected_formula.RData")

# Plot observed and predicted crossovers for regionsGR grouped into hexiles

# meanSPO11 hexiles
sort(dat$meanSPO11)[round(1:6*(nrow(dat)/6))]
#[1] -1.205848832 -0.723219146 -0.372348353  0.001650497  0.467101543 [6]  4.399179154
levels(cut(x = dat$meanSPO11,
           breaks = c(min(dat$meanSPO11, na.rm = T),
                      sort(dat$meanSPO11)[round(1:6*(nrow(dat)/6))])))
#[1] "(-4.47,-1.21]"    "(-1.21,-0.723]"   "(-0.723,-0.372]"  "(-0.372,0.00165]"
#[5] "(0.00165,0.467]"  "(0.467,4.4]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$meanSPO11,
                      breaks = c(min(dat$meanSPO11, na.rm = T),
                                 sort(dat$meanSPO11)[round(1:6*(nrow(dat)/6))])))
pdf("GLM_binomial_logit_1kbWin_1bpStep_SNPsPB_SPO11_hexiles_CO_boxplot.pdf")
boxplot(
  lapply(ssID, function(x) {
    sapply(1:100, function(ii) {
      sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
      sum( dat$width[x] ) * 1e6
    })
  }),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "SPO11-1-oligo hexiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topleft",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()

# meanSNPsPB hexiles
sort(dat$meanSNPsPB)[round(1:6*(nrow(dat)/6))]
#[1] 0.002004717 0.004002079 0.006371795 0.010021739 0.015288462 0.062791667
#[1] 0.002355150 0.003801124 0.005421569 0.007500909 0.010410427 0.038117165
levels(cut(x = dat$meanSNPsPB,
           breaks = c(min(dat$meanSNPsPB, na.rm = T),
                      sort(dat$meanSNPsPB)[round(1:6*(nrow(dat)/6))])))
#[1] "(1.13e-05,0.002]" "(0.002,0.004]"    "(0.004,0.00637]"  "(0.00637,0.01]"
#[5] "(0.01,0.0153]"    "(0.0153,0.0628]"
#[1] "(1.17e-06,0.00236]" "(0.00236,0.0038]"   "(0.0038,0.00542]"
#[4] "(0.00542,0.0075]"   "(0.0075,0.0104]"    "(0.0104,0.0381]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$meanSNPsPB,
                      breaks = c(min(dat$meanSNPsPB, na.rm = T),
                                 sort(dat$meanSNPsPB)[round(1:6*(nrow(dat)/6))])))
pdf("GLM_binomial_logit_1kbWin_1bpStep_SNPsPB_SNPsPerBase_hexiles_CO_boxplot.pdf")
boxplot(
  lapply(ssID, function(x) {
    sapply(1:100, function(ii) {
      sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
      sum( dat$width[x] ) * 1e6
    })
  }),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "SNPs per base hexiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topleft",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()

# meanMNase hexiles
sort(dat$meanMNase)[round(1:6*(nrow(dat)/6))]
#[1] -1.17183549 -0.51678188 -0.02808389  0.38526674  0.81293932  3.12292755
levels(cut(x = dat$meanMNase,
           breaks = c(min(dat$meanMNase, na.rm = T),
                      sort(dat$meanMNase)[round(1:6*(nrow(dat)/6))])))
#[1] "(-4.76,-1.17]"    "(-1.17,-0.517]"   "(-0.517,-0.0281]" "(-0.0281,0.385]"
#[5] "(0.385,0.813]"    "(0.813,3.12]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$meanMNase,
                      breaks = c(min(dat$meanMNase, na.rm = T),
                                 sort(dat$meanMNase)[round(1:6*(nrow(dat)/6))])))
pdf("GLM_binomial_logit_1kbWin_1bpStep_SNPsPB_MNase_hexiles_CO_boxplot.pdf")
boxplot(
  lapply(ssID, function(x) {
    sapply(1:100, function(ii) {
      sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
      sum( dat$width[x] ) * 1e6
    })
  }),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "Nucleosomes hexiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topright",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()

# meanH3K4me3 hexiles
sort(dat$meanH3K4me3)[round(1:6*(nrow(dat)/6))]
#[1] -1.8706652 -1.3072388 -0.7229520 -0.1216027  0.5737178  4.1235247
levels(cut(x = dat$meanH3K4me3,
           breaks = c(min(dat$meanH3K4me3, na.rm = T),
                      sort(dat$meanH3K4me3)[round(1:6*(nrow(dat)/6))])))
#[1] "(-5.03,-1.87]"   "(-1.87,-1.31]"   "(-1.31,-0.723]"  "(-0.723,-0.122]"
#[5] "(-0.122,0.574]"  "(0.574,4.12]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$meanH3K4me3,
                      breaks = c(min(dat$meanH3K4me3, na.rm = T),
                                 sort(dat$meanH3K4me3)[round(1:6*(nrow(dat)/6))])))
pdf("GLM_binomial_logit_1kbWin_1bpStep_SNPsPB_H3K4me3_hexiles_CO_boxplot.pdf")
boxplot(
  lapply(ssID, function(x) {
    sapply(1:100, function(ii) {
      sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
      sum( dat$width[x] ) * 1e6
    })
  }),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "H3K4me3 hexiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topright",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()

# meanDNAmeth quintiles
sort(dat$meanDNAmeth)[round(1:6*(nrow(dat)/6))]
#[1] 0.000000000 0.002414189 0.014161412 0.186714907 0.442316717 1.000000000
levels(cut(x = dat$meanDNAmeth,
           breaks = c(
                      sort(dat$meanDNAmeth)[round(1:6*(nrow(dat)/6))])))
#[1] "(0,0.00241]"      "(0.00241,0.0142]" "(0.0142,0.187]"   "(0.187,0.442]"
#[5] "(0.442,1]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$meanDNAmeth,
                      breaks = c(
                                 sort(dat$meanDNAmeth)[round(1:6*(nrow(dat)/6))])))
pdf("GLM_binomial_logit_1kbWin_1bpStep_SNPsPB_DNAmeth_quintiles_CO_boxplot.pdf")
boxplot(
  lapply(ssID, function(x) {
    sapply(1:100, function(ii) {
      sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
      sum( dat$width[x] ) * 1e6
    })
  }),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "DNA methylation quintiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topright",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()

# width hexiles
sort(dat$width)[round(1:6*(nrow(dat)/6))]
#[1]     57    107    181    314    634 127250
levels(cut(x = dat$width,
           breaks = c(min(dat$width, na.rm = T),
                      sort(dat$width)[round(1:6*(nrow(dat)/6))])))
#[1] "(4,57]"         "(57,107]"       "(107,181]"      "(181,314]"
#[5] "(314,634]"      "(634,1.27e+05]"
ssID <- split(x = 1:nrow(dat),
              f = cut(x = dat$width,
                      breaks = c(min(dat$width, na.rm = T),
                                 sort(dat$width)[round(1:6*(nrow(dat)/6))])))
pdf("GLM_binomial_logit_1kbWin_1bpStep_SNPsPB_width_hexiles_CO_boxplot.pdf")
boxplot(
  lapply(ssID, function(x) {
    sapply(1:100, function(ii) {
      sum( rbinom(n = length(x), size = 1, prob = glm_predict[x]) ) /
      sum( dat$width[x] ) * 1e6
    })
  }),
  main = "CO~(SPO11-1+nucleosomes+H3K4me3+DNAmeth+SNPsPB+width)^2",
  xlab = "Width hexiles",
  ylab = "Crossovers per Mb",
  cex.axis = 0.5
#  names = as.character(6:1)
)
tco <- sapply(ssID, function(x) {
  sum(dat$CO[x]) /
  sum(dat$width[x]) * 1e6
})
points(x = 1:length(ssID),
       y = tco,
       col = "red", pch = 19)
legend("topright",
       legend = c("Observed", "Predicted"),
       text.col = c("red", "black"), bty = "n")
dev.off()
