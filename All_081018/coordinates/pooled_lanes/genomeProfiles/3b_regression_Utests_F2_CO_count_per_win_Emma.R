#!/applications/R/R-3.5.0/bin/Rscript

# Use regression models to determine whether significant differences exist
# between genotypes with regard to CO frequency within genomic windows

# Resources consulted:
# http://datavoreconsulting.com/programming-tips/count-data-glms-choosing-poisson-negative-binomial-zero-inflated-poisson/
# https://www.physiology.org/doi/full/10.1152/advan.00017.2010?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%3dpubmed
# https://stats.idre.ucla.edu/r/dae/poisson-regression/
# https://stats.idre.ucla.edu/r/dae/zip/
# https://stats.idre.ucla.edu/r/dae/negative-binomial-regression/
# https://stats.idre.ucla.edu/r/dae/zinb/
# https://data.princeton.edu/wws509/r/overdispersion

# Usage from within directory containing windowed CO counts files:
# ./3b_regression_Utests_F2_CO_count_per_win_Emma.R EL_wt_ColBur_F2_2019 EL_taf4b_ColBur_F2_2019 1Mb 0.1

# Some of these packages may not already be installed
# e.g., install.packages("pscl")
library(GenomicRanges)
library(parallel)
library(MASS) # glm.nb() included
library(pscl) # zeroinfl() included

args <- commandArgs(trailingOnly = TRUE)
pop1Name <- args[1]
pop2Name <- args[2]
winName <- args[3]
FDR <- as.numeric(args[4])
FDRname <- paste0("FDR", as.character(FDR))

# Load list of tables in which each list element is a table of windowed
# crossover interval counts for one F2 individual
load(paste0("./", pop1Name,
            "_pooled_F2_CO_frequency_per_",
            winName, "_window_list.RData"))
pop1_indCOs_perWindow_list <- indCOs_perWindow_list
rm(indCOs_perWindow_list)
load(paste0("./", pop2Name,
            "_pooled_F2_CO_frequency_per_",
            winName, "_window_list.RData"))
pop2_indCOs_perWindow_list <- indCOs_perWindow_list
rm(indCOs_perWindow_list)

# Define and create new directories and subdirectories
# to contain results
outDir <- paste0("./regression_models/")
plotDir <- paste0(outDir, FDRname, "/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

UtestDir <- paste0("./U_tests/")
UtestPlotDir <- paste0(UtestDir, FDRname, "/")
system(paste0("[ -d ", UtestDir, " ] || mkdir ", UtestDir))
system(paste0("[ -d ", UtestPlotDir, " ] || mkdir ", UtestPlotDir))

# For pop1, create a list in which each list element (x) is a vector of
# F2 crossover interval counts for one window
pop1_winIndCOs_list <- mclapply(seq_along(1:dim(pop1_indCOs_perWindow_list[[1]])[1]), function(x) {
  windowedCOs <- NULL
  for(y in seq_along(pop1_indCOs_perWindow_list)) {
    windowedCOs <- c(windowedCOs, pop1_indCOs_perWindow_list[[y]]$COs[x])
  }
  windowedCOs
}, mc.cores = detectCores())
# For pop2, create a list in which each list element (x) is a vector of
# F2 crossover interval counts for one window
pop2_winIndCOs_list <- mclapply(seq_along(1:dim(pop2_indCOs_perWindow_list[[1]])[1]), function(x) {
  windowedCOs <- NULL
  for(y in seq_along(pop2_indCOs_perWindow_list)) {
    windowedCOs <- c(windowedCOs, pop2_indCOs_perWindow_list[[y]]$COs[x])
  }
  windowedCOs
}, mc.cores = detectCores())

# Calculate mean crossover interval counts for each population for plotting
pop1_winIndCOs_means <- sapply(seq_along(pop1_winIndCOs_list), function(x) {
  mean(pop1_winIndCOs_list[[x]])
})
pop2_winIndCOs_means <- sapply(seq_along(pop2_winIndCOs_list), function(x) {
  mean(pop2_winIndCOs_list[[x]])
})
# Calculate variance of crossover interval counts for each population
pop1_winIndCOs_vars <- sapply(seq_along(pop1_winIndCOs_list), function(x) {
  var(pop1_winIndCOs_list[[x]])
})
pop2_winIndCOs_vars <- sapply(seq_along(pop2_winIndCOs_list), function(x) {
  var(pop2_winIndCOs_list[[x]])
})

# Plot means vs variances
pdf(paste0(outDir, pop1Name, "_", pop2Name,
           "_crossover_means_vs_variances_in_",
           winName, "_genomic_windows.pdf"),
    height = 5, width = 10)
par(mfrow = c(1, 2))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(x = pop1_winIndCOs_means, y = pop1_winIndCOs_vars, pch = 19,
     xlab = "Mean", ylab = "Variance",
     xlim = c(min(pop1_winIndCOs_means, pop1_winIndCOs_vars),
              max(pop1_winIndCOs_means, pop1_winIndCOs_vars)),
     ylim = c(min(pop1_winIndCOs_means, pop1_winIndCOs_vars),
              max(pop1_winIndCOs_means, pop1_winIndCOs_vars)),
     main = paste0(pop1Name, " crossovers in ", winName,
                   " genomic windows"),
     cex.main = 0.8)
plot(x = pop2_winIndCOs_means, y = pop2_winIndCOs_vars, pch = 19,
     xlab = "Mean", ylab = "Variance",
     xlim = c(min(pop2_winIndCOs_means, pop2_winIndCOs_vars),
              max(pop2_winIndCOs_means, pop2_winIndCOs_vars)),
     ylim = c(min(pop2_winIndCOs_means, pop2_winIndCOs_vars),
              max(pop2_winIndCOs_means, pop2_winIndCOs_vars)),
     main = paste0(pop2Name, " crossovers in ", winName,
                   " genomic windows"),
     cex.main = 0.8)
dev.off()

# Combine above crossover interval counts lists so that each list element
# is a vector of pop1 followed by pop2 F2 crossover interval counts for one window
pops_winIndCOs_list <- lapply(seq_along(pop1_winIndCOs_list), function(x) {
  c(pop1_winIndCOs_list[[x]], pop2_winIndCOs_list[[x]])
})

# Create a factor with two levels (pop1Name and pop2Name, corresponding to genotype)
# to be used as the predictor in regression models
genotype <- factor(c(rep(pop1Name, times = length(pop1_winIndCOs_list[[1]])),
                     rep(pop2Name, times = length(pop2_winIndCOs_list[[1]]))),
                   levels = c(pop1Name, pop2Name))

# Create dataframe combining genotype and pops_winIndCOs_list
pops_winIndCOs_df <- data.frame(genotype = factor(c(rep(pop1Name, times = length(pop1_winIndCOs_list[[1]])),
                                                    rep(pop2Name, times = length(pop2_winIndCOs_list[[1]]))),
                                                  levels = c(pop1Name, pop2Name)),
                                matrix(unlist(pops_winIndCOs_list),
                                       ncol = length(pop1_winIndCOs_list),
                                       byrow = F))

# For each genomic window, fit various regression models in which genotype is the
# predictor and crossover counts are the response variable
# Determine which model has best fit

# Normal linear regression, equivalent to a t-test:
normal <- lapply(seq_along(pops_winIndCOs_list), function(x) {
  glm(formula = pops_winIndCOs_list[[x]] ~ genotype,
      family = gaussian(link = "identity"))
})
# Get estimates
normal_estimates <- lapply(seq_along(normal), function(x) {
  coef(normal[[x]])
})
# Get coefficients
normal_coef <- lapply(seq_along(normal), function(x) {
#  if( sum(pops_winIndCOs_list[[x]]) > 0 ) {
    coef(summary(normal[[x]]))
#  } else {
#    mat <- matrix(rep(NA, times = 8), ncol = 4)
#    rownames(mat) <- c("(Intercept)", paste0("genotype", pop2Name))
#    colnames(mat) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
#    mat
#  }
})
# Get P-values
normal_Pvals <- sapply(seq_along(normal), function(x) {
  coef(summary(normal[[x]]))[8]
})
# Correct for multiple testing
normal_AdjPvals <- p.adjust(p = normal_Pvals, method = "BH")
# Get Bayesian Information Criterion (BIC); the smaller, the better the fit
normal_BIC <- sapply(seq_along(normal), function(x) {
  if( !is.infinite(BIC(normal[[x]])) ) {
    BIC(normal[[x]])
  } else {
    NA
  }
})
# Get 95% confidence intervals
normal_CIs <- lapply(seq_along(normal), function(x) {
  if( sum(pops_winIndCOs_list[[x]]) > 0 ) {
   confint(normal[[x]])
  } else {
    mat <- matrix(rep(NA, times = 4), ncol = 2)
    rownames(mat) <- c("(Intercept)", paste0("genotype", pop2Name))
    colnames(mat) <- c("2.5 %", "97.5 %")
    mat
  }
})

# Evaluate model goodness-of-fit using chi-squared test based on the
# residual deviance and degrees of freedom
# A P-value > 0.05 indicates that the model fits the data 
normal_chisqPvals <- sapply(seq_along(normal), function(x) {
  1 - pchisq(summary(normal[[x]])$deviance,
             summary(normal[[x]])$df.residual)
})
print(paste0("Sum of P-values from chi-squared goodness-of-fit tests evaluating normal linear regression models for ",
             length(normal), " genomic windows:"))
print(sum(normal_chisqPvals))

# Use the model to predict mean crossover counts for each genotype
# and corresponding standard errors
# The model predicts the correct mean counts
normal_predict <- lapply(seq_along(normal), function(x) {
  cbind(data.frame(genotype = c(pop1Name, pop2Name)),
        mean = predict(normal[[x]],
                       newdata = data.frame(genotype = c(pop1Name, pop2Name)),
                       type = "response"),
        SE = predict(normal[[x]],
                     newdata = data.frame(genotype = c(pop1Name, pop2Name)),
                     type = "response",
                     se.fit = TRUE)$se.fit
       )
})


# Poisson regression:
poisson <- lapply(seq_along(pops_winIndCOs_list), function(x) {
  glm(formula = pops_winIndCOs_list[[x]] ~ genotype,
      family = poisson(link = "log"),
      maxit = 100)
})
# Get estimates
poisson_estimates <- lapply(seq_along(poisson), function(x) {
  coef(poisson[[x]])
})
# Get coefficients
poisson_coef <- lapply(seq_along(poisson), function(x) {
    coef(summary(poisson[[x]]))
})
# Get P-values
poisson_Pvals <- sapply(seq_along(poisson), function(x) {
  coef(summary(poisson[[x]]))[8]
})
# Correct for multiple testing
poisson_AdjPvals <- p.adjust(p = poisson_Pvals, method = "BH")
# Get Bayesian Information Criterion (BIC); the smaller, the better the fit
poisson_BIC <- sapply(seq_along(poisson), function(x) {
  if( !is.infinite(BIC(poisson[[x]])) ) {
    BIC(poisson[[x]])
  } else {
    NA
  }
})

# Evaluate model goodness-of-fit using chi-squared test based on the
# residual deviance and degrees of freedom
# A P-value > 0.05 indicates that the model fits the data 
poisson_chisqPvals <- sapply(seq_along(poisson), function(x) {
  1 - pchisq(summary(poisson[[x]])$deviance,
             summary(poisson[[x]])$df.residual)
})
print(paste0("Sum of P-values from chi-squared goodness-of-fit tests evaluating Poisson regression models for ",
             length(poisson), " genomic windows:"))
print(sum(poisson_chisqPvals))

# Use the model to predict mean crossover counts for each genotype
# and corresponding standard errors
# The model predicts the correct mean counts
poisson_predict <- lapply(seq_along(poisson), function(x) {
  cbind(data.frame(genotype = c(pop1Name, pop2Name)),
        mean = predict(poisson[[x]],
                       newdata = data.frame(genotype = c(pop1Name, pop2Name)),
                       type = "response"),
        SE = predict(poisson[[x]],
                     newdata = data.frame(genotype = c(pop1Name, pop2Name)),
                     type = "response",
                     se.fit = TRUE)$se.fit
       )
})

# Zero-inflated Poisson (ZIP) regression:
ZIP <- lapply(seq_along(pops_winIndCOs_list), function(x) {
  zeroinfl(formula = pops_winIndCOs_list[[x]] ~ genotype | 1,
           dist = "poisson",
           maxit = 2000)
})
print("Representative example: ZIP model for first genomic window")
print(summary(ZIP[[1]]))

# Get Bayesian Information Criterion (BIC); the smaller, the better the fit
ZIP_BIC <- sapply(seq_along(ZIP), function(x) {
  if( !is.infinite(AIC(ZIP[[x]], k = log(length(pops_winIndCOs_list[[x]])))) ) {
    AIC(ZIP[[x]], k = log(length(pops_winIndCOs_list[[x]])))
  } else {
    NA
  }
})

# Determine if normal BIC is lower than Poisson BIC (indicating better fit)
normal_BIC_vs_poisson_BIC <- sapply(seq_along(normal), function(x) {
  normal_BIC[x] < poisson_BIC[x]
})
# Determine if normal BIC is NA
normal_BIC_isNA <- sapply(seq_along(normal), function(x) {
  is.na(normal_BIC[x])
})
print(paste0("Normal BIC is < Poisson BIC for ",
             sum(normal_BIC_vs_poisson_BIC, na.rm = T), "/", length(normal), " genomic windows"))
print(paste0("Normal BIC is NA for ",
             sum(normal_BIC_isNA), "/", length(normal), " genomic windows"))

# Determine if Poisson BIC is lower than ZIP BIC (indicating better fit)
poisson_BIC_vs_ZIP_BIC <- sapply(seq_along(poisson), function(x) {
  poisson_BIC[x] < ZIP_BIC[x]
})
print(paste0("Poisson BIC is < ZIP BIC for ",
             sum(poisson_BIC_vs_ZIP_BIC, na.rm = T), "/", length(poisson), " genomic windows"))

# Note: negative binomial models could not be fit to these data
# due to absence of overdispersion 
## Negative binomial regression:
#negbin <- lapply(seq_along(pops_winIndCOs_list), function(x) {
#  glm.nb(formula = pops_winIndCOs_list[[x]] ~ genotype,
#         control = glm.control(maxit = 2000))
#})
## Zero-inflated negative binomial regression:
#ZINB <- lapply(seq_along(pops_winIndCOs_list), function(x) {
#  zeroinfl(formula = pops_winIndCOs_list[[x]] ~ genotype | 1,
#           dist = "negbin")
#})

# Plot poisson_BIC vs normal_BIC
pdf(paste0(outDir, pop1Name, "_", pop2Name,
           "_Poisson_vs_normal_BICs_for_crossovers_in_",
           winName, "_genomic_windows.pdf"),
    height = 5, width = 5)
par(mfrow = c(1, 1))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(x = 1:length(poisson_BIC), y = poisson_BIC, col = "red", type = "p", pch = 19,
     xlab = paste0("Genomic windows (", winName, ")"),
     ylab = "BIC",
     ylim = c(min(c(poisson_BIC, normal_BIC), na.rm = T),
              max(c(poisson_BIC, normal_BIC), na.rm = T)))
lines(x = 1:length(poisson_BIC), y = normal_BIC, col = "blue", type = "p", pch = 19)
legend("top",
       legend = c("Poisson regression", "Normal linear regression"),
       col = c("red", "blue"),
       text.col = c("red", "blue"),
       text.font = c(1, 1),
       ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
dev.off()


# Perfom Mann–Whitney–Wilcoxon tests comparing pop1Name and pop2Name COs
# in each genomic window
Utests <- lapply(seq_along(pop1_winIndCOs_list), function(x) {
  wilcox.test(x = pop1_winIndCOs_list[[x]],
              y = pop2_winIndCOs_list[[x]],
              alternative = "two.sided")
})
UtestPvals <- sapply(seq_along(Utests), function(x) {
  Utests[[x]]$p.val
})

# Adjust P-values by applying the Benjamini-Hochberg multiple testing correction
UtestAdjPvals <- p.adjust(p = UtestPvals, method = "BH")


# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

# Make cumulative genomes
sumchr <- cumsum(c(0, chrLens))
print(sumchr)
sumchr_tot <- sumchr[length(sumchr)]
print(sumchr_tot)

centromeres <- sapply(seq_along(centromeres), function(x) {
  centromeres[x] + sumchr[x]
})
print(centromeres)
pericenStart <- sapply(seq_along(pericenStart), function(x) {
  pericenStart[x] + sumchr[x]
})
print(pericenStart)
pericenEnd <- sapply(seq_along(pericenEnd), function(x) {
  pericenEnd[x] + sumchr[x]
})
print(pericenEnd)

# Function to plot -log10-transformed P-values and adjusted P-values derived from Fisher's exact tests
minusLog10PvalPlot <- function(xplot,
                               Pvals, PvalsCol,
                               AdjPvals, AdjPvalsCol) {
  plot(x = xplot, y = -log10(Pvals), col = PvalsCol, type = "l", lwd = 1.5,
       xlim = c(0, 1.2e8),
       ylim = c(0,
                pmax(-log10(0.05), max(-log10(Pvals), na.rm = T))),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = "")
  mtext(side = 2, line = 2.25, cex = 1, col = PvalsCol,
        text = bquote("-Log"[10]*"("*italic("P")*"-value)"))
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  abline(h = -log10(0.05), lty = 5, lwd = 1, col = PvalsCol)

  par(new = T)
  plot(x = xplot, y = -log10(AdjPvals), col = AdjPvalsCol, type = "l", lwd = 1.5,
       xlim = c(0, 1.2e8),
       ylim = c(0,
                pmax(-log10(FDR), max(-log10(AdjPvals), na.rm = T))),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(pop1Name)*" versus "*.(pop2Name)*" crossover counts per "*.(winName)),
       cex.main = 1)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1, adj = c(0.5, -3.0), xpd = NA, srt = -90, col = AdjPvalsCol,
       labels = bquote("-Log"[10]*"(BH-adjusted "*italic("P")*"-value)"))
  axis(side = 4, cex.axis = 1, lwd.tick = 1.5)
  axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"))
  mtext(side = 1, line = 2.25, cex = 1, text = "Coordinates (Mb)")
  abline(h = -log10(FDR), lty = 5, lwd = 1, col = AdjPvalsCol)
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = centromeres, lty = 3, lwd = 2, col = "black")
  box(lwd = 1.5)
}

# Function to plot Tel to Cen profile of windowed pop1 vs pop2 crossovers (one Y-axis) 
pop1Vpop2GenomePlot <- function(xplot,
                                pop1, pop1Col,
                                pop2, pop2Col,
                                Ylabel,
                                legendLoc,
                                legendLabs) {
  plot(xplot, pop1, type = "l", lwd = 1.5, col = pop1Col,
       xlim = c(0, 1.2e8),
       ylim = c(0,
                max(c(pop1, pop2))),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(pop1, pop2, method = "spearman", use = "pairwise.complete.obs"), digits = 2))))
  lines(xplot, pop2, type = "l", lwd = 1.5, col = pop2Col)
  axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"))
  mtext(side = 1, line = 2.25, cex = 1, text = "Coordinates (Mb)")
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, text = Ylabel, col = "black")
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = centromeres, lty = 3, lwd = 2, col = "black")
  box(lwd = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c(pop1Col, pop2Col),
         text.col = c(pop1Col, pop2Col),
         text.font = c(1, 1),
         ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
}

# Plot genome profiles of U-test P-values and 
# population mean crossovers per window
pdf(file = paste0(UtestPlotDir,
                  "MannWhitneyWilcoxon_Pvals_genomeProfiles_",
                  pop1Name, "_vs_", pop2Name, "_", winName, "_", FDRname, ".pdf"),
    height = 8, width = 10)
par(mfrow = c(2, 1))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))

# Plot P-values
minusLog10PvalPlot(xplot = pop1_indCOs_perWindow_list[[1]]$cumWindow,
                   Pvals = UtestPvals,
                   AdjPvals = UtestAdjPvals,
                   PvalsCol = "dodgerblue3",
                   AdjPvalsCol = "red")
# Plot population mean crossovers per window
pop1Vpop2GenomePlot(xplot = pop1_indCOs_perWindow_list[[1]]$cumWindow,
                    pop1 = pop1_winIndCOs_means,
                    pop2 = pop2_winIndCOs_means,
                    Ylabel = "Population mean crossovers",
                    legendLoc = "topleft",
                    legendLabs = c(pop1Name, pop2Name),
                    pop1Col = "grey50",
                    pop2Col = "forestgreen")
dev.off()
print(paste0("MannWhitneyWilcoxon_Pvals_genomeProfiles_", pop1Name, "_vs_", pop2Name, "_", winName, "_", FDRname, ".pdf written to ", UtestPlotDir))

# Plot genome profiles of normal linear regression P-values and 
# population mean crossovers per window
pdf(file = paste0(plotDir,
                  "NormalLinearRegression_Pvals_genomeProfiles_",
                  pop1Name, "_vs_", pop2Name, "_", winName, "_", FDRname, ".pdf"),
    height = 8, width = 10)
par(mfrow = c(2, 1))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))

# Plot P-values
minusLog10PvalPlot(xplot = pop1_indCOs_perWindow_list[[1]]$cumWindow,
                   Pvals = normal_Pvals,
                   AdjPvals = normal_AdjPvals,
                   PvalsCol = "dodgerblue3",
                   AdjPvalsCol = "red")
# Plot population mean crossovers per window
pop1Vpop2GenomePlot(xplot = pop1_indCOs_perWindow_list[[1]]$cumWindow,
                    pop1 = pop1_winIndCOs_means,
                    pop2 = pop2_winIndCOs_means,
                    Ylabel = "Population mean crossovers",
                    legendLoc = "topleft",
                    legendLabs = c(pop1Name, pop2Name),
                    pop1Col = "grey50",
                    pop2Col = "forestgreen")
dev.off()
print(paste0("NormalLinearRegression_Pvals_genomeProfiles_", pop1Name, "_vs_", pop2Name, "_", winName, "_", FDRname, ".pdf written to ", plotDir))

# Plot genome profiles of Poisson regression P-values and 
# population mean crossovers per window
pdf(file = paste0(plotDir,
                  "PoissonRegression_Pvals_genomeProfiles_",
                  pop1Name, "_vs_", pop2Name, "_", winName, "_", FDRname, ".pdf"),
    height = 8, width = 10)
par(mfrow = c(2, 1))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))

# Plot P-values
minusLog10PvalPlot(xplot = pop1_indCOs_perWindow_list[[1]]$cumWindow,
                   Pvals = poisson_Pvals,
                   AdjPvals = poisson_AdjPvals,
                   PvalsCol = "dodgerblue3",
                   AdjPvalsCol = "red")
# Plot population mean crossovers per window
pop1Vpop2GenomePlot(xplot = pop1_indCOs_perWindow_list[[1]]$cumWindow,
                    pop1 = pop1_winIndCOs_means,
                    pop2 = pop2_winIndCOs_means,
                    Ylabel = "Population mean crossovers",
                    legendLoc = "topleft",
                    legendLabs = c(pop1Name, pop2Name),
                    pop1Col = "grey50",
                    pop2Col = "forestgreen")
dev.off()
print(paste0("PoissonRegression_Pvals_genomeProfiles_", pop1Name, "_vs_", pop2Name, "_", winName, "_", FDRname, ".pdf written to ", plotDir))
