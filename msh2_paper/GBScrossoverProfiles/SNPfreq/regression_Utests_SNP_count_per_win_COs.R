#!/applications/R/R-3.5.0/bin/Rscript

# Use regression models and U-tests to determine whether significant
# differences exist between genotypes with regard to SNP frequencies
# within windows along crossover loci, and within genotypes between
# crossover loci and random loci

# Resources consulted:
# http://datavoreconsulting.com/programming-tips/count-data-glms-choosing-poisson-negative-binomial-zero-inflated-poisson/
# https://www.physiology.org/doi/full/10.1152/advan.00017.2010?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%3dpubmed
# https://stats.idre.ucla.edu/r/dae/poisson-regression/
# https://stats.idre.ucla.edu/r/dae/zip/
# https://stats.idre.ucla.edu/r/dae/negative-binomial-regression/
# https://stats.idre.ucla.edu/r/dae/zinb/
# https://data.princeton.edu/wws509/r/overdispersion

# Usage from within directory containing windowed CO counts files:
# ./regression_Utests_SNP_count_per_win_COs.R coller.filtarb coller.filtmsh2 5000 5kb 50 0.1

# Some of these packages may not already be installed
# e.g., install.packages("pscl")
library(GenomicRanges)
library(parallel)
library(MASS) # glm.nb() included
library(pscl) # zeroinfl() included

#pop1Name <- "coller.filtarb"
#pop2Name <- "coller.filtmsh2"
#flankSize <- 5000
#flankName <- "5kb"
#winSize <- 50
#FDR <- 0.1
#FDRname <- paste0("FDR", as.character(FDR))

args <- commandArgs(trailingOnly = TRUE)
pop1Name <- args[1]
pop2Name <- args[2]
flankSize <- as.numeric(args[3])
flankName <- as.character(args[4])
winSize <- as.numeric(args[5])
FDR <- as.numeric(args[6])
FDRname <- paste0("FDR", as.character(FDR))

matDir <- "matrices/"

# Load matrices in which each row corresponds to a crossover or random locus
# and each column corresponds to a window within that locus
pop1_matList <- list(read.table(paste0(matDir,
                                       pop1Name, "_COs_SNP_frequency_feature_target_and_",
                                       flankName, "_flank_dataframe.txt"),
                                header = T),
                     read.table(paste0(matDir,
                                       pop1Name, "_COs_SNP_frequency_ranLoc_target_and_",
                                       flankName, "_flank_dataframe.txt"),
                                header = T))
pop2_matList <- list(read.table(paste0(matDir,
                                       pop2Name, "_COs_SNP_frequency_feature_target_and_",
                                       flankName, "_flank_dataframe.txt"),
                                header = T),
                     read.table(paste0(matDir,
                                       pop2Name, "_COs_SNP_frequency_ranLoc_target_and_",
                                       flankName, "_flank_dataframe.txt"),
                                header = T))

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


# For each population, calculate mean and variances of SNP frequencies
# in each window along crossovers (list element [[1]]) and
# random loci (list element [[2]])
pop1_means <- list(as.vector(apply(X = pop1_matList[[1]],
                                   MARGIN = 2,
                                   FUN = mean)),
                   as.vector(apply(X = pop1_matList[[2]],
                                   MARGIN = 2,
                                   FUN = mean))) 
pop2_means <- list(as.vector(apply(X = pop2_matList[[1]],
                                   MARGIN = 2,
                                   FUN = mean)),
                   as.vector(apply(X = pop2_matList[[2]],
                                   MARGIN = 2,
                                   FUN = mean))) 
pop1_vars <- list(as.vector(apply(X = pop1_matList[[1]],
                                   MARGIN = 2,
                                   FUN = var)),
                   as.vector(apply(X = pop1_matList[[2]],
                                   MARGIN = 2,
                                   FUN = var))) 
pop2_vars <- list(as.vector(apply(X = pop2_matList[[1]],
                                  MARGIN = 2,
                                  FUN = var)),
                  as.vector(apply(X = pop2_matList[[2]],
                                  MARGIN = 2,
                                  FUN = var))) 

# For pop1, create a list in which each list element (x) is a vector of
# F2 crossover interval counts for one proportionally scaled window,
# summed across all chromosome arms
pop1_winIndCOs_list <- mclapply(seq_along(1:dim(pop1_indCOs_perWindow_list[[1]])[1]), function(x) {
  windowedCOs <- NULL
  for(y in seq_along(pop1_indCOs_perWindow_list)) {
    windowedCOs <- c(windowedCOs, sum(pop1_indCOs_perWindow_list[[y]][x,]))
  }
  windowedCOs
}, mc.cores = dim(pop1_indCOs_perWindow_list[[1]])[1])
# For pop2, create a list in which each list element (x) is a vector of
# F2 crossover interval counts for one window
pop2_winIndCOs_list <- mclapply(seq_along(1:dim(pop2_indCOs_perWindow_list[[1]])[1]), function(x) {
  windowedCOs <- NULL
  for(y in seq_along(pop2_indCOs_perWindow_list)) {
    windowedCOs <- c(windowedCOs, sum(pop2_indCOs_perWindow_list[[y]][x,]))
  }
  windowedCOs
}, mc.cores = dim(pop2_indCOs_perWindow_list[[1]])[1])

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

# Plot means vs variances as a quick check for overdispersion;
# data look overdispersed (greater variances than means)
pdf(paste0(outDir, pop1Name, "_", pop2Name,
           "_crossover_and_ranLoc_SNP_means_vs_variances.pdf"),
    height = 10, width = 10)
par(mfrow = c(2, 2))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(x = pop1_means[[1]], y = pop1_vars[[1]], pch = 19,
     xlab = "Mean", ylab = "Variance",
     xlim = c(min(pop1_means[[1]], pop1_vars[[1]]),
              max(pop1_means[[1]], pop1_vars[[1]])),
     ylim = c(min(pop1_means[[1]], pop1_vars[[1]]),
              max(pop2_means[[1]], pop1_vars[[1]])),
     main = paste0("SNPs around ", pop1Name, " crossovers"),
     cex.main = 0.8)
plot(x = pop2_means[[1]], y = pop2_vars[[1]], pch = 19,
     xlab = "Mean", ylab = "Variance",
     xlim = c(min(pop2_means[[1]], pop2_vars[[1]]),
              max(pop2_means[[1]], pop2_vars[[1]])),
     ylim = c(min(pop2_means[[1]], pop2_vars[[1]]),
              max(pop2_means[[1]], pop2_vars[[1]])),
     main = paste0("SNPs around ", pop2Name, " crossovers"),
     cex.main = 0.8)
plot(x = pop1_means[[2]], y = pop1_vars[[2]], pch = 19,
     xlab = "Mean", ylab = "Variance",
     xlim = c(min(pop1_means[[2]], pop1_vars[[2]]),
              max(pop1_means[[2]], pop1_vars[[2]])),
     ylim = c(min(pop1_means[[2]], pop1_vars[[2]]),
              max(pop1_means[[2]], pop1_vars[[2]])),
     main = paste0("SNPs around ", pop1Name, " random loci"),
     cex.main = 0.8)
plot(x = pop2_means[[2]], y = pop2_vars[[2]], pch = 19,
     xlab = "Mean", ylab = "Variance",
     xlim = c(min(pop2_means[[2]], pop2_vars[[2]]),
              max(pop2_means[[2]], pop2_vars[[2]])),
     ylim = c(min(pop2_means[[2]], pop2_vars[[2]]),
              max(pop2_means[[2]], pop2_vars[[2]])),
     main = paste0("SNPs around ", pop2Name, " random loci"),
     cex.main = 0.8)
dev.off()


# Create a factor with two levels (pop1Name and pop2Name, corresponding to genotype)
# to be used as the predictor in regression models
genotype <- factor(c(rep(pop1Name, times = dim(pop1_matList[[1]])[1]),
                     rep(pop2Name, times = dim(pop2_matList[[1]])[1])),
                   levels = c(pop1Name, pop2Name))

# For each population, create a factor with two levels ("Crossovers" and "Random loci",
# corresponding to feature) to be used as the predictor in regression models
pop1_feature <- factor(c(rep("Crossovers", times = dim(pop1_matList[[1]])[1]),
                         rep("Random loci", times = dim(pop1_matList[[2]])[1])),
                       levels = c("Crossovers", "Random loci"))
pop2_feature <- factor(c(rep("Crossovers", times = dim(pop2_matList[[1]])[1]),
                         rep("Random loci", times = dim(pop2_matList[[2]])[1])),
                       levels = c("Crossovers", "Random loci"))


## By genotype

# For each  window, fit various regression models in which genotype is the
# predictor and SNP frequency is the response variable
# Determine which model has best fit

# Normal linear regression, equivalent to a t-test:
genotype_normal <- lapply(seq_along(pop1_matList[[1]]), function(x) {
  glm(formula = c(pop1_matList[[1]][[x]], pop2_matList[[1]][[x]]) ~ genotype,
      family = gaussian(link = "identity"))
})
# Get estimates
genotype_normal_estimates <- lapply(seq_along(genotype_normal), function(x) {
  coef(genotype_normal[[x]])
})
# Get coefficients
genotype_normal_coef <- lapply(seq_along(genotype_normal), function(x) {
  coef(summary(genotype_normal[[x]]))
})
# Get P-values
genotype_normal_Pvals <- sapply(seq_along(genotype_normal), function(x) {
  coef(summary(genotype_normal[[x]]))[8]
})
# Correct for multiple testing
genotype_normal_AdjPvals <- p.adjust(p = genotype_normal_Pvals, method = "BH")
# Get Bayesian Information Criterion (BIC); the smaller, the better the fit
genotype_normal_BIC <- sapply(seq_along(genotype_normal), function(x) {
  if( !is.infinite(BIC(genotype_normal[[x]])) ) {
    BIC(genotype_normal[[x]])
  } else {
    NA
  }
})
# Get Akaike Information Criterion (AIC); the smaller, the better the fit
genotype_normal_AIC <- sapply(seq_along(genotype_normal), function(x) {
  if( !is.infinite(AIC(genotype_normal[[x]], k = 2)) ) {
    AIC(genotype_normal[[x]], k = 2)
  } else {
    NA
  }
})
# Get 95% confidence intervals
genotype_normal_CIs <- lapply(seq_along(genotype_normal), function(x) {
  confint(genotype_normal[[x]])
})

# Evaluate model goodness-of-fit using chi-squared test based on the
# residual deviance and degrees of freedom
# A P-value > 0.05 indicates that the model fits the data 
genotype_normal_chisqPvals <- sapply(seq_along(genotype_normal), function(x) {
  1 - pchisq(summary(genotype_normal[[x]])$deviance,
             summary(genotype_normal[[x]])$df.residual)
})
print(paste0("P-values from chi-squared goodness-of-fit tests evaluating genotype_normal linear regression models for ",
             length(genotype_normal), " genomic windows:"))
print(genotype_normal_chisqPvals)
print(paste0("Sum of P-values from chi-squared goodness-of-fit tests evaluating genotype_normal linear regression models for ",
             length(genotype_normal), " genomic windows:"))
print(sum(genotype_normal_chisqPvals))

# Use the model to predict mean SNP frequencies for each genotype
# and corresponding standard errors
# The model predicts the correct mean counts
genotype_normal_predict <- lapply(seq_along(genotype_normal), function(x) {
  cbind(data.frame(genotype = c(pop1Name, pop2Name)),
        mean = predict(genotype_normal[[x]],
                       newdata = data.frame(genotype = c(pop1Name, pop2Name)),
                       type = "response"),
        SE = predict(genotype_normal[[x]],
                     newdata = data.frame(genotype = c(pop1Name, pop2Name)),
                     type = "response",
                     se.fit = TRUE)$se.fit
       )
})
genotype_normal_correct_predictions <- sapply(seq_along(genotype_normal_predict), function(x) {
  c(print(all.equal(genotype_normal_predict[x][[1]]$mean[1],
                    pop1_means[[1]][x])),
    print(all.equal(genotype_normal_predict[x][[1]]$mean[2],
                    pop2_means[[1]][x])))
})
print(paste0(as.character(sum(genotype_normal_correct_predictions)),
             "/", as.character(length(genotype_normal_correct_predictions))))


# Poisson regression:
genotype_poisson <- lapply(seq_along(pop1_matList[[1]]), function(x) {
  glm(formula = c(pop1_matList[[1]][[x]], pop2_matList[[1]][[x]]) ~ genotype,
      family = poisson(link = "log"))
})
# Get estimates
genotype_poisson_estimates <- lapply(seq_along(genotype_poisson), function(x) {
  coef(genotype_poisson[[x]])
})
# Get coefficients
genotype_poisson_coef <- lapply(seq_along(genotype_poisson), function(x) {
    coef(summary(genotype_poisson[[x]]))
})
# Get P-values
genotype_poisson_Pvals <- sapply(seq_along(genotype_poisson), function(x) {
  coef(summary(genotype_poisson[[x]]))[8]
})
# Correct for multiple testing
genotype_poisson_AdjPvals <- p.adjust(p = genotype_poisson_Pvals, method = "BH")
# Get Bayesian Information Criterion (BIC); the smaller, the better the fit
genotype_poisson_BIC <- sapply(seq_along(genotype_poisson), function(x) {
  if( !is.infinite(BIC(genotype_poisson[[x]])) ) {
    BIC(genotype_poisson[[x]])
  } else {
    NA
  }
})
# Get Akaike Information Criterion (AIC); the smaller, the better the fit
genotype_poisson_AIC <- sapply(seq_along(genotype_poisson), function(x) {
  if( !is.infinite(AIC(genotype_poisson[[x]], k = 2)) ) {
    AIC(genotype_poisson[[x]], k = 2)
  } else {
    NA
  }
})
# Get 95% confidence intervals
genotype_poisson_CIs <- lapply(seq_along(genotype_poisson), function(x) {
  confint(genotype_poisson[[x]])
})

# Evaluate model goodness-of-fit using chi-squared test based on the
# residual deviance and degrees of freedom
# A P-value > 0.05 indicates that the model fits the data 
genotype_poisson_chisqPvals <- sapply(seq_along(genotype_poisson), function(x) {
  1 - pchisq(summary(genotype_poisson[[x]])$deviance,
             summary(genotype_poisson[[x]])$df.residual)
})
print(paste0("P-values from chi-squared goodness-of-fit tests evaluating Poisson regression models for ",
             length(genotype_poisson), " genomic windows:"))
print(genotype_poisson_chisqPvals)
print(paste0("Sum of P-values from chi-squared goodness-of-fit tests evaluating Poisson regression models for ",
             length(genotype_poisson), " genomic windows:"))
print(sum(genotype_poisson_chisqPvals))

# Use the model to predict mean crossover counts for each genotype
# and corresponding standard errors
# The model predicts the correct mean counts
genotype_poisson_predict <- lapply(seq_along(genotype_poisson), function(x) {
  cbind(data.frame(genotype = c(pop1Name, pop2Name)),
        mean = predict(genotype_poisson[[x]],
                       newdata = data.frame(genotype = c(pop1Name, pop2Name)),
                       type = "response"),
        SE = predict(genotype_poisson[[x]],
                     newdata = data.frame(genotype = c(pop1Name, pop2Name)),
                     type = "response",
                     se.fit = TRUE)$se.fit
       )
})
genotype_poisson_correct_predictions <- sapply(seq_along(genotype_poisson_predict), function(x) {
  c(print(all.equal(genotype_poisson_predict[x][[1]]$mean[1],
                    pop1_means[[1]][x])),
    print(all.equal(genotype_poisson_predict[x][[1]]$mean[2],
                    pop2_means[[1]][x])))
})
print(paste0(as.character(sum(genotype_poisson_correct_predictions)),
             "/", as.character(length(genotype_poisson_correct_predictions))))


# Zero-inflated Poisson (ZIP) regression:

# Check if ordinary Poisson model underestimates the probability of
# zero SNPs in each window
genotype_zobs <- sapply(seq_along(pop1_matList[[1]]), function(x) {
  mean(c(pop1_matList[[1]][[x]], pop2_matList[[1]][[x]]) == 0)
})
genotype_poisson_zpr_means <- sapply(seq_along(genotype_poisson), function(x) {
  mean(dpois(0, exp(predict(genotype_poisson[[x]]))))
})
print(paste0("genotype_poisson models predict fewer than observed occurrences of zero SNPs in a window for ",
             as.character(sum(genotype_poisson_zpr_means < genotype_zobs)),
             "/", length(genotype_zobs), " windows"))

# Make ZIP model
genotype_ZIP <- lapply(seq_along(pop1_matList[[1]]), function(x) {
  zeroinfl(formula = c(pop1_matList[[1]][[x]], pop2_matList[[1]][[x]]) ~ genotype | 1,
           dist = "poisson",
           maxit = 2000)
})
print("Representative example: ZIP model for first genomic window")
print(summary(genotype_ZIP[[1]]))
# Get estimates
genotype_ZIP_estimates <- lapply(seq_along(genotype_ZIP), function(x) {
  coef(genotype_ZIP[[x]])
})
# Get coefficients
genotype_ZIP_coef <- lapply(seq_along(genotype_ZIP), function(x) {
    coef(summary(genotype_ZIP[[x]]))
})
# Get P-values
genotype_ZIP_Pvals <- sapply(seq_along(genotype_ZIP), function(x) {
  coef(summary(genotype_ZIP[[x]]))$count[8]
})
# Correct for multiple testing
genotype_ZIP_AdjPvals <- p.adjust(p = genotype_ZIP_Pvals, method = "BH")
# Get Bayesian Information Criterion (BIC); the smaller, the better the fit
genotype_ZIP_BIC <- sapply(seq_along(genotype_ZIP), function(x) {
  if( !is.infinite(AIC(genotype_ZIP[[x]], k = log(length(pop1_matList[[1]][[x]])))) ) {
    AIC(genotype_ZIP[[x]], k = log(length(pop1_matList[[1]][[x]])))
  } else {
    NA
  }
})
# Get Akaike Information Criterion (AIC); the smaller, the better the fit
genotype_ZIP_AIC <- sapply(seq_along(genotype_ZIP), function(x) {
  if( !is.infinite(AIC(genotype_ZIP[[x]], k = 2)) ) {
    AIC(genotype_ZIP[[x]], k = 2)
  } else {
    NA
  }
})
# Get 95% confidence intervals
genotype_ZIP_CIs <- lapply(seq_along(genotype_ZIP), function(x) {
  confint(genotype_ZIP[[x]])
})

# Use the model to predict mean crossover counts for each genotype
# and corresponding standard errors
# The model predicts the correct mean counts
genotype_ZIP_predict <- lapply(seq_along(genotype_ZIP), function(x) {
  cbind(data.frame(genotype = c(pop1Name, pop2Name)),
        mean = predict(genotype_ZIP[[x]],
                       newdata = data.frame(genotype = c(pop1Name, pop2Name)),
                       type = "response"),
        SE = predict(genotype_ZIP[[x]],
                     newdata = data.frame(genotype = c(pop1Name, pop2Name)),
                     type = "response",
                     se.fit = TRUE)$se.fit
       )
})
genotype_ZIP_correct_predictions <- sapply(seq_along(genotype_ZIP_predict), function(x) {
  c(print(all.equal(genotype_ZIP_predict[x][[1]]$mean[1],
                    pop1_means[[1]][x])),
    print(all.equal(genotype_ZIP_predict[x][[1]]$mean[2],
                    pop2_means[[1]][x])))
})
print(paste0(as.character(sum(genotype_ZIP_correct_predictions)),
             "/", as.character(length(genotype_ZIP_correct_predictions))))


# Evaluate whether model solves the problem of excess zeroes by
# predicting π and μ, and calculate the combined probability of no SNPs
# Use the predict() options "zero" and "count" to obtain π and μ
genotype_ZIP_zero <- lapply(seq_along(genotype_ZIP), function(x) {
  predict(genotype_ZIP[[x]], type = "zero") #  π
})
genotype_ZIP_count <- lapply(seq_along(genotype_ZIP), function(x) {
  predict(genotype_ZIP[[x]], type = "count") # μ
})
genotype_ZIP_zero_count_means <- sapply(seq_along(genotype_ZIP), function(x) {
  mean( genotype_ZIP_zero[[x]] + 
         ( (1-genotype_ZIP_zero[[x]]) * exp(-genotype_ZIP_count[[x]]) )
  )
})
print(paste0("genotype_ZIP models predict fewer than observed occurrences of zero SNPs in a window for ",
             as.character(sum(genotype_ZIP_zero_count_means < genotype_zobs)),
             "/", length(genotype_zobs), " windows"))
# ZIP model accurately estimates the probability of zero SNPs in each window

# Compare zero-inflated Poisson with ordinary Poisson regression model
# using the Vuong test
genotype_ZIP_vuong <- lapply(seq_along(genotype_ZIP), function(x) {
  capture.output(vuong(m1 = genotype_poisson[[x]], m2 = genotype_ZIP[[x]]))
})
 


## Model comparisons
# Determine if Poisson BIC is lower than normal BIC (indicating better fit)
genotype_poisson_BIC_vs_genotype_normal_BIC <- sapply(seq_along(genotype_poisson), function(x) {
  genotype_poisson_BIC[x] < genotype_normal_BIC[x]
})
# Determine if genotype_normal BIC is NA
genotype_normal_BIC_isNA <- sapply(seq_along(genotype_normal), function(x) {
  is.na(genotype_normal_BIC[x])
})
print(paste0("Poisson BIC is < Normal BIC for ",
             sum(genotype_poisson_BIC_vs_genotype_normal_BIC, na.rm = T), "/", length(genotype_poisson), " windows"))
print(paste0("Normal BIC is NA for ",
             sum(genotype_normal_BIC_isNA), "/", length(genotype_normal), " scaled windows"))

# Determine if Poisson BIC is lower than ZIP BIC (indicating better fit)
genotype_poisson_BIC_vs_genotype_ZIP_BIC <- sapply(seq_along(genotype_poisson), function(x) {
  genotype_poisson_BIC[x] < genotype_ZIP_BIC[x]
})
print(paste0("Poisson BIC is < genotype_ZIP BIC for ",
             sum(genotype_poisson_BIC_vs_genotype_ZIP_BIC, na.rm = T), "/", length(genotype_poisson), " windows"))

# Note: negative binomial models should not be fit to these data
# due to absence of overdispersion 
## N
#negbin <- lapply(seq_along(pop1_matList[[1]]), function(x) {
#  glm.nb(formula = pops_winIndCOs_list[[x]] ~ genotype,
#         control = glm.control(maxit = 2000))
#})
## Zero-inflated negative binomial regression:
#ZINB <- lapply(seq_along(pop1_matList[[1]]), function(x) {
#  zeroinfl(formula = pops_winIndCOs_list[[x]] ~ genotype | 1,
#           dist = "negbin")
#})

# Plot genotype_poisson_BIC vs genotype_normal_BIC
pdf(paste0(outDir, pop1Name, "_", pop2Name,
           "_Poisson_vs_normal_BICs_for_crossovers_",
           propName, "_of_chromosome_arms.pdf"),
    height = 5, width = 5)
par(mfrow = c(1, 1))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(x = 1:length(genotype_poisson_BIC), y = genotype_poisson_BIC, col = "red", type = "p", pch = 19,
     xlab = paste0("Scaled windows (", propName, ")"),
     ylab = "BIC",
     ylim = c(min(genotype_poisson_BIC, genotype_normal_BIC),
              max(genotype_poisson_BIC, genotype_normal_BIC)))
lines(x = 1:length(genotype_poisson_BIC), y = genotype_normal_BIC, col = "blue", type = "p", pch = 19)
legend("top",
       legend = c("Poisson regression", "Normal linear regression"),
       col = c("red", "blue"),
       text.col = c("red", "blue"),
       text.font = c(1, 1),
       ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
dev.off()

# Plot genotype_poisson_AIC vs genotype_normal_AIC
pdf(paste0(outDir, pop1Name, "_", pop2Name,
           "_Poisson_vs_normal_AICs_for_crossovers_",
           propName, "_of_chromosome_arms.pdf"),
    height = 5, width = 5)
par(mfrow = c(1, 1))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(x = 1:length(genotype_poisson_AIC), y = genotype_poisson_AIC, col = "red", type = "p", pch = 19,
     xlab = paste0("Scaled windows (", propName, ")"),
     ylab = "AIC",
     ylim = c(min(genotype_poisson_AIC, genotype_normal_AIC),
              max(genotype_poisson_AIC, genotype_normal_AIC)))
lines(x = 1:length(genotype_poisson_AIC), y = genotype_normal_AIC, col = "blue", type = "p", pch = 19)
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


# Function to plot TEL-CEN profile of -log10-transformed P-values and adjusted P-values
minusLog10PvalPlot <- function(xplot,
                               Pvals, PvalsCol,
                               AdjPvals, AdjPvalsCol) {
  plot(x = xplot, y = -log10(Pvals), col = PvalsCol, type = "l", lwd = 1.5,
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
       ylim = c(0,
                pmax(-log10(FDR), max(-log10(AdjPvals), na.rm = T))),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(pop1Name)*" versus "*.(pop2Name)*" crossover counts"),
       cex.main = 1)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1, adj = c(0.5, -3.0), xpd = NA, srt = -90, col = AdjPvalsCol,
       labels = bquote("-Log"[10]*"(BH-adjusted "*italic("P")*"-value)"))
  axis(side = 4, cex.axis = 1, lwd.tick = 1.5)
  if(prop <= 15) {
    axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
         at = c(1, seq(2, prop, by = 1)),
         labels = c(expression(italic("TEL")),
                    seq(2, prop-1, by = 1),
                    expression(italic("CEN"))))
  } else {
    axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
         at = c(1, pretty(1:length(xplot))[2:(length(pretty(1:length(xplot)))-1)], prop),
         labels = c(expression(italic("TEL")),
                    pretty(1:length(xplot))[2:(length(pretty(1:length(xplot)))-1)],
                    expression(italic("CEN"))))
  }
  mtext(side = 1, line = 2.25, cex = 1, text = paste0("Scaled windows (", propName, ")"))
  abline(h = -log10(FDR), lty = 5, lwd = 1, col = AdjPvalsCol)

  box(lwd = 1.5)
}

# Function to plot TEL-CEN profile of windowed pop1 vs pop2 crossovers (one Y-axis) 
pop1Vpop2GenomePlot <- function(xplot,
                                pop1, pop1Col,
                                pop2, pop2Col,
                                Ylabel,
                                legendLoc,
                                legendLabs) {
  plot(xplot, pop1, type = "l", lwd = 1.5, col = pop1Col,
       ylim = c(0,
                max(c(pop1, pop2))),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(italic("r"[s]) ~ " = " ~ .(round(cor(pop1, pop2, method = "spearman", use = "pairwise.complete.obs"), digits = 2))))
  lines(xplot, pop2, type = "l", lwd = 1.5, col = pop2Col)
  if(prop <= 15) {
    axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
         at = c(1, seq(2, prop, by = 1)),
         labels = c(expression(italic("TEL")),
                    seq(2, prop-1, by = 1),
                    expression(italic("CEN"))))
  } else {
    axis(side = 1, cex.axis = 1, lwd.tick = 1.5,
         at = c(1, pretty(1:length(xplot))[2:(length(pretty(1:length(xplot)))-1)], prop),
         labels = c(expression(italic("TEL")),
                    pretty(1:length(xplot))[2:(length(pretty(1:length(xplot)))-1)],
                    expression(italic("CEN"))))
  }
  mtext(side = 1, line = 2.25, cex = 1, text = paste0("Scaled windows (", propName, ")"))
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2.25, cex = 1, text = Ylabel, col = "black")
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
                  "MannWhitneyWilcoxon_Pvals_TelCenProfiles_",
                  pop1Name, "_vs_", pop2Name, "_", propName, "_", FDRname, ".pdf"),
    height = 7, width = 7)
par(mfrow = c(2, 1))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
# Plot P-values
minusLog10PvalPlot(xplot = 1:length(pop1_winIndCOs_means),
                   Pvals = UtestPvals,
                   AdjPvals = UtestAdjPvals,
                   PvalsCol = "dodgerblue3",
                   AdjPvalsCol = "red")
# Plot population mean crossovers per window
pop1Vpop2GenomePlot(xplot = 1:length(pop1_winIndCOs_means),
                    pop1 = pop1_winIndCOs_means,
                    pop2 = pop2_winIndCOs_means,
                    Ylabel = "Population mean crossovers",
                    legendLoc = "top",
                    legendLabs = c(pop1Name, pop2Name),
                    pop1Col = "grey50",
                    pop2Col = "forestgreen")
dev.off()
print(paste0("MannWhitneyWilcoxon_Pvals_TelCenProfiles_", pop1Name, "_vs_", pop2Name, "_", propName, "_", FDRname, ".pdf written to ", UtestPlotDir))

# Plot genome profiles of genotype_normal linear regression P-values and 
# population mean crossovers per window
pdf(file = paste0(plotDir,
                  "NormalLinearRegression_Pvals_TelCenProfiles_",
                  pop1Name, "_vs_", pop2Name, "_", propName, "_", FDRname, ".pdf"),
    height = 7, width = 7)
par(mfrow = c(2, 1))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
# Plot P-values
minusLog10PvalPlot(xplot = 1:length(pop1_winIndCOs_means),
                   Pvals = genotype_normal_Pvals,
                   AdjPvals = genotype_normal_AdjPvals,
                   PvalsCol = "dodgerblue3",
                   AdjPvalsCol = "red")
# Plot population mean crossovers per window
pop1Vpop2GenomePlot(xplot = 1:length(pop1_winIndCOs_means),
                    pop1 = pop1_winIndCOs_means,
                    pop2 = pop2_winIndCOs_means,
                    Ylabel = "Population mean crossovers",
                    legendLoc = "top",
                    legendLabs = c(pop1Name, pop2Name),
                    pop1Col = "grey50",
                    pop2Col = "forestgreen")
dev.off()
print(paste0("NormalLinearRegression_Pvals_TelCenProfiles_", pop1Name, "_vs_", pop2Name, "_", propName, "_", FDRname, ".pdf written to ", plotDir))

# Plot genome profiles of Poisson regression P-values and 
# population mean crossovers per window
pdf(file = paste0(plotDir,
                  "PoissonRegression_Pvals_TelCenProfiles_",
                  pop1Name, "_vs_", pop2Name, "_", propName, "_", FDRname, ".pdf"),
    height = 7, width = 7)
par(mfrow = c(2, 1))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
# Plot P-values
minusLog10PvalPlot(xplot = 1:length(pop1_winIndCOs_means),
                   Pvals = genotype_poisson_Pvals,
                   AdjPvals = genotype_poisson_AdjPvals,
                   PvalsCol = "dodgerblue3",
                   AdjPvalsCol = "red")
# Plot population mean crossovers per window
pop1Vpop2GenomePlot(xplot = 1:length(pop1_winIndCOs_means),
                    pop1 = pop1_winIndCOs_means,
                    pop2 = pop2_winIndCOs_means,
                    Ylabel = "Population mean crossovers",
                    legendLoc = "top",
                    legendLabs = c(pop1Name, pop2Name),
                    pop1Col = "grey50",
                    pop2Col = "forestgreen")
dev.off()
