#!/applications/R/R-3.5.0/bin/Rscript

# Use regression models to determine whether significant
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
# ./regression_SNP_count_per_win_COs.R coller.filtarb coller.filtmsh2 5000 5kb 200 200bp 0.1

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
#winSize <- 200
#winName <- "200bp"
#FDR <- 0.1
#FDRname <- paste0("FDR", as.character(FDR))

args <- commandArgs(trailingOnly = TRUE)
pop1Name <- args[1]
pop2Name <- args[2]
flankSize <- as.numeric(args[3])
flankName <- as.character(args[4])
winSize <- as.numeric(args[5])
winName <- as.character(args[6])
FDR <- as.numeric(args[7])
FDRname <- paste0("FDR", as.character(FDR))

matDir <- "matrices/"

# Load matrices in which each row corresponds to a crossover or random locus
# and each column corresponds to a window within that locus
pop1_matList <- list(read.table(paste0(matDir,
                                       pop1Name, "_COs_SNP_frequency_feature_target_and_",
                                       flankName, "_flank_", winName, "_win_dataframe.txt"),
                                header = T),
                     read.table(paste0(matDir,
                                       pop1Name, "_COs_SNP_frequency_ranLoc_target_and_",
                                       flankName, "_flank_", winName, "_win_dataframe.txt"),
                                header = T))
pop2_matList <- list(read.table(paste0(matDir,
                                       pop2Name, "_COs_SNP_frequency_feature_target_and_",
                                       flankName, "_flank_", winName, "_win_dataframe.txt"),
                                header = T),
                     read.table(paste0(matDir,
                                       pop2Name, "_COs_SNP_frequency_ranLoc_target_and_",
                                       flankName, "_flank_", winName, "_win_dataframe.txt"),
                                header = T))

# Define and create new directories and subdirectories
# to contain results
sizeDir <- paste0(flankName, "_flank_", winName, "_win/")
outDir <- paste0(sizeDir, "regression_models/")
plotDir <- paste0(outDir, FDRname, "/")
system(paste0("[ -d ", sizeDir, " ] || mkdir ", sizeDir))
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

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

# Plot means vs variances as a quick check for overdispersion;
# data look overdispersed (greater variances than means)
pdf(paste0(outDir, pop1Name, "_", pop2Name,
           "_crossover_and_ranLoc_SNP_means_vs_variances_",
           flankName, "_flank_", winName, "_win.pdf"),
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

# Get R2 
genotype_normal_R2 <- sapply(seq_along(genotype_normal), function(x) {
  1 - (genotype_normal[[x]]$deviance/genotype_normal[[x]]$null.deviance)
})
# Get Cohen's F2 as a measure of effect size for linear regression
genotype_normal_F2 <- sapply(seq_along(genotype_normal), function(x) {
  genotype_normal_R2[x] / (1 - genotype_normal_R2[x])
})


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

# Get R2 
genotype_poisson_R2 <- sapply(seq_along(genotype_poisson), function(x) {
  1 - (genotype_poisson[[x]]$deviance/genotype_poisson[[x]]$null.deviance)
})
# Get Cohen's F2 as a measure of effect size for linear regression
genotype_poisson_F2 <- sapply(seq_along(genotype_poisson), function(x) {
  genotype_poisson_R2[x] / (1 - genotype_poisson_R2[x])
})


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

# Make model
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
                       type = "response")
         ## Doesn't work for this type of model for some reason
#        SE = predict(genotype_ZIP[[x]],
#                     newdata = data.frame(genotype = c(pop1Name, pop2Name)),
#                     type = "response",
#                     se.fit = TRUE)$se.fit
       )
})
genotype_ZIP_correct_predictions <- sapply(seq_along(genotype_ZIP_predict), function(x) {
  c(print(all.equal(genotype_ZIP_predict[x][[1]]$mean[1],
                    pop1_means[[1]][x])),
    print(all.equal(genotype_ZIP_predict[x][[1]]$mean[2],
                    pop2_means[[1]][x])))
})
#print(paste0(as.character(sum(genotype_ZIP_correct_predictions)),
#             "/", as.character(length(genotype_ZIP_correct_predictions))))


# Evaluate whether ZIP model solves the problem of excess zeroes by
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

# Compare ZIP with ordinary Poisson regression model
# using the Vuong test
genotype_ZIP_vuong <- lapply(seq_along(genotype_ZIP), function(x) {
  capture.output(vuong(m1 = genotype_poisson[[x]], m2 = genotype_ZIP[[x]]))
})
 

## Zero-inflated negative binomial regression:
# Make model
genotype_ZINB <- mclapply(seq_along(pop1_matList[[1]]), function(x) {
  zeroinfl(formula = c(pop1_matList[[1]][[x]], pop2_matList[[1]][[x]]) ~ genotype | 1,
           dist = "negbin")
}, mc.cores = detectCores())
print("Representative example: ZINB model for first genomic window")
print(summary(genotype_ZINB[[1]]))
# Get estimates
genotype_ZINB_estimates <- lapply(seq_along(genotype_ZINB), function(x) {
  coef(genotype_ZINB[[x]])
})
# Get coefficients
genotype_ZINB_coef <- lapply(seq_along(genotype_ZINB), function(x) {
    coef(summary(genotype_ZINB[[x]]))
})
# Get P-values
genotype_ZINB_Pvals <- sapply(seq_along(genotype_ZINB), function(x) {
  coef(summary(genotype_ZINB[[x]]))$count[11]
})
# Correct for multiple testing
genotype_ZINB_AdjPvals <- p.adjust(p = genotype_ZINB_Pvals, method = "BH")
# Get Bayesian Information Criterion (BIC); the smaller, the better the fit
genotype_ZINB_BIC <- sapply(seq_along(genotype_ZINB), function(x) {
  if( !is.infinite(AIC(genotype_ZINB[[x]], k = log(length(pop1_matList[[1]][[x]])))) ) {
    AIC(genotype_ZINB[[x]], k = log(length(pop1_matList[[1]][[x]])))
  } else {
    NA
  }
})
# Get Akaike Information Criterion (AIC); the smaller, the better the fit
genotype_ZINB_AIC <- sapply(seq_along(genotype_ZINB), function(x) {
  if( !is.infinite(AIC(genotype_ZINB[[x]], k = 2)) ) {
    AIC(genotype_ZINB[[x]], k = 2)
  } else {
    NA
  }
})
# Get 95% confidence intervals
genotype_ZINB_CIs <- lapply(seq_along(genotype_ZINB), function(x) {
  confint(genotype_ZINB[[x]])
})

# Use the model to predict mean crossover counts for each genotype
# and corresponding standard errors
# The model predicts the correct mean counts
genotype_ZINB_predict <- lapply(seq_along(genotype_ZINB), function(x) {
  cbind(data.frame(genotype = c(pop1Name, pop2Name)),
        mean = predict(genotype_ZINB[[x]],
                       newdata = data.frame(genotype = c(pop1Name, pop2Name)),
                       type = "response")
         ## Doesn't work for this type of model for some reason
#        SE = predict(genotype_ZINB[[x]],
#                     newdata = data.frame(genotype = c(pop1Name, pop2Name)),
#                     type = "response",
#                     se.fit = TRUE)$se.fit
       )
})
genotype_ZINB_correct_predictions <- sapply(seq_along(genotype_ZINB_predict), function(x) {
  c(print(all.equal(genotype_ZINB_predict[x][[1]]$mean[1],
                    pop1_means[[1]][x])),
    print(all.equal(genotype_ZINB_predict[x][[1]]$mean[2],
                    pop2_means[[1]][x])))
})
#print(paste0(as.character(sum(genotype_ZINB_correct_predictions)),
#             "/", as.character(length(genotype_ZINB_correct_predictions))))


# Evaluate whether ZINB model solves the problem of excess zeroes by
# predicting π and μ, and calculate the combined probability of no SNPs
# Use the predict() options "zero" and "count" to obtain π and μ
genotype_ZINB_zero <- lapply(seq_along(genotype_ZINB), function(x) {
  predict(genotype_ZINB[[x]], type = "zero") #  π
})
genotype_ZINB_count <- lapply(seq_along(genotype_ZINB), function(x) {
  predict(genotype_ZINB[[x]], type = "count") # μ
})
genotype_ZINB_zero_count_means <- sapply(seq_along(genotype_ZINB), function(x) {
  mean( genotype_ZINB_zero[[x]] + 
         ( (1-genotype_ZINB_zero[[x]]) * exp(-genotype_ZINB_count[[x]]) )
  )
})
print(paste0("genotype_ZINB models predict fewer than observed occurrences of zero SNPs in a window for ",
             as.character(sum(genotype_ZINB_zero_count_means < genotype_zobs)),
             "/", length(genotype_zobs), " windows"))
# ZINB model more accurately estimates the probability of zero SNPs in each window
# than ordinary poisson model

# Compare ZINB with ZIP model
# using the Vuong test
genotype_ZINB_vuong <- lapply(seq_along(genotype_ZINB), function(x) {
  capture.output(vuong(m1 = genotype_ZIP[[x]], m2 = genotype_ZINB[[x]]))
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
print(paste0("genotype_poisson BIC is < Normal BIC for ",
             sum(genotype_poisson_BIC_vs_genotype_normal_BIC, na.rm = T), "/", length(genotype_poisson), " windows"))
print(paste0("Normal BIC is NA for ",
             sum(genotype_normal_BIC_isNA), "/", length(genotype_normal), " scaled windows"))

# Determine if Poisson BIC is lower than ZIP BIC (indicating better fit)
genotype_poisson_BIC_vs_genotype_ZIP_BIC <- sapply(seq_along(genotype_poisson), function(x) {
  genotype_poisson_BIC[x] < genotype_ZIP_BIC[x]
})
print(paste0("genotype_poisson BIC is < genotype_ZIP BIC for ",
             sum(genotype_poisson_BIC_vs_genotype_ZIP_BIC, na.rm = T), "/", length(genotype_poisson), " windows"))

# Determine if ZINB BIC is lower than ZIP BIC (indicating better fit)
genotype_ZINB_BIC_vs_genotype_ZIP_BIC <- sapply(seq_along(genotype_ZINB), function(x) {
  genotype_ZINB_BIC[x] < genotype_ZIP_BIC[x]
})
print(paste0("genotype_ZINB BIC is < genotype_ZIP BIC for ",
             sum(genotype_ZINB_BIC_vs_genotype_ZIP_BIC, na.rm = T), "/", length(genotype_ZINB), " windows"))


# Note: negative binomial models should not be fit to these data
# due to absence of overdispersion 
## N
#negbin <- lapply(seq_along(pop1_matList[[1]]), function(x) {
#  glm.nb(formula = pops_winIndCOs_list[[x]] ~ genotype,
#         control = glm.control(maxit = 2000))
#})

# Plot genotype_poisson_BIC vs genotype_normal_BIC
pdf(paste0(outDir,
           "genotype_poisson_vs_normal_BICs_for_SNPs_around_",
           pop1Name, "_and_", pop2Name, "_crossovers_",
           flankName, "_flank_", winName, "_win.pdf"),
    height = 5, width = 8)
par(mfrow = c(1, 1))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(x = 1:length(genotype_poisson_BIC), y = genotype_poisson_BIC, col = "red", type = "p", pch = 19,
     xlab = paste0("Windows (", as.character(winSize), " bp)"),
     ylab = "BIC",
     main = paste0("genotype regression BICs for SNPs around \n",
                   pop1Name, " and ", pop2Name, " crossovers"),
     ylim = c(min(genotype_poisson_BIC, genotype_normal_BIC),
              max(genotype_poisson_BIC, genotype_normal_BIC)))
lines(x = 1:length(genotype_poisson_BIC), y = genotype_normal_BIC, col = "blue", type = "p", pch = 19)
legend("bottomleft",
       legend = c("Poisson regression", "Normal linear regression"),
       col = c("red", "blue"),
       text.col = c("red", "blue"),
       text.font = c(1, 1),
       ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
dev.off()

# Plot genotype_poisson_AIC vs genotype_normal_AIC
pdf(paste0(outDir,
           "genotype_poisson_vs_normal_AICs_for_SNPs_around_",
           pop1Name, "_and_", pop2Name, "_crossovers_",
           flankName, "_flank_", winName, "_win.pdf"),
    height = 5, width = 8)
par(mfrow = c(1, 1))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(x = 1:length(genotype_poisson_AIC), y = genotype_poisson_AIC, col = "red", type = "p", pch = 19,
     xlab = paste0("Windows (", as.character(winSize), " bp)"),
     ylab = "AIC",
     main = paste0("genotype regression AICs for SNPs around \n",
                   pop1Name, " and ", pop2Name, " crossovers"),
     ylim = c(min(genotype_poisson_AIC, genotype_normal_AIC),
              max(genotype_poisson_AIC, genotype_normal_AIC)))
lines(x = 1:length(genotype_poisson_AIC), y = genotype_normal_AIC, col = "blue", type = "p", pch = 19)
legend("bottomleft",
       legend = c("Poisson regression", "Normal linear regression"),
       col = c("red", "blue"),
       text.col = c("red", "blue"),
       text.font = c(1, 1),
       ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
dev.off()

# Plot genotype_poisson_BIC vs genotype_ZIP_BIC
pdf(paste0(outDir,
           "genotype_poisson_vs_ZIP_BICs_for_SNPs_around_",
           pop1Name, "_and_", pop2Name, "_crossovers_",
           flankName, "_flank_", winName, "_win.pdf"),
    height = 5, width = 8)
par(mfrow = c(1, 1))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(x = 1:length(genotype_poisson_BIC), y = genotype_poisson_BIC, col = "red", type = "p", pch = 19,
     xlab = paste0("Windows (", as.character(winSize), " bp)"),
     ylab = "BIC",
     main = paste0("genotype regression BICs for SNPs around \n",
                   pop1Name, " and ", pop2Name, " crossovers"),
     ylim = c(min(genotype_poisson_BIC, genotype_ZIP_BIC),
              max(genotype_poisson_BIC, genotype_ZIP_BIC)))
lines(x = 1:length(genotype_poisson_BIC), y = genotype_ZIP_BIC, col = "blue", type = "p", pch = 19)
legend("bottomleft",
       legend = c("Poisson regression", "ZIP regression"),
       col = c("red", "blue"),
       text.col = c("red", "blue"),
       text.font = c(1, 1),
       ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
dev.off()

# Plot genotype_poisson_AIC vs genotype_ZIP_AIC
pdf(paste0(outDir,
           "genotype_poisson_vs_ZIP_AICs_for_SNPs_around_",
           pop1Name, "_and_", pop2Name, "_crossovers_",
           flankName, "_flank_", winName, "_win.pdf"),
    height = 5, width = 8)
par(mfrow = c(1, 1))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(x = 1:length(genotype_poisson_AIC), y = genotype_poisson_AIC, col = "red", type = "p", pch = 19,
     xlab = paste0("Windows (", as.character(winSize), " bp)"),
     ylab = "AIC",
     main = paste0("genotype regression AICs for SNPs around \n",
                   pop1Name, " and ", pop2Name, " crossovers"),
     ylim = c(min(genotype_poisson_AIC, genotype_ZIP_AIC),
              max(genotype_poisson_AIC, genotype_ZIP_AIC)))
lines(x = 1:length(genotype_poisson_AIC), y = genotype_ZIP_AIC, col = "blue", type = "p", pch = 19)
legend("bottomleft",
       legend = c("Poisson regression", "ZIP regression"),
       col = c("red", "blue"),
       text.col = c("red", "blue"),
       text.font = c(1, 1),
       ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
dev.off()

# Plot genotype_ZINB_BIC vs genotype_ZIP_BIC
pdf(paste0(outDir,
           "genotype_ZINB_vs_ZIP_BICs_for_SNPs_around_",
           pop1Name, "_and_", pop2Name, "_crossovers_",
           flankName, "_flank_", winName, "_win.pdf"),
    height = 5, width = 8)
par(mfrow = c(1, 1))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(x = 1:length(genotype_ZINB_BIC), y = genotype_ZINB_BIC, col = "red", type = "p", pch = 19,
     xlab = paste0("Windows (", as.character(winSize), " bp)"),
     ylab = "BIC",
     main = paste0("genotype regression BICs for SNPs around \n",
                   pop1Name, " and ", pop2Name, " crossovers"),
     ylim = c(min(genotype_ZINB_BIC, genotype_ZIP_BIC),
              max(genotype_ZINB_BIC, genotype_ZIP_BIC)))
lines(x = 1:length(genotype_ZINB_BIC), y = genotype_ZIP_BIC, col = "blue", type = "p", pch = 19)
legend("bottomleft",
       legend = c("ZINB regression", "ZIP regression"),
       col = c("red", "blue"),
       text.col = c("red", "blue"),
       text.font = c(1, 1),
       ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
dev.off()

# Plot genotype_ZINB_AIC vs genotype_ZIP_AIC
pdf(paste0(outDir,
           "genotype_ZINB_vs_ZIP_AICs_for_SNPs_around_",
           pop1Name, "_and_", pop2Name, "_crossovers_",
           flankName, "_flank_", winName, "_win.pdf"),
    height = 5, width = 8)
par(mfrow = c(1, 1))
par(mar = c(4.1, 4.1, 3.1, 4.1))
par(mgp = c(3, 1, 0))
plot(x = 1:length(genotype_ZINB_AIC), y = genotype_ZINB_AIC, col = "red", type = "p", pch = 19,
     xlab = paste0("Windows (", as.character(winSize), " bp)"),
     ylab = "AIC",
     main = paste0("genotype regression AICs for SNPs around \n",
                   pop1Name, " and ", pop2Name, " crossovers"),
     ylim = c(min(genotype_ZINB_AIC, genotype_ZIP_AIC),
              max(genotype_ZINB_AIC, genotype_ZIP_AIC)))
lines(x = 1:length(genotype_ZINB_AIC), y = genotype_ZIP_AIC, col = "blue", type = "p", pch = 19)
legend("bottomleft",
       legend = c("ZINB regression", "ZIP regression"),
       col = c("red", "blue"),
       text.col = c("red", "blue"),
       text.font = c(1, 1),
       ncol = 1, cex = 0.7, lwd = 1.5, bty = "n")
dev.off()


# Function for plotting average SNP frequency profiles
plotAvgSNPfreq <- function(dat1, dat2,
                           ran1, ran2,
                           col1, col2,
                           mainTitle,
                           flankSize, winSize,
                           flankLabL1, flankLabR1,
                           flankLabL2, flankLabR2,
                           legendLoc, legendLabs) {
  plot(x = 1:(((flankSize*2)+winSize)/winSize),
       y = dat1, col = col1,
       type = "l", lwd = 3, ann = F,
       ylim = c(min(dat1, dat2, ran1, ran2),
                max(dat1, dat2, ran1, ran2)),
       xaxt = "n", yaxt = "n")
  lines(x = 1:(((flankSize*2)+winSize)/winSize),
        y = dat2, col = col2,
        type = "l", lwd = 3)
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle)
  axis(side = 2, at = pretty(c(dat1, dat2, ran1, ran2)), lwd = 2, cex.axis = 0.7)
  mtext(side = 2, line = 2.0, cex = 0.8,
        text = bquote("SNP frequency" ~
                      .(as.character(winSize)) ~ "bp"^-1))
  axis(side = 1, lwd = 2,
       at = c(1,
              ((flankSize/winSize)/2),
              ((flankSize+winSize)/winSize),
              (((flankSize*2)+winSize)/winSize)-((flankSize/winSize)/2),
              (((flankSize*2)+winSize)/winSize)),
       labels = c("", "", "", "", ""))
  mtext(side = 1, line = 0.5, cex = 0.7,
        at = c(1,
               ((flankSize/winSize)/2),
               ((flankSize+winSize)/winSize),
               (((flankSize*2)+winSize)/winSize)-((flankSize/winSize)/2),
               (((flankSize*2)+winSize)/winSize)),
        text = c(flankLabL1, flankLabL2, "Midpoint", flankLabR2, flankLabR1))
  abline(v = (flankSize+winSize)/winSize, lty = 3, lwd = 2)
  box(lwd = 2)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col1, col2),
         text.col = c(col1, col2),
         text.font = c(1, 1),
         ncol = 1, cex = 0.7, lwd = 2, bty = "n")
}

# Function to plot profile of -log10-transformed P-values and adjusted P-values
minusLog10PvalPlot <- function(Pvals, AdjPvals,
                               PvalsCol, AdjPvalsCol,
                               mainTitle,
                               flankSize, winSize,
                               flankLabL1, flankLabR1,
                               flankLabL2, flankLabR2) {
  plot(x = 1:(((flankSize*2)+winSize)/winSize),
       y = -log10(Pvals), col = PvalsCol,
       type = "l", lwd = 3, ann = F,
       ylim = c(0,
                pmax(-log10(0.05), max(-log10(Pvals), na.rm = T))),
       xaxt = "n", yaxt = "n")
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle)
  mtext(side = 2, line = 2.0, cex = 0.8, col = PvalsCol,
        text = bquote("-Log"[10]*"("*italic("P")*"-value)"))
  axis(side = 2, cex.axis = 0.7, lwd = 2)
  abline(h = -log10(0.05), lty = 5, lwd = 2, col = PvalsCol)
  par(new = T)
  plot(x = 1:(((flankSize*2)+winSize)/winSize),
       y = -log10(AdjPvals), col = AdjPvalsCol,
       type = "l", lwd = 3,
       ylim = c(0,
                pmax(-log10(FDR), max(-log10(AdjPvals), na.rm = T))),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = "",
       cex.main = 1)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1, adj = c(0.5, -2.5), xpd = NA, srt = -90, col = AdjPvalsCol,
       labels = bquote("-Log"[10]*"(BH-adjusted "*italic("P")*"-value)"))
  axis(side = 4, cex.axis = 0.7, lwd = 2)
  abline(h = -log10(0.1), lty = 5, lwd = 2, col = AdjPvalsCol)
  axis(side = 1, lwd = 2,
       at = c(1,
              ((flankSize/winSize)/2),
              ((flankSize+winSize)/winSize),
              (((flankSize*2)+winSize)/winSize)-((flankSize/winSize)/2),
              (((flankSize*2)+winSize)/winSize)),
       labels = c("", "", "", "", ""))
  mtext(side = 1, line = 0.5, cex = 0.7,
        at = c(1,
               ((flankSize/winSize)/2),
               ((flankSize+winSize)/winSize),
               (((flankSize*2)+winSize)/winSize)-((flankSize/winSize)/2),
               (((flankSize*2)+winSize)/winSize)),
        text = c(flankLabL1, flankLabL2, "Midpoint", flankLabR2, flankLabR1))
  abline(v = (flankSize+winSize)/winSize, lty = 3, lwd = 2)
  box(lwd = 2)
}

# Function to plot profile of -log10-transformed adjusted P-values and estimates
minusLog10AdjPvalEstimatePlot <- function(AdjPvals, ests,
                                          AdjPvalsCol, estsCol,
                                          mainTitle,
                                          flankSize, winSize,
                                          flankLabL1, flankLabR1,
                                          flankLabL2, flankLabR2) {
  plot(x = 1:(((flankSize*2)+winSize)/winSize),
       y = -log10(AdjPvals), col = AdjPvalsCol,
       type = "l", lwd = 3, ann = F,
       ylim = c(0,
                pmax(-log10(0.05), max(-log10(AdjPvals), na.rm = T))),
       xaxt = "n", yaxt = "n")
  mtext(side = 3, line = 0.5, cex = 0.8, text = mainTitle)
  mtext(side = 2, line = 2.0, cex = 0.8, col = AdjPvalsCol,
        text = bquote("-Log"[10]*"(BH-adjusted "*italic("P")*"-value)"))
  axis(side = 2, cex.axis = 0.7, lwd = 2)
  abline(h = -log10(0.1), lty = 5, lwd = 2, col = AdjPvalsCol)
  par(new = T)
  plot(x = 1:(((flankSize*2)+winSize)/winSize),
       y = ests, col = estsCol,
       type = "l", lwd = 3,
       ylim = c(min(ests), max(ests)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = "",
       cex.main = 1)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 1, adj = c(0.5, -3), xpd = NA, srt = -90, col = estsCol,
       labels = bquote("Effect size estimate"))
  axis(side = 4, cex.axis = 0.7, lwd = 2)
  axis(side = 1, lwd = 2,
       at = c(1,
              ((flankSize/winSize)/2),
              ((flankSize+winSize)/winSize),
              (((flankSize*2)+winSize)/winSize)-((flankSize/winSize)/2),
              (((flankSize*2)+winSize)/winSize)),
       labels = c("", "", "", "", ""))
  mtext(side = 1, line = 0.5, cex = 0.7,
        at = c(1,
               ((flankSize/winSize)/2),
               ((flankSize+winSize)/winSize),
               (((flankSize*2)+winSize)/winSize)-((flankSize/winSize)/2),
               (((flankSize*2)+winSize)/winSize)),
        text = c(flankLabL1, flankLabL2, "Midpoint", flankLabR2, flankLabR1))
  abline(v = (flankSize+winSize)/winSize, lty = 3, lwd = 2)
  box(lwd = 2)
}

# Get estimates
genotype_ZINB_estimates_vector <- sapply(seq_along(genotype_ZINB_estimates), function(x) {
  genotype_ZINB_estimates[[x]][[2]]
})

# Plot
pdf(paste0(plotDir, "SNP_frequency_around_", pop1Name, "_v_", pop2Name,
           "_crossovers_and_ranLoc_",
           flankName, "_flank_", winName, "_win.pdf"),
    height = 8, width = 9)
par(mfcol = c(2, 2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 0.75, 0))
plotAvgSNPfreq(dat1 = colMeans(pop1_matList[[1]]),
               dat2 = colMeans(pop1_matList[[2]]),
               ran1 = colMeans(pop2_matList[[1]]),
               ran2 = colMeans(pop2_matList[[2]]),
               col1 = "forestgreen",
               col2 = "grey50",
               mainTitle = bquote("SNP frequency around"~.(pop1Name)~
                                  "crossovers and random loci"),
               flankSize = flankSize,
               winSize = winSize,
               flankLabL1 = paste0("-", as.character(flankSize/1000), " kb"),
               flankLabR1 = paste0(as.character(flankSize/1000), " kb"),
               flankLabL2 = paste0("-", as.character((flankSize/1000)/2), " kb"),
               flankLabR2 = paste0(as.character((flankSize/1000)/2), " kb"),
               legendLoc = "topleft",
               legendLabs = c("Crossovers", "Random"))
plotAvgSNPfreq(dat1 = colMeans(pop2_matList[[1]]),
               dat2 = colMeans(pop2_matList[[2]]),
               ran1 = colMeans(pop1_matList[[1]]),
               ran2 = colMeans(pop1_matList[[2]]),
               col1 = "magenta3",
               col2 = "grey30",
               mainTitle = bquote("SNP frequency around"~.(pop2Name)~
                                  "crossovers and random loci"),
               flankSize = flankSize,
               winSize = winSize,
               flankLabL1 = paste0("-", as.character(flankSize/1000), " kb"),
               flankLabR1 = paste0(as.character(flankSize/1000), " kb"),
               flankLabL2 = paste0("-", as.character((flankSize/1000)/2), " kb"),
               flankLabR2 = paste0(as.character((flankSize/1000)/2), " kb"),
               legendLoc = "topleft",
               legendLabs = c("Crossovers", "Random"))
plotAvgSNPfreq(dat1 = colMeans(pop1_matList[[1]]),
               dat2 = colMeans(pop2_matList[[1]]),
               ran1 = colMeans(pop1_matList[[2]]),
               ran2 = colMeans(pop2_matList[[2]]),
               col1 = "forestgreen",
               col2 = "magenta3",
               mainTitle = bquote("SNP frequency around"~.(pop1Name)~
                                  "vs"~.(pop2Name)*" crossovers"),
               flankSize = flankSize,
               winSize = winSize,
               flankLabL1 = paste0("-", as.character(flankSize/1000), " kb"),
               flankLabR1 = paste0(as.character(flankSize/1000), " kb"),
               flankLabL2 = paste0("-", as.character((flankSize/1000)/2), " kb"),
               flankLabR2 = paste0(as.character((flankSize/1000)/2), " kb"),
               legendLoc = "topleft",
               legendLabs = c(pop1Name, pop2Name))
minusLog10AdjPvalEstimatePlot(AdjPvals = genotype_ZINB_AdjPvals,
                              ests = genotype_ZINB_estimates_vector,
                              AdjPvalsCol = "red",
                              estsCol = "dodgerblue3",
                              mainTitle = bquote("SNP frequency around"~.(pop1Name)~
                                                 "vs"~.(pop2Name)*" crossovers"),
                              flankSize = flankSize,
                              winSize = winSize,
                              flankLabL1 = paste0("-", as.character(flankSize/1000), " kb"),
                              flankLabR1 = paste0(as.character(flankSize/1000), " kb"),
                              flankLabL2 = paste0("-", as.character((flankSize/1000)/2), " kb"),
                              flankLabR2 = paste0(as.character((flankSize/1000)/2), " kb"))
#minusLog10PvalPlot(Pvals = genotype_ZINB_Pvals,
#                   AdjPvals = genotype_ZINB_AdjPvals,
#                   PvalsCol = "dodgerblue3",
#                   AdjPvalsCol = "red",
#                   mainTitle = bquote("SNP frequency around"~.(pop1Name)~
#                                  "vs"~.(pop2Name)*" crossovers"),
#                   flankSize = flankSize,
#                   winSize = winSize,
#                   flankLabL1 = paste0("-", as.character(flankSize/1000), " kb"),
#                   flankLabR1 = paste0(as.character(flankSize/1000), " kb"),
#                   flankLabL2 = paste0("-", as.character((flankSize/1000)/2), " kb"),
#                   flankLabR2 = paste0(as.character((flankSize/1000)/2), " kb"))
dev.off()

