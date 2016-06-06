# Purpose: estimate difference in K functions at h = 0.15 across all simulated datasets

# Author: Matthew Shane Loop

# Preliminaries
# Load arguments and set variables
if(TRUE){
    args=(commandArgs(TRUE))
    if(length(args)==0){
        print('no args')
        i=1; full=TRUE
    }else{
        print('eval agrs')
        for(i in 1:length(args)){
            eval(parse(text=args[[i]]))
        }
}}

## Packages
library(spatstat)
library(plyr)

## User-defined functions
source(file = '../../r-scripts/kdiff.R')
source(file = '../../r-scripts/d_alternative.R')

## Parameters
nsim <- 199
nrank <- 5
form <- c('true', 'snapped', 'removed', 'jittered_0_05', 'jittered_0_1', 'jittered_0_15')
B <- 1000
n <- 100
q <- 0.05
v <- 10
sigma <- 0.1


## Simulation conditions
conditions <- read.table(file = '../data/conditions.txt', sep = '\t', header = TRUE)

for(j in 1:nrow(conditions)){

     u <- 322478 + i*j*200
     set.seed(u, kind = "Mersenne-Twister", normal.kind = "Inversion")

     # Import datasets and format for spatstat

     ## True
     true_df <- read.table(file = paste("../data/datasets/c-", conditions$condition[j], "/", i, "/true.txt", sep = ""), sep = "\t", header = T)
     true <- ppp(x = true_df$x, y = true_df$y, marks = as.factor(true_df$marks))

     ## Snapped
     snapped_df <- read.table(file = paste("../data/datasets/c-", conditions$condition[j], "/", i, "/snapped.txt", sep = ""), sep = "\t", header = T)
     snapped <- ppp(x = snapped_df$x, y = snapped_df$y, marks = as.factor(snapped_df$marks))

     ## Removed
     removed_df <- read.table(file = paste("../data/datasets/c-", conditions$condition[j], "/", i, "/removed.txt", sep = ""), sep = "\t", header = T)
     removed <- ppp(x = removed_df$x, y = removed_df$y, marks = as.factor(removed_df$marks))

     ## Jittered, 0.05
     jittered_0_05_df <- read.table(file = paste("../data/datasets/c-", conditions$condition[j], "/", i, "/jittered-0-05.txt", sep = ""), sep = "\t", header = T)
     jittered_0_05 <- ppp(x = jittered_0_05_df$x, y = jittered_0_05_df$y, marks = as.factor(jittered_0_05_df$marks))

     ## Jittered, 0.1
     jittered_0_1_df <- read.table(file = paste("../data/datasets/c-", conditions$condition[j], "/", i, "/jittered-0-1.txt", sep = ""), sep = "\t", header = T)
     jittered_0_1 <- ppp(x = jittered_0_1_df$x, y = jittered_0_1_df$y, marks = as.factor(jittered_0_1_df$marks))

     ## Jittered, 0.15
     jittered_0_15_df <- read.table(file = paste("../data/datasets/c-", conditions$condition[j], "/", i, "/jittered-0-15.txt", sep = ""), sep = "\t", header = T)
     jittered_0_15 <- ppp(x = jittered_0_15_df$x, y = jittered_0_15_df$y, marks = as.factor(jittered_0_15_df$marks))


     # Perform the random labeling nsim times determine the high and low ranks based upon crit
     r <- c(0, conditions$h[j])

     ## True
     true_kdiff <- envelope(true, kdiff, r = r, cr = "iso", nsim = nsim, nrank = nrank, verbose = FALSE, savefuns = TRUE, simulate = expression(rlabel(true)))

     ## Snapped
     snapped_kdiff <- envelope(snapped, kdiff, r = r, cr = "iso", nsim = nsim, nrank = nrank, verbose = FALSE, savefuns = TRUE, simulate = expression(rlabel(snapped)))

     ## Removed
     removed_kdiff <- envelope(removed, kdiff, r = r, cr = "iso", nsim = nsim, nrank = nrank, verbose = FALSE, savefuns = TRUE, simulate = expression(rlabel(removed)))

     ## Jittered, 0.05
     jittered_0_05_kdiff <- envelope(jittered_0_05, kdiff, r = r, cr = "iso", nsim = nsim, nrank = nrank, verbose = FALSE, savefuns = TRUE, simulate = expression(rlabel(jittered_0_05)))

     ## Jittered, 0.1
     jittered_0_1_kdiff <- envelope(jittered_0_1, kdiff, r = r, cr = "iso", nsim = nsim, nrank = nrank, verbose = FALSE, savefuns = TRUE, simulate = expression(rlabel(jittered_0_1)))

     ## Jittered, 0.15
     jittered_0_15_kdiff <- envelope(jittered_0_15, kdiff, r = r, cr = "iso", nsim = nsim, nrank = nrank, verbose = FALSE, savefuns = TRUE, simulate = expression(rlabel(jittered_0_15)))





     reject <- rep(0, times = ((length(r) - 1)*length(form)))
     condition <- rep(conditions$condition[j], times = ((length(r) - 1)*length(form)))
     mismeasurement_probability <- rep(conditions$mismeasurement_probabilities[j], times = ((length(r) - 1)*length(form)))
     jitter_radius <- rep(conditions$jitter_radius[j], times = ((length(r) - 1)*length(form)))
     iteration <- rep(i, times = ((length(r) - 1)*length(form)))
     true_d <- rep(d_alternative(h = conditions$h[j], n = n, q = q, v =
v,
sigma = sigma), times = ((length(r) - 1)*length(form)))
     d_hat <- c(true_kdiff[["obs"]][-1], snapped_kdiff[["obs"]][-1], removed_kdiff[["obs"]][-1], jittered_0_05_kdiff[["obs"]][-1], jittered_0_1_kdiff[["obs"]][-1], jittered_0_15_kdiff[["obs"]][-1])

     lo <- c(true_kdiff[["lo"]][-1], snapped_kdiff[["lo"]][-1], removed_kdiff[["lo"]][-1], jittered_0_05_kdiff[["lo"]][-1], jittered_0_1_kdiff[["lo"]][-1], jittered_0_15_kdiff[["lo"]][-1])

     hi <- c(true_kdiff[["hi"]][-1], snapped_kdiff[["hi"]][-1], removed_kdiff[["hi"]][-1], jittered_0_05_kdiff[["hi"]][-1], jittered_0_1_kdiff[["hi"]][-1], jittered_0_15_kdiff[["hi"]][-1])
     range <- rep(conditions$h[j], times = ((length(r) - 1)*length(form)))

     data <- data.frame(condition = condition, mismeasurement_probability = mismeasurement_probability, range = range, iteration = iteration, true_d = true_d, form = form, d_hat = d_hat, lo = lo, hi = hi, reject = reject)
     data$reject[which(data$d_hat < data$lo | data$d_hat > data$hi)] <- 1
     write.table(data, paste("results-", conditions$condition[j], "-", i, ".txt", sep = ""), sep = "\t", row.names = FALSE, append = FALSE, eol = '\n')

}

# Reproducibility
sessionInfo()
