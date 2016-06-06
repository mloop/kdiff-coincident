# Purpose: determine whether the snapped datasets satisfy sanity checks

# Author: Matthew Shane Loop

# Preliminaries

## Packages
library(dplyr)
library(spatstat)

## Load conditions
conditions <- read.table(file = 'conditions.txt', sep = '\t', header = TRUE)

## Set parameters
B <- 1000  ## Setting the full 6,200 number of iterations makes computation too long

# Determine whether proportion of cases and controls being snapped is equal
snapped_datasets <- data.frame()
for(j in 3:nrow(conditions)){  ## First 2 conditions have no snapped points (p = 0)
 for(i in 1:B){
     snapped <- read.table(file = paste("../data/datasets/c-", conditions$condition[j], "/", i, "/snapped.txt", sep = ""), sep = "\t", header = T)
     snapped_pp <- ppp(x = snapped$x, y = snapped$y, marks = as.factor(snapped$marks))
     duplicates <- snapped_pp[(duplicated(snapped_pp, rule = 'unmark') == TRUE)]
     proportion_cases <- duplicates[duplicates$marks == 'case']$n/duplicates$n
     data <- data.frame(condition = conditions$condition[j], hypothesis = conditions$hypothesis[j], mismeasurement_probability = conditions$mismeasurement_probabilities[j], iteration = i, proportion_cases = proportion_cases)
     snapped_datasets <- rbind(snapped_datasets, data)
 }
}
snapped_datasets <- tbl_df(snapped_datasets)

summarize(snapped_datasets,
          proportion_cases_snapped = mean(proportion_cases))

## Conclusion: cases and controls are being snapped in approximately equal proportions

# Determine whether number of points snapped matches specified snapping probability
snapped_datasets <- data.frame()
for(j in 3:nrow(conditions)){  ## First 2 conditions have no snapped points (p = 0)
    for(i in 1:B){
        snapped <- read.table(file = paste("../data/datasets/c-", conditions$condition[j], "/", i, "/snapped.txt", sep = ""), sep = "\t", header = T)
        snapped_pp <- ppp(x = snapped$x, y = snapped$y, marks = as.factor(snapped$marks))
        duplicates <- snapped_pp[(duplicated(snapped_pp, rule = 'unmark') == TRUE)]
        proportion_snapped <- duplicates$n/snapped_pp$n
        data <- data.frame(condition = conditions$condition[j], hypothesis = conditions$hypothesis[j], mismeasurement_probability = conditions$mismeasurement_probabilities[j], iteration = i, proportion_snapped = proportion_snapped)
        snapped_datasets <- rbind(snapped_datasets, data)
    }
}
snapped_datasets <- tbl_df(snapped_datasets)
snapped_grouped <- group_by(snapped_datasets, mismeasurement_probability)

summarize(snapped_grouped,
          estimated_snap_probability = mean(proportion_snapped))

## Conclusion: points are bein
