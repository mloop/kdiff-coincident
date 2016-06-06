# Purpose: generate datasets for simulation

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

## User-defined functions
source(file = '../../r-scripts/cluster_genesis.R')
source(file = '../../r-scripts/create_alternative.R')
source(file = '../../r-scripts/snap.R')
source(file = '../../r-scripts/jitter_overlapping.R')


## Import conditions
conditions <- read.table(file = 'conditions.txt', sep = '\t', header = TRUE)

## Parameters for dataset generation
centroid_coords <- data.frame(x = 0.5, y = 0.5)
n <- 100
q <- 0.05
v <- 10
sigma <- 0.1
B <- 1000

# Simulate datasets
troublemakers <- data.frame()
for(j in 1:nrow(conditions)){
    dir.create(c(paste('datasets/c-', conditions$condition[j], sep = "")), recursive = TRUE)

        dir.create(c(paste('datasets/c-', conditions$condition[j], '/', i, sep = "")), recursive = TRUE)
        u <- 5345343 + i*j*200
        set.seed(u, kind = "Mersenne-Twister", normal.kind = "Inversion")

        ## Create true dataset
        {tryCatch({true <- create_alternative(n = n, q = q, v = v, sigma = sigma)}, error = function(e) {
            u <- u + 1
            set.seed(u, kind = "Mersenne-Twister", normal.kind = "Inversion")
            true <- create_alternative(n = n, q = q, v = v, sigma = sigma)
            t <- data.frame(condition = j, iteration = i, old_seed = (u - 1), new_seed = u)
            write.table(t, file = 'troublemakers.txt', append = TRUE, sep = '\t', row.names = FALSE, col.names = FALSE)
        }, finally = function(f){
            return(true)
        })}

        ## Snap dataset
        snapped <- snap(true, centroid_coordinates = centroid_coords, prob = conditions$mismeasurement_probabilities[j])

        ## Remove dataset
        removed <- snapped[duplicated(snapped, rule = 'unmark') == FALSE]

    ## Jitter datasets

    ### Radius 0.05
        jittered_0.05 <- jitter_overlapping(snapped, r = 0.05)

        ### Radius 0.1
        jittered_0.1 <- jitter_overlapping(snapped, r = 0.1)

        ### Radius 0.15
        jittered_0.15 <- jitter_overlapping(snapped, r = 0.15)

    ## Write out datasets
        write.table(data.frame(true), paste("datasets/c-", conditions$condition[j], "/", i, '/true.txt', sep = ""), sep = "\t", row.names = FALSE)
        write.table(data.frame(snapped), paste("datasets/c-", conditions$condition[j], "/", i, '/snapped.txt', sep = ""), sep = "\t", row.names = FALSE)
        write.table(data.frame(removed), paste("datasets/c-", conditions$condition[j], "/", i, '/removed.txt', sep = ""), sep = "\t", row.names = FALSE)
    write.table(data.frame(jittered_0.05), paste("datasets/c-", conditions$condition[j], "/", i, '/jittered-0-05.txt', sep = ""), sep = "\t", row.names = FALSE)
    write.table(data.frame(jittered_0.1), paste("datasets/c-", conditions$condition[j], "/", i, '/jittered-0-1.txt', sep = ""), sep = "\t", row.names = FALSE)
    write.table(data.frame(jittered_0.15), paste("datasets/c-", conditions$condition[j], "/", i, '/jittered-0-15.txt', sep = ""), sep = "\t", row.names = FALSE)


}

# Reproducibility
sessionInfo()
