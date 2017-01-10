rm(list=ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

main_dir <- file.path(CORE_DIR, "work/mQTL_analysis_example/")
setwd(main_dir)

## -----------------------------------------------------------------------------
## Joint inference on large mQTL data
## ! Recommended RAM 512G, if not enough memory is available the code will crash !
##
## The locus package need to be installed: 
## devtools::install_github("hruffieux/locus", ref = "v0.1.1") 
## -----------------------------------------------------------------------------

Rdata_obj_dir <- "path/to/simulated/data"
load(file.path(Rdata_obj_dir, "simulated_data_from_real_SNPs.RData"))

results_dir <- paste(CORE_DIR, "results/Repl_mQTL_analysis_seed_", my_seed, "/", sep = "")
load(file.path(results_dir, "hyper_init_vb.RData"))


require(locus)

my_seed <- 111

batch <- T
tol <- 1e-3
maxit <- 5000

verbose <- T

n_perm <- 25
n_cpus <- 1

list_blocks <- NULL

results_dir <- paste(CORE_DIR, "results/Permute_repl_seed_", my_seed, "/", sep = "")
dir.create(results_dir)

sink(paste(results_dir, "out.txt", sep=""), append = F, split = T,
     type = "output")
sink(file(paste(results_dir, "err.txt", sep=""), open = "wt"), type = "message")

generate_null(n_perm = n_perm, Y = dat$phenos, X = dat$snps, p0_av = p0_av, 
              Z = NULL, list_hyper = list_hyper, list_init = list_init, 
              list_blocks = list_blocks, user_seed = my_seed, tol = tol, 
              maxit = maxit, batch = batch, verbose = verbose,
              results_dir = results_dir, n_cpus = n_cpus)

## Reproducibility
#
sessionInfo()
#
# R version 3.2.0 (2015-04-16)
# Platform: x86_64-unknown-linux-gnu (64-bit)
# Running under: Red Hat Enterprise Linux ComputeNode 7.2 (Maipo)
#
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] locus_0.1.0
