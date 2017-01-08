rm(list=ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

main_dir <- file.path(CORE_DIR, "work/mQTL_analysis_example/")
setwd(main_dir)

## -----------------------------------------------------------------------------
## the following packages need to be installed: 
## parallel, varbvs
## -----------------------------------------------------------------------------

Rdata_obj_dir <- "path/to/simulated/data"
load(file.path(Rdata_obj_dir, "simulated_data_from_real_SNPs.RData"))

results_dir <- paste(CORE_DIR, "results/Repl_mQTL_analysis_seed_", my_seed, "/", sep = "")
load(file.path(results_dir, "hyper_varbvs.RData"))

my_seed <- 111
RNGkind("L'Ecuyer-CMRG") # to ensure reproducibility when using parallel processes
set.seed(my_seed)

n_perm <- 25

n_cpus <- 12
verbose <- F

bool_save <- T

if( bool_save ){
  results_dir <- paste(CORE_DIR, "Permute_repl_varbvs_seed_", vec_seed, "/", sep = "")
  dir.create(results_dir)
  
  sink(paste(results_dir, "out.txt", sep=""), append = F, split = T,
       type = "output")
  sink(file(paste(results_dir, "err.txt", sep=""), open = "wt"), type = "message")
}

require(varbvs)
require(parallel)

X <- apply(dat$snps, 2, as.double)
Y <- dat$phenos
rownames(Y) <- NULL

for (i in 1:n_perm) {
  
  ind_perm <- sample(1:n)
  
  ppi <- mclapply(1:d, function(k) {
    out <- varbvs(X, Z = NULL, Y[ind_perm, k], family ="gaussian",
                  sigma = list_hyper$sigma[,k], sa = list_hyper$sa[,k],
                  logodds = list_hyper$logodds[,k], update.sigma = F,
                  update.sa = F, verbose = verbose)
    w <- normalizelogweights(out$logw)
    out$alpha %*% c(w)
  }, mc.cores=n_cpus)
  ppi <- do.call(cbind, ppi)
  
  rownames(ppi) <- colnames(X)
  colnames(ppi) <- colnames(dat$phenos)
  
  save(ind_perm, ppi,
       file = paste(results_dir, "varbvs_real_data_", i, ".RData", sep=""))
}


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
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods
# [8] base
#
# other attached packages:
#   [1] varbvs_2.0.0
#
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-2  latticeExtra_0.6-28 grid_3.2.0
# [4] lattice_0.20-31
