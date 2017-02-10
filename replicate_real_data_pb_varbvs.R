rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

main_dir <- file.path(CORE_DIR, "work/mQTL_analysis_example/")
setwd(main_dir)

## -----------------------------------------------------------------------------
## the following packages need to be installed: 
## gplots, parallel, ROCR, varbvs
## -----------------------------------------------------------------------------


Rdata_obj_dir <- "path/to/simulated/data"
load(file.path(Rdata_obj_dir, "simulated_data_from_real_SNPs.RData"))

results_dir <- paste(CORE_DIR, "results/Repl_mQTL_analysis_seed_", my_seed, "/", sep = "")

my_seed <- 123
RNGkind("L'Ecuyer-CMRG") # to ensure reproducibility when using parallel processes
set.seed(my_seed)


require(varbvs)
require(parallel)

n_cpus <- 12

bool_save <- TRUE

X <- apply(dat$snps, 2, as.double)
Y <- dat$phenos
d <- ncol(Y)

names_out_varbvs <- c("ppi", "sigma", "sa", "logodds")


rt_varbvs <- system.time(out_varbvs <- mclapply(1:d, function(k) {
  out <- varbvs(X, Z = NULL, Y[, k], family = "gaussian", update.sigma = TRUE,
                update.sa = TRUE, verbose = FALSE)
  w <- normalizelogweights(out$logw)
  
  out <- list(out$alpha %*% c(w), out$sigma, out$sa, out$logodds)
  names(out) <- names_out_varbvs
  out
}, mc.cores = n_cpus))


list_hyper <- lapply(names_out_varbvs[-1],
                     function(key) do.call(cbind, lapply(out_varbvs, `[[`, key)))
names(list_hyper) <- names_out_varbvs[-1]

ppi_varbvs <- do.call(cbind, lapply(out_varbvs, `[[`, "ppi"))
colnames(ppi_varbvs) <- colnames(Y)
rm(out_varbvs)

cat("varbvs run done in [seconds]: \n")
print(rt_varbvs)


if (bool_save) {
  save(list_hyper, file = paste(results_dir, "hyper_varbvs.RData", sep = ""))
  save.image(file = paste(results_dir, "output_varbvs.RData", sep = ""))
  png(paste(results_dir,"hm_varbvs_ind_p0.png", sep = ""), w = 2500, h = 2750,
      res = 500, type = "cairo")
}
require(gplots)
hmcols<-colorRampPalette(c("grey98", "black"))(256)

par(cex.main = 0.8, cex.lab = 0.8)
heatmap.2(ppi_varbvs[ind_p0,], dendrogram = "none", col = hmcols, Rowv = FALSE,
          Colv = FALSE, main = "varbvs, PPI",
          density.info = "none", trace = "none", key = FALSE,
          xlab = "Y", ylab = "X (Only rows with at least one active covariate)",
          lhei = c(1, 3), lwid = c(0.2, 1),
          labRow = FALSE, labCol = FALSE, margins = c(3, 3))
if (bool_save) {
  dev.off()
}

compute_roc <- function(ppi, pat) {
  vec_rank <- as.numeric(as.factor(ppi))
  vec_pat <- as.numeric(pat)
  
  require(ROCR)
  pred <- prediction(vec_rank, vec_pat)
  performance(pred, measure = "tpr", x.measure = "fpr")
}

perf_varbvs <- compute_roc(ppi_varbvs, dat$pat)

if (bool_save) {
  png(paste(results_dir,"ROC_curve_varbvs.png", sep = ""), w = 2750, h = 2750,
      res = 500, type = "cairo")
}
plot(perf_varbvs, col = "pink", type = "l", lwd = 2,
     main = paste("ROC curve varbvs\n p = ", format(p, scientific = FALSE), ", d = ",
                  d, ", n = ", n, sep = ""))
abline(0,1)

if (bool_save) {
  dev.off()
  save.image(file = paste(results_dir, "output_varbvs.RData", sep = ""))
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
#   [1] ROCR_1.0-7    gplots_2.17.0 varbvs_2.0.0
#
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-2  KernSmooth_2.23-14  latticeExtra_0.6-28
# [4] gdata_2.17.0        grid_3.2.0          caTools_1.17.1
# [7] bitops_1.0-6        gtools_3.5.0        lattice_0.20-31

