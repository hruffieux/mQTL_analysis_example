rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

main_dir <- file.path(CORE_DIR, "work/mQTL_analysis_example/")
setwd(main_dir)

## -----------------------------------------------------------------------------
## Joint inference on large mQTL data
## ! Recommended RAM 512G, if not enough memory is available the code will crash !
##
## The following packages need to be installed:
## gplots, ROCR, locus.
## The version of locus used in this repository can be installed using:
## devtools::install_github("hruffieux/locus", ref = "v0.1.1")
## -----------------------------------------------------------------------------


require(locus)

RNGkind("L'Ecuyer-CMRG") # to ensure reproducibility when using parallel processes
my_seed <- 123
set.seed(my_seed)



## ---------------------------------------------------------------------
## --- THIS CODE CHUNK WAS USED TO SIMULATE DATA FROM REAL SNP DATA ----
## --  it can't be executed as we do not provide the real SNPs but the -
## --- obtained simulated SNPs can be downloaded and loaded below ------
#

## simulates SNPs by emulating the correlation structure and minor allele
## frequencies of real tag SNPs
#
snp_path <- file.path(CORE_DIR, "data/mQTL_baseline_data.RData")
load(snp_path)

n <- 350

p <- NULL # total number of SNPs available
bl_lgth <- min(1000, p) # block width (for dependence structure replication)

n_cpus <- 12

list_snps <- replicate_real_snps(n = n, real_snps = SNPs, bl_lgth = bl_lgth,
                                 p = p, n_cpus = n_cpus)
p <- ncol(list_snps$snps)


## simulates molecular outcomes with correlation structure by blocks, similarly
## to the real metabolite expression levels
#
d <- 250 # number of outcomes

# 7 blocks of equicorrelated outcomes
cor_type_ph <- "equicorrelated"
vec_rho_ph <- sample(seq(0.5, 0.8, by = 0.1), 7, replace = TRUE)
var_err <- 1

list_phenos <- generate_phenos(n = n, d = d, var_err = var_err,
                               cor_type = cor_type_ph, vec_rho = vec_rho_ph,
                               n_cpus = n_cpus)


## simulates assocation pattern between the SNPs and the outcomes
#
p0 <- 750  # number of active SNPs
ind_p0 <- sample(1:p, p0, replace = FALSE)

d0 <- 175  # number of active outcomes
ind_d0 <- sample(1:d, d0, replace = FALSE)


# associations distributed across blocks of correlated outcomes so that a given
# active SNP
n_freq <- 200
sh2_freq <- 500

n_rare <- 5
sh2_rare <- 5

# for each active snp, vec_prob_sh will be sampled with replacement 7x
# (= # of correlated outcome blocks)
#
vec_prob_sh <- c(rbeta(n_freq, shape1 = 1, shape2 = sh2_freq),
                 rbeta(n_rare, shape1 = 1, shape2 = sh2_rare))

max_tot_pve <- 0.4 # maximum total proportion of outcome variance explained by the SNPs.

dat <- generate_dependence(list_snps, list_phenos, ind_d0, ind_p0, vec_prob_sh,
                           family = "gaussian", max_tot_pve = max_tot_pve)

Rdata_obj_dir <- file.path(CORE_DIR, "data/")
save(bl_lgth, cor_type_ph, d, d0, dat, ind_d0, ind_p0, list_phenos, list_snps,
     max_tot_pve, n, n_freq, n_rare, p, p0, sh2_freq, sh2_rare, vec_prob_sh,
     vec_rho_ph, var_err,
     file = file.path(Rdata_obj_dir, "simulated_data_from_real_SNPs.RData"))

## ------------------------------- END ---------------------------------



## ---------------------------------------------------------------------
## ---------------- LOAD HERE THE DATA SIMULATED ABOVE  ----------------
## -------- the .RData file can be downloaded from figshare ------------
## ------- https://dx.doi.org/10.6084/m9.figshare.4509755.v1 -----------


Rdata_obj_dir <- "path/to/simulated/data"
load(file.path(Rdata_obj_dir, "simulated_data_from_real_SNPs.RData"))

## ------------------------------- END ---------------------------------


## --------------------------- VB INFERENCE ----------------------------


tol <- 1e-6

cv <- TRUE

if( cv ){
  p0_av <- NULL
  n_folds <- 3
  size_p0_av_grid <- 5
  n_cpus <- 3
  tol_cv <- 1e-3
  list_cv <- set_cv(n, p, n_folds, size_p0_av_grid, n_cpus, tol_cv)
} else {
  p0_av <- min(1000, p / 3)
  list_cv <- NULL
}

bool_save <- TRUE

if( bool_save ){
  results_dir <- paste(CORE_DIR, "results/Repl_mQTL_analysis_seed_", my_seed,
                       "/", sep = "")
  dir.create(results_dir)

  sink(paste(results_dir, "out.txt", sep = ""), append = FALSE, split = TRUE,
       type = "output")
  sink(file(paste(results_dir, "err.txt", sep = ""), open = "wt"), type = "message")
}

if (bool_save) {
  save.image(file = paste(results_dir, "output_initial.RData", sep = ""))
  png(paste(results_dir,"hm_pat_ind_p0.png", sep = ""), w = 2500, h = 2750,
      res = 500, type = "cairo")
}
require(gplots)
hmcols<-colorRampPalette(c("grey98", "black"))(256)

par(cex.main = 0.8, cex.lab = 0.8)
heatmap.2(1*dat$pat[ind_p0,], dendrogram = "none", col = hmcols, Rowv = FALSE,
          Colv = FALSE, main = "Simulated pattern",
          density.info = "none", trace = "none", key = FALSE,
          xlab = "Y", ylab = "X (Only rows with at least one active covariate)",
          lhei = c(1, 3), lwid = c(0.2, 1),
          labRow = FALSE, labCol = FALSE, margins = c(3, 3))
if (bool_save) {
  dev.off()
  png(paste(results_dir,"hm_simulated_abs_beta_ind_p0.png", sep = ""),
      w = 2500, h = 2750, res = 500, type = "cairo")
}
par(cex.main = 0.8, cex.lab = 0.8)
heatmap.2(abs(dat$beta[ind_p0,]), dendrogram = "none", col = hmcols, Rowv = FALSE,
          Colv = FALSE, main = "Simulated |beta|",
          density.info = "none", trace = "none", key = FALSE,
          xlab = "Y", ylab = "X (Only rows with at least one active covariate)",
          lhei = c(1, 3), lwid = c(0.2, 1),
          labRow = FALSE, labCol = FALSE, margins = c(3, 3))
if (bool_save) {
  dev.off()
}


rt_vb <- system.time(out_vb <- locus(dat$phenos, dat$snps, p0_av, 
                                     list_cv = list_cv, user_seed = my_seed, 
                                     tol = tol, save_hyper = TRUE, 
                                     save_init = TRUE))

cat("VB run done in [seconds]: \n")
print(rt_vb)

if (bool_save) {
  save.image(file = paste(results_dir, "output_vb.RData", sep = ""))
  list_hyper <- out_vb$list_hyper
  list_init <- out_vb$list_init
  p0_av <- out_vb$p_star
  save(list_hyper, list_init, p0_av,
       file = paste(results_dir, "hyper_init_vb.RData", sep = ""))
  png(paste(results_dir,"hm_vb_ind_p0.png", sep = ""), w = 2500, h = 2750,
      res = 500, type = "cairo")
}
par(cex.main = 0.8, cex.lab = 0.8)
heatmap.2(out_vb$gam_vb[ind_p0,], dendrogram = "none", col = hmcols, Rowv = FALSE,
          Colv = FALSE, main = "VB, PPI",
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

perf_vb <- compute_roc(out_vb$gam_vb, dat$pat)

if (bool_save) {
  png( paste(results_dir,"ROC_curves_vb.png", sep = ""), w = 2750, h = 2750,
       res = 500, type = "cairo")
}
plot(perf_vb, col = "darkblue", type = "l", lwd = 2,
     main = paste("ROC curve VB \n p = ", format(p, scientific = FALSE),
                  ", d = ", d, ", n = ", n, sep = ""))
abline(0,1)
if (bool_save) {
  dev.off()
  save.image(file = paste(results_dir, "output_vb.RData", sep = ""))
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
#   [1] locus_0.1.0   ROCR_1.0-7    gplots_2.17.0 Matrix_1.2-2
#
# loaded via a namespace (and not attached):
#   [1] KernSmooth_2.23-14 gdata_2.17.0       grid_3.2.0         caTools_1.17.1
# [5] bitops_1.0-6       gtools_3.5.0       lattice_0.20-31

