rm(list=ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

main_dir <- file.path(CORE_DIR, "work/mQTL_analysis_example/")
setwd(main_dir)


## -----------------------------------------------------------------------------
## 2 x (200 + 1) posterior inclusion probability matrices each of size 
## ~200,000 x 250 need to be loaded
## recommended RAM 200G, if not enough memory is available the code will crash.
##
## the following packages need to be installed: 
## ggplot2, grid, gridExtra, parallel, xtable
## -----------------------------------------------------------------------------

source("utils_permutations.R")

bool_save <- TRUE

main_res_dir <- file.path(CORE_DIR, "results/Repl_mQTL_analysis_seed_123/")



######
# VB #
######


## loads the ppi matrices
## ----------------------
#
# VB runs on permuted datasets
#
vec_seed <- 1:9 * 111
vec_dir <- paste(CORE_DIR, "results/Permute_repl_seed_", vec_seed, "/", sep = "")

counts_perm_runs(vec_dir)

meth <- "vb"
list_ppi_vb <- load_perm_ppi(vec_dir, meth, n_cpus = 1)


# actual VB run
#
load(file.path(main_res_dir, "output_vb.RData"))


## estimates FDR for given ppi thresholds
## --------------------------------------
#
vec_thres_vb <- seq(0.001, 0.999, by = 0.001)
FDR_vb_med <- estimate_FDR(list_ppi_vb, out_vb$gam_vb, vec_thres_vb, type = "median", 
                           n_cpus = 1)
rm(list_ppi_vb, vec_thres_vb)

if(bool_save) save(FDR_vb_med, file = file.path(main_res_dir, "FDR_vb.RData"))



##########
# VARBVS #
##########


## loads the ppi matrices
## ----------------------
#
# varbvs runs on permuted datasets
#
vec_seed <- 1:9 * 111
vec_dir <- paste(CORE_DIR, "Permute_repl_varbvs_seed_", vec_seed, "/", sep = "")

counts_perm_runs(vec_dir)

meth <- "varbvs"
list_ppi_varbvs <- load_perm_ppi(vec_dir, meth, n_cpus = 1)


# actual varbvs run
#
load(file.path(main_res_dir, "output_varbvs.RData"))


## estimates FDR for given ppi thresholds (adaptive thresholds)
## --------------------------------------
#
vec_thres_varbvs <- seq(0.001, 0.999, by = 0.001)

FDR_varbvs_med <- estimate_FDR(list_ppi_varbvs, ppi_varbvs, vec_thres_varbvs,
                               type = "median", n_cpus = 1, bool_adapt_thres = TRUE)

rm(list_ppi_varbvs, vec_thres_varbvs)

if(bool_save) save(FDR_varbvs_med, file = file.path(main_res_dir, "FDR_varbvs.RData"))




##############################################################
# SPLINE FIT AND FDR-BASED COMPARISION BETWEEN VB AND VARBVS #
##############################################################

## computes thresholds for desired FDR values
## ------------------------------------------
#
vec_FDR_desired <- seq(0.2, 0.05, by=-0.05) # here we used 0.2 instead of 0.25
# (as FDR of 25 % produces PPI maps
# that are too large to be displayed)
# but 0.25 for the table below
# spline fit
title <- "vb"
hm_vb <- prepare_hm(FDR_vb_med$vec_FDR, names(FDR_vb_med$vec_FDR), vec_FDR_desired,
                    out_vb$gam_vb, title, main_res_dir, spar = 0.7, bool_save)


title <- "varbvs"
# normalized posterior inclusion probabilities since adaptive threshold for varbvs
ppi_varbvs_adapt_thres <- apply(ppi_varbvs, 2,
                                function(ppi_t) ppi_t / median(ppi_t)) * median(ppi_varbvs)

hm_varbvs <- prepare_hm(FDR_varbvs_med$vec_FDR, names(FDR_varbvs_med$vec_FDR),
                        vec_FDR_desired, ppi_varbvs_adapt_thres, title, out_dir,
                        spar = 0.7, bool_save)

# indices of snps (ind_s) and metabolites (ind_t) having significant associations
# either for VB or for varbvs
ind_s <- sort(unique(c(hm_vb$tab[,1], hm_varbvs$tab[,1])))
ind_t <- sort(unique(c(hm_vb$tab[,2], hm_varbvs$tab[,2])))

# maximum values for side panels
x_max <- max(c(rowSums(ppi_varbvs_adapt_thres[ind_s, ind_t] > hm_varbvs$vec_thres_select[1]),
               rowSums(out_vb$gam_vb[ind_s, ind_t] > hm_vb$vec_thres_select[1])))
y_max <- max(c(colSums(ppi_varbvs_adapt_thres[ind_s, ind_t] > hm_varbvs$vec_thres_select[1]),
               colSums(out_vb$gam_vb[ind_s, ind_t] > hm_vb$vec_thres_select[1])))


## comparison with the simulated pattern
## -------------------------------------
#
true_pat <- (!dat$pat)*1
thres <- 0
text <- "X = simulated \n pattern"
symb <- "X"
cexc <- 5
color <- "red"
obj_comp <- make_obj_compare(true_pat, thres, text, symb, color, cexc)


title <- "VB"
title_file <- "VB_FDR_map"

make_FDR_heatmap(out_vb$gam_vb, ind_s, ind_t, hm_vb$FDR_select, x_max, y_max,
                 title, out_dir, title_file, bool_save, bool_repl = TRUE, obj_comp)

title <- "varbvs"
title_file <- "varbvs_FDR_map"

make_FDR_heatmap(ppi_varbvs_adapt_thres, ind_s, ind_t, hm_varbvs$FDR_select,
                 x_max, y_max, title, out_dir, title_file, bool_save,
                 bool_repl = TRUE, obj_comp)


## number of correct associations recovered at prescribed FDR values
## -----------------------------------------------------------------
#
list_tab_tp_vb <- sapply(hm_vb$vec_thres_select,
                         function(thres) which(out_vb$gam_vb > thres & dat$pat,
                                               arr.ind = TRUE))
n_tp_vb <- unlist(lapply(list_tab_tp_vb, nrow))

list_tab_tp_varbvs <- sapply(hm_varbvs$vec_thres_select,
                             function(thres) which(ppi_varbvs_adapt_thres > thres & dat$pat,
                                                   arr.ind = TRUE))
n_tp_varbvs <- unlist(lapply(list_tab_tp_varbvs, nrow))

n_inter_tp <- sapply(1:length(vec_FDR_desired),
                     function(i) nrow(merge(list_tab_tp_vb[[i]], list_tab_tp_varbvs[[i]],
                                            by = c("row", "col"))))

tb <- apply(cbind(vec_FDR_desired * 100, n_tp_vb, n_tp_varbvs, n_inter_tp), 2, rev)
colnames(tb) <- c("Permutation-based FDR (%)", "# TP, VB", "# TP, varbvs",
                  "# TP, VB inter varbvs")

require(xtable)
print(xtable(tb, digits = 0, align=rep("c", 5)), include.rownames = FALSE)


## Reproducibility
#
sessionInfo()
#
# R version 3.2.3 (2015-12-10)
# Platform: x86_64-redhat-linux-gnu (64-bit)
# Running under: Red Hat Enterprise Linux Server 7.2 (Maipo)
#
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
#   [1] parallel  grid      stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] gridExtra_2.2.1 xtable_1.8-2    ggplot2_2.2.0
#
# loaded via a namespace (and not attached):
#   [1] colorspace_1.2-6 scales_0.4.1     assertthat_0.1   lazyeval_0.2.0   plyr_1.8.4
# [6] tools_3.2.3      gtable_0.1.2     tibble_1.2-12    Rcpp_0.12.7      munsell_0.4.2
