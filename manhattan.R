rm(list=ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

main_dir <- file.path(CORE_DIR, "work/mQTL_analysis_example/")
setwd(main_dir)

## loads results from the replicated mQTL data analysis.
#
results_dir <- file.path(CORE_DIR, "results/Repl_mQTL_analysis_seed_123/")

load(file.path(results_dir, "output_vb.RData"))
load(file.path(results_dir, "output_varbvs.RData"))


## generic function to create manhattan plots from any type of association summaries.
## The function can also use annotation data (if available) to arrange SNPs by
## chromosome. Specific SNPs can be spotted with colors and their names can be
## displayed.
#
make_manhattan <- function(scores_x, name, thres = 5e-7, results_dir = NULL,
                           log_10 = TRUE, minus = FALSE, ytitle = "score",
                           size_text = 0.5, path_chr = NULL, vec_pat = NULL,
                           vec_col = c("black", "red"), leg = NULL, add = 0.0) {
  
  if (log_10) scores_x <- log10(scores_x)
  if (minus) scores_x <- - scores_x
  
  if (is.null(path_chr)) {
    
    data <- as.data.frame(scores_x, stringsAsFactors = FALSE)
    if (is.null(vec_pat)) spot <- FALSE
    else spot <- vec_pat
    
  } else {
    
    load(path_chr)
    
    pos <- positions[names(positions) %in% names(scores_x)]
    if (sum(names(pos) != names(scores_x)))
      stop("Names of SNP positions and names of scores do not match. \n")
    
    chr <- chromosomes[names(chromosomes)%in%names(scores_x)]
    if (sum(names(chr) != names(scores_x)))
      stop("Names in chromosome vector and names of scores do not match. \n")
    
    chr_string <- paste("chr", as.character(chr), sep = "")
    rsID <- names(scores_x)
    
    data <- as.data.frame(cbind(chr_string, pos, rsID, scores_x),
                          stringsAsFactors = FALSE)
    
    spot <- as.numeric(gsub(data$chr_string[data$chr_string != "chrX"],
                            pattern = "chr", replacement = "")) %% 2
    spot <- c(spot, rep(1, sum(data$chr=="chrX")))
    
    t_chr <- table(chr)
    incr <- 0
    vec_at <- NULL
    for( ch in 1:22 ){
      incr <- incr + floor(t_chr[names(t_chr) == ch] / 2)
      vec_at <- c(vec_at, incr)
      incr <- incr + ceiling(t_chr[names(t_chr) == ch] / 2)
    }
    incr <- incr + floor(t_chr[names(t_chr) == "X"] / 2)
    vec_at <- c(vec_at, incr)
  }
  
  spot_thres <- ifelse(as.numeric(as.character(data[,"scores_x"])) > thres, TRUE, FALSE)
  
  at_y <- seq(min(as.numeric(as.character(data$scores_x))),
              max(as.numeric(as.character(data$scores_x))) + 0.1, length.out = 4)
  if(minus)
    sign_at_y <- - at_y
  else
    sign_at_y <- at_y
  
  if(log_10)
    lab_y <- format(10^sign_at_y, digits = 3, scientific = TRUE)
  else
    lab_y <- format(sign_at_y, digits = 3, scientific = TRUE)
  
  if (!is.null(results_dir)) {
    png(paste(results_dir,"manhattan_plot_", name, ".png", sep=""), w=3250, h=2500,
        res=500, type="cairo")
  }
  
  plot(x = 1:nrow(data), y = as.numeric(as.character(data$scores_x)),
       ylim = c(min(as.numeric(as.character(data$scores_x))),
                max(as.numeric(as.character(data$scores_x))) + add),
       col = vec_col[as.factor(spot)], pch=20, axes=FALSE,
       main = paste("Manhattan plot ", name, sep = ""), xlab="SNPs", ylab=ytitle)
  
  if (any(spot_thres)) {
    if (is.null(path_chr))
      text(labels = names(scores_x)[spot_thres], x = c(1:nrow(data))[spot_thres],
           y = as.numeric(as.character(data$scores_x))[spot_thres], pos = 3,
           cex = size_text)
    else
      text(labels = as.character(data$rsID[spot_thres]), x = c(1:nrow(data))[spot_thres],
           y = as.numeric(as.character(data$scores_x))[spot_thres], pos = 3,
           cex = size_text)
  }
  if (!is.null(path_chr))
    axis(1, at=vec_at, labels=c( paste("chr", 1:22, sep = ""), "chrX"), cex.axis=0.8, las=2)
  axis(2, at=at_y, labels= lab_y,  cex.axis=0.8)
  if (!is.null(leg))
    legend("topleft",leg, pch=20, col=vec_col[2], bg=vec_col[2], bty="n")
  if (!is.null(results_dir))
    dev.off()
}



## Manhattan plot for LOCUS
#
make_manhattan(out_vb$om_vb, name = "VB", thres = log10(3.75e-5),
               results_dir = results_dir, ytitle = "E_q(omega | y)",
               vec_pat = rowSums(dat$pat) > 0, leg = "active", add = 0.05)


## Manhattan plot for varbvs
#
make_manhattan(rowSums(ppi_varbvs), name = "varbvs", thres = log10(1.1),
               results_dir = results_dir, ytitle = "rowsums(PPI)",
               vec_pat = rowSums(dat$pat) > 0, leg = "active", add = 0.35)


## Reproducibility
#
sessionInfo()
#
# R version 3.3.1 (2016-06-21)
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
#   [1] stats     graphics  grDevices utils     datasets  methods   base
#
# loaded via a namespace (and not attached):
#   [1] tools_3.3.1

