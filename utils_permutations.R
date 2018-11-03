# counts the number of permutation analyses available for empirical null estimation
#
counts_perm_runs <- function(vec_dir) {
  n_perm <- 0
  for(curr_dir in vec_dir) {
    n_perm <- n_perm + length(list.files(path = curr_dir, pattern = "*.RData"))
    print(n_perm)
  }
}

# loads the posterior inclusion probability matrices obtained by applying LOCUS 
# or varbvs on permuted data
#
load_perm_ppi <- function(vec_dir, meth, n_cpus = 1) {
  
  stopifnot(meth %in% c("vb", "varbvs"))
  
  list_ppi <- NULL
  
  for(curr_dir in vec_dir) {
    
    cat(paste("Loads files in ", curr_dir, "\n", sep = ""))
    cat("...\n")
    
    files <- list.files(path = curr_dir, pattern = "*.RData")
    
    require(parallel)
    list_perm <- mclapply(files, function(f) { 
      load(file.path(curr_dir, f))
      if (meth == "vb") gam_vb
      else ppi
    }, mc.cores = n_cpus)
    
    list_ppi <- c(list_ppi, list_perm)
    
  }
  
  names(list_ppi) <- paste("perm_", 1:length(list_ppi), sep = "")
  list_ppi
  
} 



# estimates permutation-based false positive rates (FDR) from a grid of thresholds
# on the posterior inclusion probabilities
#
estimate_FDR <- function(list_perm_ppi, actual_ppi, vec_thres_ppi, type, 
                         n_cpus = 1, bool_tab = TRUE, bool_adapt_thres = FALSE) {
  
  stopifnot(type %in% c("mean", "median"))
  
  if (bool_adapt_thres) {
    
    med <- median(actual_ppi)
    vec_med_k <- apply(actual_ppi, 2, median)
    actual_ppi <- sweep(actual_ppi, 2, vec_med_k, FUN = "/") * med
    list_perm_ppi <- lapply(list_perm_ppi, 
                            function(perm_ppi) sweep(perm_ppi, 2, vec_med_k, FUN = "/") * med)
    
  } 
  
  list_tab_declared <-  vec_n_declared <- vec_FDR <- NULL
  for (thres_ppi in vec_thres_ppi) {
    
    cat(paste("Computes estimated FDR for threshold ", thres_ppi, "\n", sep = ""))
    cat("...\n")
    
    require(parallel)
    list_n_FP <- mclapply(list_perm_ppi, 
                          function(perm_ppi) sum(perm_ppi > thres_ppi), mc.cores = n_cpus)
    
    if (type == "mean") m_n_FP <- mean(unlist(list_n_FP))
    else m_n_FP <- median(unlist(list_n_FP))
    
    if (bool_tab) {
      if (is.null(dim(actual_ppi))) {
        tab_declared <- which(actual_ppi > thres_ppi)
        n_declared <- length(tab_declared)
      } else { 
        tab_declared <- which(actual_ppi > thres_ppi, arr.ind = TRUE)
        n_declared <- nrow(tab_declared)
      }
    } else {
      n_declared <- sum(actual_ppi > thres_ppi)
    }
    
    FDR <- m_n_FP / n_declared
    
    if(bool_tab) list_tab_declared <- c(list_tab_declared, list(tab_declared))
    
    vec_n_declared <- c(vec_n_declared, n_declared)
    vec_FDR <- c(vec_FDR, FDR)
  }
  
  # excludes impossible FDR caused for instance by no discovery in the actual
  # dataset and, possibly, a mean/median number of FP equal to 0
  bool_rm <- is.nan(vec_FDR) | is.infinite(vec_FDR) | vec_FDR > 1
  vec_thres_ppi <- vec_thres_ppi[!bool_rm]
  vec_FDR <- vec_FDR[!bool_rm]
  vec_n_declared <- vec_n_declared[!bool_rm]
  if(bool_tab) list_tab_declared <- list_tab_declared[!bool_rm]
  
  if (length(vec_FDR) > 0) {
    names(vec_n_declared) <- names(vec_FDR) <- vec_thres_ppi
    if(bool_tab) names(list_tab_declared) <- vec_thres_ppi
  }
  
  list("list_tab_declared" = list_tab_declared, "vec_n_declared" = vec_n_declared, 
       "vec_FDR" = vec_FDR)
  
}


# in case a spline cannot be fitted to the estimated FDR, perpare_hm finds the 
# threshold corresponding to the nearest estimated FDR value
#
find_closest_FDR_thres <- function(vec_FDR_avail, vec_FDR_desired) {
  sapply(vec_FDR_desired, function(FDR) which.min(abs(vec_FDR_avail - FDR)))
}


# finds posterior inclusion probabilities thresholds corresponding to specific FDR
# values by fitting a spline to previously estimated FDR-threshold values
prepare_hm <- function(vec_FDR, vec_thres, vec_FDR_desired, ppi, title, 
                       out_dir = NULL, spar = NULL, bool_save = FALSE, bool_log = FALSE) {
  
  if (bool_save) stopifnot(!is.null(out_dir))
  if (!is.null(spar) & length(unique(vec_FDR)) < 4) { 
    spar <- NULL
    warning("At least 4 unique values in vec_FDR are needed to fit a spline. Finds closest FDR value instead.")
  }
  
  if (is.null(spar)) { 
    # matches the values of vec_FDR_desired with the closest ones in vec_FDR and
    # takes the corresponding thresholds
    ind_select <- find_closest_FDR_thres(vec_FDR, vec_FDR_desired)
    FDR_select <- vec_FDR[ind_select]
    tab <- which(ppi > vec_thres[ind_select[1]], arr.ind = TRUE)
    vec_thres_select <- vec_thres[ind_select]
    
  } else { # fits a spline
    
    if (bool_log) {
      vec_thres <- log(vec_thres)
      ylabel <- "Threshold on log PPI"
    } else {
      ylabel <- "Threshold on PPI"
    } 
    
    spl <- smooth.spline(x = vec_FDR, y = vec_thres, spar = spar, tol = 1e-6) 
    spl_pred <- predict(spl, vec_FDR_desired)
    
    if (bool_log) vec_thres_select <- exp(spl_pred$y)
    else vec_thres_select <- spl_pred$y
    
    FDR_select <- vec_FDR_desired 
    names(FDR_select) <- vec_thres_select
    tab <- which(ppi > vec_thres_select[1], arr.ind = TRUE)
    
    if (bool_save)
      png(paste(out_dir, "/FDR_spline_", title, ".png", sep=""), w=3500, h=3000, 
          res=500, type="cairo")
    
    plot(vec_FDR, vec_thres, main = "Permutation-based FDR, cubic spline fit", 
         xlab = "FDR", ylab = ylabel)
    lines(spl, col="red", lwd = 2)
    
    if (bool_save) dev.off()
  }
  
  print("FDR_select: ")
  print(FDR_select)
  
  list("vec_thres_select"= vec_thres_select, "FDR_select" = FDR_select, "tab" = tab)
}



# gathers settings for the comparison pattern overlaid on the FDR heatmap
#
make_obj_compare <- function(mat, thres, text, symb = "x", color = "red", cexc = 6) {
  
  create_named_list <- function(...) { 
    setNames(list(...), as.character(match.call()[-1])) 
  }
  
  obj <- create_named_list(mat, thres, text, symb, color, cexc)
  class(obj) <- "obj_comp"
  obj
}



# plots a heatmap for association evidence based estimated FDR. 
# vec_FDR must contain FDRs in decreasing order and must be named with the 
# corresponding PPI thresholds.
#
make_FDR_heatmap <- function(ppi, ind_s, ind_t, vec_FDR, x_max, y_max, title,
                             out_dir = NULL, title_file = NULL, bool_save = FALSE,
                             bool_repl = FALSE, obj_comp1 = NULL, obj_comp2 = NULL) {
  
  if(bool_save) stopifnot(!is.null(out_dir) & !is.null(title_file))
  
  if(!all(c(class(obj_comp1), class(obj_comp2)) %in% c("obj_comp", "NULL")))
    stop("obj_comp1 and obj_comp2, if provided, must be of class `obj_comp'.")
  
  if(class(obj_comp1) == "NULL" & class(obj_comp2) == "obj_comp")
    stop("if only one `obj_comp' object is provided, it must be passed to obj_comp1.")
  
  ppi <- ppi[ind_s, ind_t]
  n_thres <- length(vec_FDR)
  vec_FDR_thres <- as.numeric(names(vec_FDR))
  mat_fdr <- matrix(0, nrow = nrow(ppi), ncol = ncol(ppi))
  
  # assigns associations based on the FDR threshold at which they are declared
  vec_pres <- rep(TRUE, n_thres+1) 
  for(i in 2:n_thres) {
    bool_int <- ppi > vec_FDR_thres[i-1] & ppi <= vec_FDR_thres[i]
    mat_fdr[bool_int] <- i-1
    if (all(!bool_int)) vec_pres[i] <- FALSE
  }
  
  bool_int <- ppi > vec_FDR_thres[n_thres]
  mat_fdr[bool_int] <- n_thres
  if (all(!bool_int)) vec_pres[n_thres+1] <- FALSE
  
  rownames(mat_fdr) <- rownames(ppi)
  colnames(mat_fdr) <- colnames(ppi)
  
  mat_fdr_melt <- as.data.frame.table(mat_fdr)
  if (is.null(obj_comp1)) {
    
    colnames(mat_fdr_melt) <- c("X", "Y", "values") 
    
  } else { # one or two additional patterns to overlay
    
    mat_comp <- obj_comp1$mat[ind_s, ind_t]
    mat_comp_melt <- as.data.frame.table(mat_comp)
    mat_fdr_melt <- cbind(mat_fdr_melt, 
                          as.numeric(as.character(rep("", nrow(mat_fdr_melt)))))
    colnames(mat_fdr_melt) <- c("X", "Y", "values", "comp1") 
    mat_fdr_melt$comp1[mat_comp_melt$Freq <= obj_comp1$thres] <- obj_comp1$symb
    mat_fdr_melt$comp1[mat_comp_melt$Freq > obj_comp1$thres] <- ""
    
    if (!is.null(obj_comp2)) {
      mat_comp <- obj_comp2$mat[ind_s, ind_t]
      mat_comp_melt <- as.data.frame.table(mat_comp)
      mat_fdr_melt <- cbind(mat_fdr_melt, 
                            as.numeric(as.character(rep("", nrow(mat_fdr_melt)))))
      colnames(mat_fdr_melt)[5] <- "comp2" 
      mat_fdr_melt$comp2[mat_comp_melt$Freq <= obj_comp2$thres] <- obj_comp2$symb
      mat_fdr_melt$comp2[mat_comp_melt$Freq > obj_comp2$thres] <- ""
    }
  } 
  
  mat_fdr_melt$X = with(mat_fdr_melt, factor(X, levels = rev(levels(X))))
  
  lab_leg <- c(paste("> ", vec_FDR[1] * 100, sep = ""), 
               paste(vec_FDR[-n_thres] * 100," - ", vec_FDR[-1] * 100, sep=""), 
               paste("< ", vec_FDR[n_thres] * 100, sep = ""))[vec_pres]
  col_val <- gray.colors(n_thres + 1, start = 1.0, end = 0.0, gamma = 2.2, alpha = NULL)[vec_pres]
  
  require(ggplot2)
  require(grid)
  
  # central panel
  #
  if(bool_repl) leg_pos <- c(-0.07, 1.075)
  else leg_pos <- c(-0.05, -0.15)
  
  hmap <- ggplot(data = mat_fdr_melt, aes(x = X, y = Y)) + theme_bw() +
    geom_tile(aes(fill = factor(values))) +
    scale_fill_manual(values = col_val, name = "perm. FDR %", labels=lab_leg) + 
    coord_flip() + ggtitle("") + labs(y = "Declared metabolites" , x = "Declared SNPs") + 
    theme(title = element_text(size = 17, face = "bold"), 
          axis.text.x = element_text(angle=90, vjust=0.5, size=10.8), 
          axis.text.y = element_text(angle=0, vjust=0.5, size=8.5),
          axis.title = element_text(size=12, face="bold"), axis.line = element_blank(), 
          panel.background = element_blank(), panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.background = element_blank(), plot.margin=unit(c(1.1,1,0.4,1), "cm"),
          legend.title = element_text(size=11, face= "bold")) + 
    theme(legend.position = leg_pos, legend.key.size = unit(0.8, "lines"), 
          legend.direction="vertical")
  if (!is.null(obj_comp1)) 
    hmap <- hmap + annotate("text", x = mat_fdr_melt$X, y = mat_fdr_melt$Y, 
                            label = mat_fdr_melt$comp1, col = obj_comp1$color, 
                            size = obj_comp1$cexc, fontface = "bold") 
  if(!is.null(obj_comp2))
    hmap <- hmap + annotate("text", x = mat_fdr_melt$X, y = mat_fdr_melt$Y, 
                            label = mat_fdr_melt$comp2, col = obj_comp2$color, 
                            size = obj_comp2$cexc, fontface = "bold")
  
  # upper panel
  #
  c_y <- colSums(mat_fdr > 0)
  counts_y <- as.data.frame.table(as.matrix(c_y))[-2]
  counts_y <- cbind(counts_y)
  colnames(counts_y) <- c("X", "values")
  counts_y$X = with(counts_y, factor(X))
  
  if (bool_repl) add_margin <- 0.7
  else add_margin <- 0
  
  q_marg_y <- ggplot(counts_y, aes(X, values)) + 
    ggtitle(paste(title, ", FDR-based evidence", sep="")) + 
    geom_bar(stat="identity", position="dodge") +
    scale_y_continuous(breaks = c(0, floor(y_max/2), y_max), limits = c(0, y_max)) + 
    theme_bw() + labs(y = "", x = "# declared SNPs per metabolite") + 
    theme(title = element_text(size = 18, face="bold"), axis.line=element_blank(),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_blank(), axis.text.y = element_text(size = 11),
          axis.ticks = element_blank(), panel.border = element_rect(colour = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.background = element_blank(), legend.title = element_text(size = 12),
          plot.margin = unit(c(1.5, 0.97, -1.5, 2.9-add_margin), "cm"))
  
  
  # right panel
  #
  c_x <- rowSums(mat_fdr > 0)
  counts_x <- as.data.frame.table(as.matrix(c_x))[-2]
  counts_x <- cbind(counts_x)
  colnames(counts_x) <- c("X", "values")
  counts_x$X = with(counts_x, factor(X, levels = rev(levels(X))))
  
  if (bool_repl) plot_margin <- c(1.13, 1.6, 2.05,-0.5)
  else plot_margin <- c(1.18,1.6,3.05,-0.5)
  
  q_marg_x <- ggplot(counts_x, aes(X, values)) + coord_flip() + ggtitle("") + 
    geom_bar(stat="identity", position="dodge") + 
    scale_y_continuous(breaks = c(0, floor(x_max/2), x_max), limits = c(0, x_max)) + 
    theme_bw() + labs(y = "", x = "# declared metabolites per SNP") + 
    theme(title = element_text(size = 15, face = "bold"), axis.line = element_blank(),
          axis.text.x = element_text(size=11), 
          axis.title=element_text(size = 12, face="bold"),
          axis.text.y = element_blank(), axis.title.y = element_text(angle = -90),
          axis.ticks = element_blank(), panel.border = element_rect(colour = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.background = element_blank(), legend.title = element_text(size=12),
          plot.margin = unit(plot_margin, "cm"))
  
  
  require(gridExtra)
  
  if (bool_save) {
    if (bool_repl) hei <- 6800
    else hei <- 4250
    png( paste(out_dir, "/FDR_", title_file, ".png", sep=""), w=5900, h=hei, res=500, type="cairo" )
  }
  r <- rectGrob(gp=gpar(fill="white", col="white"), height=0.01)
  if (bool_repl) wid <- 0.25
  else wid <- 0.3
  grid.arrange(q_marg_y, r, hmap, q_marg_x, ncol=2, widths = c(1.8, wid), heights = c(wid-0.15, 1))
  if (bool_repl) add <- 0.065
  else add <- 0
  if (!is.null(obj_comp1)) 
    grid.text(obj_comp1$text, x = unit(0.838 + add / 2, "npc"), 
              y = unit(0.11 - add, "npc"), just = c("left", "bottom"), 
              gp = gpar(fontsize = 11, col = obj_comp1$color, fontface="bold"))
  if (!is.null(obj_comp2)) 
    grid.text(obj_comp2$text, x = unit(0.838 + add / 2, "npc"), 
              y = unit(0.05 - add, "npc"), just = c("left", "bottom"), 
              gp = gpar(fontsize = 11, col = obj_comp2$color, fontface="bold"))
  
  if (bool_save) dev.off()
  
}
