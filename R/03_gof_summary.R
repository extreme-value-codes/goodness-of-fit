# 03_gof_summary.R -- average GoF p-values across all saved runs in results/
# Creates: gof_avg (list of 3 matrices, 10 x 4), gof_tab (data frame) -- or NULL if no results found

fn <- list.files("results", pattern = "^Results.*\\.txt$", full.names = TRUE)
gof_avg <- NULL; gof_tab <- NULL

if (length(fn) > 0) {
  tot <- list(m7 = 0, m8 = 0, m9 = 0)
  for (f in fn) {
    M <- matrix(scan(f, quiet = TRUE), ncol = 8, byrow = TRUE)
    tot$m7 <- tot$m7 + M[1:10,  5:8]
    tot$m8 <- tot$m8 + M[11:20, 5:8]
    tot$m9 <- tot$m9 + M[21:30, 5:8]
  }
  gof_avg <- lapply(tot, function(x) x / length(fn))
  names(gof_avg) <- c("Model 7: constant", "Model 8: linear trend", "Model 9: trend + 11-yr cycle")
  
  nm <- c("R", "KS", "CvM", "AD")
  gof_tab <- data.frame(k = 1:10, signif(gof_avg[[1]], 2), signif(gof_avg[[2]], 2), signif(gof_avg[[3]], 2))
  names(gof_tab) <- c("k", paste0("M7.", nm), paste0("M8.", nm), paste0("M9.", nm))
}

# LaTeX-ready rows (optional):
# if (!is.null(gof_tab)) for (i in 1:10) cat(paste(gof_tab[i, -1], collapse = "&"), "\\\\\n")
