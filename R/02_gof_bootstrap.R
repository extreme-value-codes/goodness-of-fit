# 02_gof_bootstrap.R -- bootstrap goodness-of-fit tests (SLOW)
#
# Usage (from the repository root):
#   Rscript R/02_gof_bootstrap.R 508        # argument = random seed
# or set seed0 / nsim in an R session before sourcing.
#
# Because the Venice data contain many ties, the GoF statistics vary from run to run.
# The paper repeats this with 10 different seeds (123, 234, 501-508) and averages the
# p-values (see 03_gof_summary.R). Each run saves results/Results<seed>.txt
# (30 rows x 8 columns: models 7, 8, 9 stacked; columns 1-4 statistics, 5-8 p-values).

source("R/00_data.R")
source("R/rlarg_test.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) seed0 <- as.integer(args[1])
if (!exists("seed0")) seed0 <- 508
if (!exists("nsim"))  nsim  <- 5000

# Perturb by uniform noise on (-0.5, 0.5) (the rounding interval) to break ties,
# then re-sort within each year
V3 <- V2
set.seed(seed0)
for (i in seq_len(N)) {
  rr <- sum(!is.na(V2[i, ]))
  V3[i, 1:rr] <- V2[i, 1:rr] - 0.5 + runif(rr)
  V3[i, 1:rr] <- sort(as.numeric(V3[i, 1:rr]), decreasing = TRUE)
}

gof <- function(ydat, mul, muinit, data, siginit)
  rlarg.test(V3, r = 10, ydat = ydat, mul = mul, sigl = NULL, shl = NULL,
             mulink = identity, siglink = identity, shlink = identity,
             muinit = muinit, siginit = siginit, shinit = -0.001, show = FALSE,
             method = "Nelder-Mead", maxit = 10000, nsim = nsim)

# Model 7: constant location
mod7 <- gof(matrix(tt, ncol = 1), NULL, mean(V3[, 1]), V3, sqrt(var(V3[, 1])))
# Model 8: linear trend
mod8 <- gof(matrix(tt, ncol = 1), 1, c(mean(V2[, 1]), 0), V3, sqrt(var(V2[, 1])))
# Model 9: linear trend + 11-year cycle
mod9 <- gof(matrix(c(tt, cos(2*pi*tt/11), sin(2*pi*tt/11)), ncol = 3), 1:3,
            c(mean(V2[, 1]), 0, 0, 0), V3, sqrt(var(V2[, 1])))

dir.create("results", showWarnings = FALSE)
write(t(rbind(mod7, mod8, mod9)), file.path("results", paste0("Results", seed0, ".txt")), ncol = 8)
