# 01_fit_models.R -- fit r-largest models, likelihood ratio tests, fitted curves
# Requires 00_data.R to have been run. Creates: mod1..mod5, tab_params, tab_lrt, pred2, pred3, pred5
#
# Models (location parameter mu_t; sigma and xi constant):
#   1: constant   2: linear trend   3: quadratic trend
#   4: linear + 18.62-year sinusoid   5: linear + 11-year sinusoid

fit <- function(ydat = NULL, mul = NULL, muinit)
  rlarg.fit(V2, r = 10, ydat = ydat, mul = mul, sigl = NULL, shl = NULL,
            mulink = identity, siglink = identity, shlink = identity,
            muinit = muinit, siginit = sqrt(var(V2[, 1])), shinit = -0.001,
            show = FALSE, method = "Nelder-Mead", maxit = 10000)

m0 <- mean(V2[, 1])
mod1 <- fit(muinit = m0)
mod2 <- fit(matrix(tt, ncol = 1), 1, c(m0, 0))
mod3 <- fit(matrix(c(tt, tt^2), ncol = 2), 1:2, c(m0, 0, 0))
mod4 <- fit(matrix(c(tt, cos(2*pi*tt/18.62), sin(2*pi*tt/18.62)), ncol = 3), 1:3, c(m0, 0, 0, 0))
mod5 <- fit(matrix(c(tt, cos(2*pi*tt/11),    sin(2*pi*tt/11)),    ncol = 3), 1:3, c(m0, 0, 0, 0))

# Parameter estimates, Model 2
tab_params <- round(cbind(mod2$mle, mod2$se, mod2$mle / mod2$se,
                          2 * pnorm(-abs(mod2$mle / mod2$se))), 4)
dimnames(tab_params) <- list(c("mu (intercept)", "mu (trend)", "sigma", "xi"),
                             c("Estimate", "Std. error", "z", "p-value"))

# Likelihood ratio tests
D <- c(2 * (mod1$nllh - mod2$nllh), 2 * (mod2$nllh - mod3$nllh),
       2 * (mod2$nllh - mod4$nllh), 2 * (mod2$nllh - mod5$nllh))
df <- c(1, 1, 2, 2)
tab_lrt <- data.frame(
  Comparison = c("2 vs 1 (linear trend)", "3 vs 2 (quadratic term)",
                 "4 vs 2 (18.62-yr cycle)", "5 vs 2 (11-yr cycle)"),
  Deviance = round(D, 3), df = df,
  `p-value` = signif(pchisq(D, df, lower.tail = FALSE), 3), check.names = FALSE)

# Fitted median curves for models 2, 3, 5
pred2 <- matrix(c(rep(1, N), tt, rep(0, N), rep(0, N)), ncol = 4) %*% mod2$mle +
  (mod2$mle[3] * (log(2)^(-mod2$mle[4]) - 1) / mod2$mle[4])
pred3 <- matrix(c(rep(1, N), tt, tt^2, rep(0, N), rep(0, N)), ncol = 5) %*% mod3$mle +
  (mod3$mle[4] * (log(2)^(-mod3$mle[5]) - 1) / mod3$mle[5])
pred5 <- matrix(c(rep(1, N), tt, cos(2*pi*tt/11), sin(2*pi*tt/11), rep(0, N), rep(0, N)), ncol = 6) %*% mod5$mle +
  (mod5$mle[5] * (log(2)^(-mod5$mle[6]) - 1) / mod5$mle[6])

# Standard ismev diagnostics (not used in the paper):
# rlarg.diag(mod2)
