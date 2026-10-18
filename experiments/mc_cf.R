#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Monte Carlo con cross-fitting K=2 para la varianza.
#
# Uso:
#   Rscript experiments/mc_cf.R <dgp 1|2|3> <normal|uniform> <n> <B> <outdir> [cores=1] [seed0=111001] [w_cap] [nuis=est|orac_gw] [c1=6] [c2=9] [n_knots_Q=10] [k_Q=auto|0|1] [sdA=2]
#   Guarda se (cross-fitting) y se_insample (residuos en muestra) del MISMO ajuste: comparacion pareada.

suppressMessages(library(thalnative))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Uso: Rscript mc_cf.R <dgp> <normal|uniform> <n> <B> <outdir> [cores] [seed0] [w_cap] [nuis]")
dgp <- as.integer(args[1]); trt <- args[2]; n <- as.integer(args[3]); B <- as.integer(args[4])
outdir <- args[5]
cores <- if (length(args) >= 6) as.integer(args[6]) else 1L
seed0 <- if (length(args) >= 7) as.integer(args[7]) else 111001L
wcap_arg <- if (length(args) >= 8) args[8] else "3"
nuis_arg <- if (length(args) >= 9) args[9] else "est"   # est | orac_gw (g y omega oraculo; Q se ajusta)
stopifnot(nuis_arg %in% c("est", "orac_gw"))
c1v <- if (length(args) >= 10) as.numeric(args[10]) else 6    # ventana del sieve: [ceil(c1 n^(1/5)), floor(c2 n^(1/5))] con k_targeting = 1
c2v <- if (length(args) >= 11) as.numeric(args[11]) else 9
stopifnot(is.finite(c1v), is.finite(c2v), c1v < c2v)
nkq <- if (length(args) >= 12) as.integer(args[12]) else 10L   # nudos por variable del HAL del resultado (n_knots_Q)
stopifnot(is.finite(nkq), nkq >= 2L)
kq_arg <- if (length(args) >= 13) args[13] else "auto"   # orden de suavidad del HAL del resultado: auto (0, salvo DGP 3) | 0 | 1
stopifnot(kq_arg %in% c("auto", "0", "1"))
sdA <- if (length(args) >= 14) as.numeric(args[14]) else 2   # desviacion tipica de A | W en el diseno normal (por defecto 2, como el paper)
stopifnot(is.finite(sdA), sdA > 0, nuis_arg == "est" || sdA == 2)   # los oraculos g0_dgp/omega_oracle asumen sd = 2
gf <- if (nuis_arg == "orac_gw") function(a, w) g0_dgp(a, w, trt) else NULL
of <- if (nuis_arg == "orac_gw") function(a) omega_oracle(a, trt) else NULL
stopifnot(dgp == 1 || dgp == 2 || dgp == 3, trt %in% c("normal", "uniform"))
wcap <- if (wcap_arg == "defecto") NULL else if (wcap_arg == "obs_max") "obs_max" else as.numeric(wcap_arg)
suffix <- paste0("_cf_cap", wcap_arg, if (seed0 != 111001L) paste0("_s", seed0) else "",
                 if (nuis_arg != "est") paste0("_", nuis_arg) else "",
                 if (c1v != 6 || c2v != 9) paste0("_c", c1v, "-", c2v) else "",
                 if (nkq != 10L) paste0("_nkq", nkq) else "",
                 if (kq_arg != "auto") paste0("_kq", kq_arg) else "",
                 if (sdA != 2) paste0("_sd", sdA) else "")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

g_w <- function(w, dgp) switch(dgp,
  5 * (w > 0),
  4 * ((w > -1) - (w > 0) + (w > 1)),
  3 * sin(7 * pi * w))

sim_paper <- function(n, dgp, trt) {
  W <- rnorm(n)
  mu <- 5 + 1.5 * W
  lo <- pnorm((0 - mu) / sdA); hi <- pnorm((10 - mu) / sdA)
  A <- if (trt == "uniform") runif(n, 0, 10) else mu + sdA * qnorm(runif(n, lo, hi))
  gw <- g_w(W, dgp)
  Y <- 3 * sin(1.5 * A) + 0.5 * A + gw + (A - 5) * gw + rnorm(n)
  list(A = A, W = W, Y = Y)
}

Eg <- switch(dgp, 5 * 0.5, 4 * (pnorm(1) - 0.5 + (1 - pnorm(1))), 0)
truth_fun <- function(a) 3 * sin(1.5 * a) + 0.5 * a + Eg * (a - 4)
grid <- seq(1, 9, length.out = 25)
tr <- truth_fun(grid)
k_Q <- if (kq_arg != "auto") as.integer(kq_arg) else if (dgp == 3) 1L else 0L

one_rep <- function(b) {
  set.seed(seed0 + b)
  d <- sim_paper(n, dgp, trt)
  t0 <- Sys.time()
  r <- tryCatch(suppressWarnings(
    thal_fit(d$A, d$W, d$Y, grid, k_targeting = 1L, k_Q = k_Q, n_proj = 500L, 
             c1 = c1v, c2 = c2v, n_knots_Q = nkq, w_cap = wcap, n_cf = 2L, g_fun = gf, omega_fun = of, seed = seed0 + b)),
    error = function(e) conditionMessage(e))
  if (is.character(r)) return(list(ok = FALSE, msg = r, b = b))
  m <- length(grid)
  Qm <- matrix(r$nuisance$Q_fun(rep(grid, each = n), rep(d$W, times = m)), n, m)
  list(ok = TRUE, b = b, psi_hat = r$psi_hat, psi_plugin = r$psi_plugin, psi_hal = colMeans(Qm),
       se = r$ci$se, se_insample = r$ci_insample$se,
       var_R = r$ic_var_parts$var_R, var_S = r$ic_var_parts$var_S, cov_RS = r$ic_var_parts$cov_RS,
       sel = r$selected, knots = r$knots, size = r$size, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

f <- file.path(outdir, sprintf("dgp%d_%s_n%d%s.rds", dgp, trt, n, suffix))
if (file.exists(f)) {
  prev <- readRDS(f)
  done <- prev$done
  results <- prev$results
  cat(sprintf("[%s] retomando desde %d / %d\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), done, B))
} else {
  done <- 0L
  results <- vector("list", B)
}

for (b in seq.int(done + 1L, B)) {
  results[[b]] <- one_rep(b)
  done <- b
  if (b %% 10 == 0 || b == B) {
    saveRDS(list(done = done, results = results, truth = tr, grid = grid, dgp = dgp, 
                 trt = trt, n = n, B = B, wcap = wcap, seed0 = seed0), f)
    cat(sprintf("[%s] %d / %d (fallos: %d)\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 
                done, B, sum(!sapply(results[1:done], function(x) x$ok))))
  }
}
cat(sprintf("[%s] FIN %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), f))
