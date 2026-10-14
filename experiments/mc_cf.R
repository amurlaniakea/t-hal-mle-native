#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Monte Carlo con cross-fitting K=2 para la varianza.
#
# Uso:
#   Rscript experiments/mc_cf.R <dgp 1|2|3> <normal> <n> <B> <outdir> [cores=1] [seed0=111001] [w_cap]

suppressMessages(library(thalnative))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Uso: Rscript mc_cf.R <dgp> <normal> <n> <B> <outdir> [cores] [seed0] [w_cap]")
dgp <- as.integer(args[1]); trt <- args[2]; n <- as.integer(args[3]); B <- as.integer(args[4])
outdir <- args[5]
cores <- if (length(args) >= 6) as.integer(args[6]) else 1L
seed0 <- if (length(args) >= 7) as.integer(args[7]) else 111001L
wcap_arg <- if (length(args) >= 8) args[8] else "3"
stopifnot(dgp == 1 || dgp == 2 || dgp == 3, trt == "normal")
wcap <- if (wcap_arg == "defecto") NULL else if (wcap_arg == "obs_max") "obs_max" else as.numeric(wcap_arg)
suffix <- paste0("_cf_cap", wcap_arg, if (seed0 != 111001L) paste0("_s", seed0) else "")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

g_w <- function(w, dgp) switch(dgp,
  5 * (w > 0),
  4 * ((w > -1) - (w > 0) + (w > 1)),
  3 * sin(7 * pi * w))

sim_paper <- function(n, dgp, trt) {
  W <- rnorm(n)
  mu <- 5 + 1.5 * W
  lo <- pnorm((0 - mu) / 2); hi <- pnorm((10 - mu) / 2)
  A <- mu + 2 * qnorm(runif(n, lo, hi))
  gw <- g_w(W, dgp)
  Y <- 3 * sin(1.5 * A) + 0.5 * A + gw + (A - 5) * gw + rnorm(n)
  list(A = A, W = W, Y = Y)
}

Eg <- switch(dgp, 5 * 0.5, 4 * (pnorm(1) - 0.5 + (1 - pnorm(1))), 0)
truth_fun <- function(a) 3 * sin(1.5 * a) + 0.5 * a + Eg * (a - 4)
grid <- seq(1, 9, length.out = 25)
tr <- truth_fun(grid)
k_Q <- if (dgp == 3) 1L else 0L

one_rep <- function(b) {
  set.seed(seed0 + b)
  d <- sim_paper(n, dgp, trt)
  t0 <- Sys.time()
  r <- tryCatch(suppressWarnings(
    thal_fit(d$A, d$W, d$Y, grid, k_targeting = 1L, k_Q = k_Q, n_proj = 500L, 
             w_cap = wcap, n_cf = 2L, seed = seed0 + b)),
    error = function(e) conditionMessage(e))
  if (is.character(r)) return(list(ok = FALSE, msg = r, b = b))
  m <- length(grid)
  Qm <- matrix(r$nuisance$Q_fun(rep(grid, each = n), rep(d$W, times = m)), n, m)
  list(ok = TRUE, b = b, psi_hat = r$psi_hat, psi_plugin = r$psi_plugin, psi_hal = colMeans(Qm),
       se = r$ci$se, size = r$size, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
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
