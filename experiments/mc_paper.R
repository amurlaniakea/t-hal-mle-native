#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Monte Carlo con los DGP EXACTOS del paper (arXiv:2607.03824, seccion 7.1).
# Un escenario = (dgp, tratamiento, n). Reanudable: si el .rds de salida ya existe,
# continua desde la ultima replica guardada.
#
# Uso:
#   Rscript experiments/mc_paper.R <dgp 1|2|3> <uniform|normal> <n> <B> <outdir> [cores=1] [seed0=20261012] [w_cap]
#   w_cap: "defecto" (sqrt(n) log(n/5), por defecto), "obs_max" o un numero (p. ej. 3). Con w_cap o seed0 distintos
#   del defecto, el nombre del .rds lleva el sufijo _cap<x> y/o _s<seed0>, asi no se pisan escenarios.
#
# Requiere el paquete instalado:  R CMD INSTALL thalnative_0.0.1.tar.gz
suppressMessages(library(thalnative))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Uso: Rscript mc_paper.R <dgp> <uniform|normal> <n> <B> <outdir> [cores] [seed0]")
dgp <- as.integer(args[1]); trt <- args[2]; n <- as.integer(args[3]); B <- as.integer(args[4])
outdir <- args[5]
cores <- if (length(args) >= 6) as.integer(args[6]) else 1L
seed0 <- if (length(args) >= 7) as.integer(args[7]) else 20261012L
wcap_arg <- if (length(args) >= 8) args[8] else "defecto"
stopifnot(wcap_arg %in% c("defecto", "obs_max") || !is.na(suppressWarnings(as.numeric(wcap_arg))))
wcap <- if (wcap_arg == "defecto") NULL else if (wcap_arg == "obs_max") "obs_max" else as.numeric(wcap_arg)
suffix <- paste0(if (wcap_arg != "defecto") paste0("_cap", wcap_arg) else "",
                 if (seed0 != 20261012L) paste0("_s", seed0) else "")
stopifnot(dgp %in% 1:3, trt %in% c("uniform", "normal"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- DGP del paper (seccion 7.1) ------------------------------------------------
# W ~ N(0,1); Y = 3 sin(1.5 a) + 0.5 a + g(w) + (a - 5) g(w) + N(0,1)
# A ~ U(0,10)  o  A | W ~ N(5 + 1.5 W, 2^2) truncada a [0,10]
g_w <- function(w, dgp) switch(dgp,
  5 * (w > 0),
  4 * ((w > -1) - (w > 0) + (w > 1)),
  3 * sin(7 * pi * w))

sim_paper <- function(n, dgp, trt) {
  W <- rnorm(n)
  A <- if (trt == "uniform") runif(n, 0, 10) else {
    mu <- 5 + 1.5 * W
    lo <- pnorm((0 - mu) / 2); hi <- pnorm((10 - mu) / 2)
    mu + 2 * qnorm(runif(n, lo, hi))
  }
  gw <- g_w(W, dgp)
  Y <- 3 * sin(1.5 * A) + 0.5 * A + gw + (A - 5) * gw + rnorm(n)
  list(A = A, W = W, Y = Y)
}

# Verdad: psi0(a) = 3 sin(1.5 a) + 0.5 a + E[g(W)] (a - 4)   [DERIVADO, se comprueba abajo]
Eg <- switch(dgp, 5 * 0.5, 4 * (pnorm(1) - 0.5 + (1 - pnorm(1))), 0)
truth_fun <- function(a) 3 * sin(1.5 * a) + 0.5 * a + Eg * (a - 4)
set.seed(1); wchk <- rnorm(2e6)
stopifnot("E[g(W)] analitica no coincide con Monte Carlo" = abs(mean(g_w(wchk, dgp)) - Eg) < 0.02)

grid <- seq(1, 9, length.out = 25)
tr <- truth_fun(grid)
k_Q <- if (dgp == 3) 1L else 0L          # el paper: orden 0 en DGP 1 y 2, orden 1 en DGP 3

one_rep <- function(b) {
  set.seed(seed0 + b)
  d <- sim_paper(n, dgp, trt)
  t0 <- Sys.time()
  r <- tryCatch(suppressWarnings(
    thal_fit(d$A, d$W, d$Y, grid, k_targeting = 1L, k_Q = k_Q, n_proj = 500L, w_cap = wcap, seed = seed0 + b)),
    error = function(e) conditionMessage(e))
  if (is.character(r)) return(list(ok = FALSE, msg = r, b = b))
  m <- length(grid)
  Qm <- matrix(r$nuisance$Q_fun(rep(grid, each = n), rep(d$W, times = m)), n, m)
  list(ok = TRUE, b = b, psi_hat = r$psi_hat, psi_plugin = r$psi_plugin, psi_hal = colMeans(Qm),
       se = r$ci$se, size = r$size, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

f <- file.path(outdir, sprintf("dgp%d_%s_n%d%s.rds", dgp, trt, n, suffix))
res <- if (file.exists(f)) readRDS(f)$res else list()
cat(sprintf("[%s] dgp=%d trt=%s n=%d B=%d cores=%d | ya hechas: %d\n", format(Sys.time(), "%F %T"), dgp, trt, n, B, cores, length(res)))
todo <- setdiff(seq_len(B), vapply(res, function(x) x$b, 1L))
chunks <- split(todo, ceiling(seq_along(todo) / max(cores, 1L)))
for (ch in chunks) {
  out <- if (cores > 1L) parallel::mclapply(ch, one_rep, mc.cores = cores, mc.set.seed = FALSE) else lapply(ch, one_rep)
  res <- c(res, out)
  saveRDS(list(res = res, grid = grid, truth = tr, dgp = dgp, trt = trt, n = n, seed0 = seed0), f)
  cat(sprintf("[%s] hechas %d / %d (fallos: %d)\n", format(Sys.time(), "%F %T"), length(res), B,
              sum(!vapply(res, function(x) isTRUE(x$ok), TRUE))))
}
