#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Experimento de diagnostico: ¿por que T-HAL falla con tratamiento normal?
# En cada replica (DGP exactos del paper, requiere el paquete con simulate_data() del paper):
#   - se ajusta Q (HAL nativo) UNA vez;
#   - se evalua T-HAL con g y omega ESTIMADAS ("est") y con g y omega ORACULO ("orac"),
#     cada una con distintos topes de peso w_cap: defecto (sqrt(n) log(n/5)), "obs_max", 10 y 5.
#   Si con g oraculo tambien falla, el problema no es la estimacion de g_n sino el targeting
#   con pesos omega/g enormes en las colas.
#
# Uso: Rscript experiments/mc_wcap.R <dgp 1|2|3> <n> <B> <outdir> [cores=1] [seed0=20261012] [caps] [nuis]
#   caps: lista separada por comas; "defecto", "obs_max" o un numero. Por defecto: defecto,obs_max,10,5
#   nuis: "both" (g/omega estimadas y oraculo; por defecto) o "est" (solo estimadas, mas rapido)
suppressMessages(library(thalnative))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Uso: Rscript mc_wcap.R <dgp> <n> <B> <outdir> [cores] [seed0]")
dgp <- as.integer(args[1]); n <- as.integer(args[2]); B <- as.integer(args[3]); outdir <- args[4]
cores <- if (length(args) >= 5) as.integer(args[5]) else 1L
seed0 <- if (length(args) >= 6) as.integer(args[6]) else 20261012L
trt <- "normal"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
grid <- seq(1, 9, length.out = 25); tr <- true_drc(grid, dgp)
k_Q <- if (dgp == 3) 1L else 0L
caps_arg <- if (length(args) >= 7) args[7] else "defecto,obs_max,10,5"
nuis_arg <- if (length(args) >= 8) args[8] else "both"
stopifnot(nuis_arg %in% c("both", "est"))
tok <- strsplit(caps_arg, ",")[[1]]
caps <- setNames(vector("list", length(tok)), ifelse(grepl("^[0-9.]+$", tok), paste0("cap", tok), tok))
for (i in seq_along(tok)) caps[i] <- list(if (tok[i] == "defecto") NULL else if (tok[i] == "obs_max") "obs_max" else as.numeric(tok[i]))

one_rep <- function(b) {
  set.seed(seed0 + b)
  d <- simulate_data(n, dgp, trt)
  out <- tryCatch(suppressWarnings({
    Q_fun <- fit_Q_hal_native(d$A, d$W, d$Y, k = k_Q, n_knots = 10L, nfolds = 5L, seed = seed0 + b)
    nuis <- list(
      est  = list(g = fit_g_native(d$A, d$W, n_bins = 5L), om = fit_omega_kde_native(d$A)),
      orac = list(g = function(a, w) g0_dgp(a, w, trt),    om = function(a) omega_oracle(a, trt)))
    if (nuis_arg == "est") nuis$orac <- NULL
    m <- length(grid)
    hal <- colMeans(matrix(Q_fun(rep(grid, each = n), rep(d$W, times = m)), n, m))
    est <- list()
    for (nn in names(nuis)) for (cc in names(caps)) {
      r <- thal_core_native(d$A, d$W, d$Y, Q_fun, nuis[[nn]]$g, nuis[[nn]]$om, grid,
                            k = 1L, n_knots = 50L, n_proj = 500L, w_cap = caps[[cc]])
      est[[paste(nn, cc, sep = "_")]] <- r$psi_hat
    }
    list(ok = TRUE, b = b, hal = hal, est = est)
  }), error = function(e) list(ok = FALSE, b = b, msg = conditionMessage(e)))
  out
}

f <- file.path(outdir, sprintf("wcap_dgp%d_n%d_s%d.rds", dgp, n, seed0))
res <- if (file.exists(f)) readRDS(f)$res else list()
cat(sprintf("[%s] wcap dgp=%d n=%d B=%d cores=%d | ya hechas: %d\n", format(Sys.time(), "%F %T"), dgp, n, B, cores, length(res)))
todo <- setdiff(seq_len(B), vapply(res, function(x) x$b, 1L))
for (ch in split(todo, ceiling(seq_along(todo) / max(cores, 1L)))) {
  out <- if (cores > 1L) parallel::mclapply(ch, one_rep, mc.cores = cores, mc.set.seed = FALSE) else lapply(ch, one_rep)
  res <- c(res, out)
  saveRDS(list(res = res, grid = grid, truth = tr, dgp = dgp, n = n, seed0 = seed0), f)
  cat(sprintf("[%s] hechas %d / %d (fallos: %d)\n", format(Sys.time(), "%F %T"), length(res), B,
              sum(!vapply(res, function(x) isTRUE(x$ok), TRUE))))
}
