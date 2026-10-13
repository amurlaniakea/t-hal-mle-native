#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
# Resume los .rds de mc_wcap.R. Uso: Rscript experiments/mc_wcap_summary.R <outdir>
args <- commandArgs(trailingOnly = TRUE); outdir <- if (length(args) >= 1) args[1] else "mc_out"
files <- sort(list.files(outdir, pattern = "^wcap_dgp.*\\.rds$", full.names = TRUE))
if (length(files) == 0) stop("No hay wcap_*.rds en ", outdir)
summ <- function(f) {
  o <- readRDS(f); ok <- Filter(function(x) isTRUE(x$ok), o$res); B <- length(ok); if (B < 2) return(NULL)
  tr <- o$truth; g <- o$grid; ctr <- g >= 2.5 & g <= 7
  mat <- function(get) do.call(rbind, lapply(ok, get))
  rep_rmse <- function(X) sqrt(rowMeans(sweep(X, 2, tr)^2))
  H <- mat(function(x) x$hal); dh <- rep_rmse(H)
  row <- function(nm, X) {
    E <- sweep(X, 2, tr); d <- rep_rmse(X) - dh
    data.frame(escenario = sub("\\.rds$", "", basename(f)), variante = nm, B = B,
      RMSE = sqrt(mean(E^2)), RMSE_centro = sqrt(mean(E[, ctr]^2)), RMSE_extremos = sqrt(mean(E[, !ctr]^2)),
      RMSE_a1 = sqrt(mean(E[, 1]^2)), RMSE_a9 = sqrt(mean(E[, ncol(E)]^2)),
      mediana_RMSE_rep = median(rep_rmse(X)), reps_RMSE_gt3 = sum(rep_rmse(X) > 3),
      dif_vs_HAL = mean(d), EE = sd(d) / sqrt(B), mejor_que_HAL = sum(d < 0))
  }
  vars <- names(ok[[1]]$est)
  do.call(rbind, c(list(row("HAL_plugin", H)), lapply(vars, function(v) row(v, mat(function(x) x$est[[v]])))))
}
tab <- do.call(rbind, lapply(files, summ)); options(width = 220); print(format(tab, digits = 3), row.names = FALSE)
utils::write.csv(tab, file.path(outdir, "wcap_summary.csv"), row.names = FALSE)
cat("\nGuardado:", file.path(outdir, "wcap_summary.csv"), "\n")
