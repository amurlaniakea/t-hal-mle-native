#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Resume los .rds de mc_paper.R con las metricas de la seccion 7.3 del paper.
# Uso: Rscript experiments/mc_summary.R <outdir>
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "mc_out"
files <- sort(list.files(outdir, pattern = "^dgp.*\\.rds$", full.names = TRUE))
if (length(files) == 0) stop("No hay .rds en ", outdir)

summ <- function(f) {
  o <- readRDS(f); ok <- Filter(function(x) isTRUE(x$ok), o$res); B <- length(ok)
  if (B < 2) return(NULL)
  tr <- o$truth; n <- o$n
  M <- function(nm) do.call(rbind, lapply(ok, `[[`, nm))                  # B x 25
  H <- M("psi_hat"); P <- M("psi_plugin"); L <- M("psi_hal"); SE <- M("se")
  sdmc <- apply(H, 2, sd)
  rmse <- function(X) sqrt(colMeans(sweep(X, 2, tr)^2))
  bias <- function(X) colMeans(X) - rep(0, ncol(X)) - tr
  rep_rmse <- function(X) sqrt(rowMeans(sweep(X, 2, tr)^2))
  d <- rep_rmse(H) - rep_rmse(L)                                          # T-HAL - HAL por replica
  ok_se <- stats::complete.cases(SE)
  cov_wald <- colMeans(abs(sweep(H, 2, tr)) <= 1.96 * SE, na.rm = TRUE)
  cov_mcse <- colMeans(abs(sweep(H, 2, tr)) <= 1.96 * rep(sdmc, each = B))
  data.frame(
    escenario = sub("\\.rds$", "", basename(f)), B_validas = B, B_fallidas = length(o$res) - B,
    RMSE_HAL = mean(rmse(L)), RMSE_plugin_dirigido = mean(rmse(P)), RMSE_THAL = mean(rmse(H)),
    dif_RMSE_THAL_menos_HAL = mean(d), EE_dif = sd(d) / sqrt(B), THAL_mejor_en = sum(d < 0),
    sesgo_abs_HAL = mean(abs(bias(L))), sesgo_abs_THAL = mean(abs(bias(H))),
    frac_sesgo_despreciable = mean(abs(bias(H) / sdmc) <= 1 / log(n)),
    cobertura_MCSE = mean(cov_mcse), cobertura_Wald = mean(cov_wald),
    ancho_MCSE = mean(3.92 * sdmc), ancho_Wald = mean(3.92 * colMeans(SE, na.rm = TRUE)),
    tam_mediano = median(vapply(ok, function(x) x$size, 1)), seg_por_rep = mean(vapply(ok, function(x) x$secs, 1)),
    check.names = FALSE)
}
tab <- do.call(rbind, lapply(files, summ))
options(width = 200); print(format(tab, digits = 3), row.names = FALSE)
utils::write.csv(tab, file.path(outdir, "summary.csv"), row.names = FALSE)
cat("\nGuardado:", file.path(outdir, "summary.csv"), "\n")
cat("Referencia del paper (7.4-7.5): T-HAL mejor RMSE que HAL; cobertura MC-SE casi nominal;\n",
    "Wald: cerca de 0.95 con tratamiento uniforme, infracobertura con normal (sobre todo DGP 2).\n", sep = "")
