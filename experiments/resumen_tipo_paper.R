#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Resumen de los .rds de mc_paper.R / mc_cf.R con las metricas de la seccion 7.3 del paper (arXiv:2607.03824):
# RMSE y sesgo PUNTUALES en cada dosis, cociente sesgo / MC-SE (banda +-1/log(n)), cobertura con el MC-SE (sd Monte Carlo)
# y cobertura de Wald (solo T-HAL), resumidos como en los diagramas de caja del paper: mediana y cuartiles sobre las 25 dosis.
# Estimadores: HAL-MLE (psi_hal), T-HAL-MLE (Plugin) (psi_plugin) y T-HAL-MLE (psi_hat).
# Uso: Rscript experiments/resumen_tipo_paper.R <outdir>
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "mc_out_cf"
files <- sort(list.files(outdir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) == 0) stop("No hay .rds en ", outdir)
q3 <- function(x) quantile(x, c(.25, .5, .75), names = FALSE)
rows <- list()
for (f in files) {
  o <- readRDS(f); res <- if (!is.null(o$results)) o$results else o$res; tr <- o$truth; n <- o$n
  ok <- Filter(function(x) isTRUE(x$ok) && !is.null(x$psi_hal), res)
  ok <- Filter(function(x) is.null(o$B) || x$b <= o$B, ok)
  ok <- ok[!duplicated(vapply(ok, function(x) x$b, 1L))]
  M <- function(nm) do.call(rbind, lapply(ok, `[[`, nm))
  est <- list(`HAL-MLE` = M("psi_hal"), `T-HAL-MLE (Plugin)` = M("psi_plugin"), `T-HAL-MLE` = M("psi_hat"))
  SE <- M("se"); B <- nrow(est[[1]]); band <- 1 / log(n)
  for (nm in names(est)) {
    H <- est[[nm]]; e <- sweep(H, 2, tr); sdmc <- apply(H, 2, sd)
    rmse <- sqrt(colMeans(e^2)); bias <- colMeans(e); ratio <- bias / sdmc
    cov_mc <- colMeans(abs(e) <= 1.96 * rep(sdmc, each = B))
    wald <- if (nm == "T-HAL-MLE") colMeans(abs(e) <= 1.96 * SE, na.rm = TRUE) else rep(NA_real_, ncol(H))
    r <- q3(rmse); rows[[length(rows) + 1]] <- data.frame(
      escenario = sub("\\.rds$", "", basename(f)), estimador = nm, B = B,
      RMSE_Q1 = r[1], RMSE_mediana = r[2], RMSE_Q3 = r[3], RMSE_max = max(rmse),
      sesgo_abs_mediana = median(abs(bias)), prop_dosis_sesgo_despreciable = mean(abs(ratio) <= band),
      cobMCSE_mediana = median(cov_mc), cobMCSE_min = min(cov_mc),
      cobWald_mediana = if (nm == "T-HAL-MLE") median(wald, na.rm = TRUE) else NA, cobWald_min = if (nm == "T-HAL-MLE") min(wald, na.rm = TRUE) else NA)
  }
}
out <- do.call(rbind, rows); rownames(out) <- NULL
print(out, digits = 3, row.names = FALSE)
write.csv(out, file.path(outdir, "resumen_tipo_paper.csv"), row.names = FALSE)
cat("\nGuardado:", file.path(outdir, "resumen_tipo_paper.csv"), "\n")
