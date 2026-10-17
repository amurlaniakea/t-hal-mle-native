#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Compara T-HAL (psi_hat) y el plug-in HAL (psi_hal) DOSIS A DOSIS en los .rds de mc_paper.R / mc_cf.R: RMSE, sesgo y
# desviacion tipica de cada uno, y la diferencia de error cuadratico medio. Tambien resume el RMSE medio con y sin los
# extremos de la rejilla, solo como analisis de sensibilidad (no como metrica de resultado).
# Uso: Rscript experiments/diag_rmse_dosis.R <outdir> [dosis_indices = todas de 3 en 3: 1,4,...,25] [extremos_excluidos = 3]
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "mc_out_cf"
idx <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else seq(1, 25, by = 3)
nex <- if (length(args) >= 3) as.integer(args[3]) else 3
files <- sort(list.files(outdir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) == 0) stop("No hay .rds en ", outdir)
for (f in files) {
  o <- readRDS(f); res <- if (!is.null(o$results)) o$results else o$res; tr <- o$truth; g <- o$grid
  ok <- Filter(function(x) isTRUE(x$ok) && !is.null(x$psi_hal), res)
  ok <- Filter(function(x) is.null(o$B) || x$b <= o$B, ok)
  ok <- ok[!duplicated(vapply(ok, function(x) x$b, 1L))]
  T <- do.call(rbind, lapply(ok, `[[`, "psi_hat")); Hh <- do.call(rbind, lapply(ok, `[[`, "psi_hal")); m <- ncol(T)
  eT <- sweep(T, 2, tr); eH <- sweep(Hh, 2, tr)
  mse <- function(e) colMeans(e^2)
  rT <- sqrt(mse(eT)); rH <- sqrt(mse(eH))
  d <- mse(eT) - mse(eH)                                            # >0: T-HAL peor en esa dosis
  se_d <- apply(eT^2 - eH^2, 2, sd) / sqrt(nrow(T))
  ii <- idx[idx <= m]
  cat("\n", sub("\\.rds$", "", basename(f)), " (B = ", nrow(T), ")\n", sep = "")
  tab <- rbind(dosis = g[ii], RMSE_HAL = rH[ii], RMSE_THAL = rT[ii], sesgo_HAL = colMeans(eH)[ii], sesgo_THAL = colMeans(eT)[ii],
               sd_HAL = apply(Hh, 2, sd)[ii], sd_THAL = apply(T, 2, sd)[ii], dMSE = d[ii], z_dMSE = d[ii] / se_d[ii])
  colnames(tab) <- NULL; print(round(tab, 2))
  keep <- (nex + 1):(m - nex)
  cat(sprintf("RMSE medio por dosis, todas: HAL %.3f, T-HAL %.3f | sin las %d dosis de cada extremo: HAL %.3f, T-HAL %.3f\n",
              mean(rH), mean(rT), nex, mean(rH[keep]), mean(rT[keep])))
  cat(sprintf("reparto del exceso de MSE de T-HAL sobre HAL: dosis del extremo inferior (%d) %.0f %%, centrales %.0f %%, extremo superior (%d) %.0f %%\n",
              nex, 100 * sum(d[1:nex]) / sum(d), 100 * sum(d[keep]) / sum(d), nex, 100 * sum(d[(m - nex + 1):m]) / sum(d)))
}
