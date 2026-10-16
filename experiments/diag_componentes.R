#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Descompone, por dosis, la varianza estimada del IC en su componente de residuos (Q-barra, R) y de variacion de W
# (Q_W, S), y la compara con la varianza Monte Carlo de psi_hat. Requiere .rds de mc_cf.R con var_R, var_S, cov_RS.
# Uso: Rscript experiments/diag_componentes.R <outdir> [dosis_indices, p. ej. 1,13,25]
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "mc_out_cf"
idx <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else c(1, 5, 9, 13, 17, 21, 25)
files <- sort(list.files(outdir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) == 0) stop("No hay .rds en ", outdir)
for (f in files) {
  o <- readRDS(f); res <- if (!is.null(o$results)) o$results else o$res
  ok <- Filter(function(x) isTRUE(x$ok) && !is.null(x$var_R), res)
  ok <- Filter(function(x) is.null(o$B) || x$b <= o$B, ok)
  ok <- ok[!duplicated(vapply(ok, function(x) x$b, 1L))]
  if (length(ok) < 2) { cat("\n", basename(f), ": sin componentes (var_R) en los .rds\n", sep = ""); next }
  M <- function(nm) do.call(rbind, lapply(ok, `[[`, nm))
  H <- M("psi_hat"); vR <- colMeans(M("var_R")); vS <- colMeans(M("var_S")); cRS <- colMeans(M("cov_RS"))
  vMC <- apply(H, 2, var); ii <- idx[idx <= ncol(H)]
  tab <- rbind(dosis = o$grid[ii], var_MC = vMC[ii],
               total_sobre_MC = (vR + vS + 2 * cRS)[ii] / vMC[ii],
               R_sobre_MC = vR[ii] / vMC[ii], S_sobre_MC = vS[ii] / vMC[ii], dosCov_sobre_MC = 2 * cRS[ii] / vMC[ii])
  cat("\n", sub("\\.rds$", "", basename(f)), " (B = ", nrow(H), ")\n", sep = "")
  colnames(tab) <- NULL; print(round(tab, 3))
}
cat("\nLectura: total_sobre_MC = (var_R + var_S + 2 cov_RS) / varianza Monte Carlo (1 = calibrado).",
    "\nSi R_sobre_MC domina y total_sobre_MC < 1, falla la componente de residuos; si S_sobre_MC domina, la de W.\n")
