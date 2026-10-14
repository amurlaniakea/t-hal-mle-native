#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Diagnostico por dosis de los .rds de mc_paper.R: cobertura de Wald, cociente SE medio / sd Monte Carlo
# y sesgo / sd Monte Carlo (9 dosis equiespaciadas). Uso: Rscript experiments/diag_cobertura.R <outdir>
args <- commandArgs(trailingOnly = TRUE)
fs <- sort(list.files(args[1], pattern = "^dgp.*\\.rds$", full.names = TRUE))
idx <- c(1, 4, 7, 10, 13, 16, 19, 22, 25)
for (f in fs) {
  o <- readRDS(f); ok <- Filter(function(x) isTRUE(x$ok), o$res); tr <- o$truth
  M <- function(nm) do.call(rbind, lapply(ok, `[[`, nm))
  H <- M("psi_hat"); SE <- M("se")
  ii <- idx[idx <= ncol(H)]
  cov <- colMeans(abs(sweep(H, 2, tr)) <= 1.96 * SE, na.rm = TRUE)
  ratio <- colMeans(SE, na.rm = TRUE) / apply(H, 2, sd)
  bias_sd <- (colMeans(H) - tr) / apply(H, 2, sd)
  cat("\n", sub("\\.rds$", "", basename(f)), " (B =", nrow(H), ")\n", sep = "")
  print(round(rbind(dosis = o$grid[ii], cob_Wald = cov[ii], SE_medio_sobre_sd_MC = ratio[ii], sesgo_sobre_sd = bias_sd[ii]), 2))
}
