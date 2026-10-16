#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Separa "SE medio mal calibrado" de "SE ruidoso entre replicas" en los .rds de mc_paper.R / mc_cf.R.
# Por dosis: cobertura de Wald real, cobertura si se usara el SE MEDIO en todas las replicas (SE fijo),
# coeficiente de variacion del SE entre replicas y correlacion entre |error| y SE.
# Uso: Rscript experiments/diag_se_ruido.R <outdir> [dosis_indices, p. ej. 1,13,25]
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "mc_out"
idx <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else c(1, 4, 7, 10, 13, 16, 19, 22, 25)
files <- sort(list.files(outdir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) == 0) stop("No hay .rds en ", outdir)
for (f in files) {
  o <- readRDS(f); res <- if (!is.null(o$results)) o$results else o$res; tr <- o$truth
  ok <- Filter(function(x) isTRUE(x$ok), res)
  ok <- Filter(function(x) is.null(o$B) || x$b <= o$B, ok)
  ok <- ok[!duplicated(vapply(ok, function(x) x$b, 1L))]
  M <- function(nm) do.call(rbind, lapply(ok, `[[`, nm))
  H <- M("psi_hat"); SE <- M("se"); B <- nrow(H)
  ii <- idx[idx <= ncol(H)]
  tab <- sapply(ii, function(j) {
    e <- H[, j] - tr[j]; s <- SE[, j]
    c(dosis = o$grid[j],
      cob_Wald = mean(abs(e) <= 1.96 * s, na.rm = TRUE),
      cob_SE_fijo = mean(abs(e) <= 1.96 * mean(s, na.rm = TRUE), na.rm = TRUE),
      cob_SE_ideal = mean(abs(e) <= 1.96 * sd(H[, j])),
      CV_SE = sd(s, na.rm = TRUE) / mean(s, na.rm = TRUE),
      cor_abs_error_SE = suppressWarnings(cor(abs(e), s, use = "complete.obs")))
  })
  cat("\n", sub("\\.rds$", "", basename(f)), " (B = ", B, ")\n", sep = "")
  colnames(tab) <- NULL
  print(round(tab, 2))
}
cat("\nLectura: cob_SE_fijo (SE medio en todas las replicas) vs cob_Wald (SE de cada replica): si cob_SE_fijo es mucho mayor,",
    "\nel ruido del SE entre replicas (CV_SE alto) o su correlacion con el error cuesta cobertura.",
    "\ncob_SE_ideal usa la sd Monte Carlo (incluye el sesgo): es el techo con un SE perfecto.\n")
