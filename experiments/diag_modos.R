#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Estructura del error de psi_hat entre dosis (sin nuevas simulaciones): correlaciones entre dosis clave, componentes
# principales de la matriz de errores (replicas x dosis) y colas del error en cada dosis. Sirve para ver si el exceso de
# varianza en una dosis viene de un unico "modo" compartido (p. ej. un desplazamiento o cambio de pendiente) o de ruido local.
# Uso: Rscript experiments/diag_modos.R <outdir> [dosis_indices = 5,9,13,17,21]
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "mc_out_cf"
idx <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else c(5, 9, 13, 17, 21)
files <- sort(list.files(outdir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) == 0) stop("No hay .rds en ", outdir)
for (f in files) {
  o <- readRDS(f); res <- if (!is.null(o$results)) o$results else o$res; tr <- o$truth
  ok <- Filter(function(x) isTRUE(x$ok), res)
  ok <- Filter(function(x) is.null(o$B) || x$b <= o$B, ok)
  ok <- ok[!duplicated(vapply(ok, function(x) x$b, 1L))]
  E <- sweep(do.call(rbind, lapply(ok, `[[`, "psi_hat")), 2, tr); B <- nrow(E); ii <- idx[idx <= ncol(E)]
  pc <- prcomp(E, center = TRUE, scale. = FALSE)
  sh <- pc$sdev^2 / sum(pc$sdev^2)
  cat("\n", sub("\\.rds$", "", basename(f)), " (B = ", B, ")\n", sep = "")
  cat("varianza explicada por las 3 primeras componentes principales:", paste0(round(100 * sh[1:3], 1), "%", collapse = " / "), "\n")
  cat("cargas de la CP1 en las dosis elegidas (signo y tamano relativo):\n")
  print(round(rbind(dosis = o$grid[ii], carga_CP1 = pc$rotation[ii, 1] / max(abs(pc$rotation[, 1])),
                    sd_error = apply(E[, ii, drop = FALSE], 2, sd)), 2))
  cat("correlacion de los errores entre dosis elegidas:\n")
  C <- cor(E[, ii, drop = FALSE]); dimnames(C) <- list(round(o$grid[ii], 2), round(o$grid[ii], 2)); print(round(C, 2))
  cat("replicas con |error| > 3 sd robustas por dosis:",
      paste(round(o$grid[ii], 2), vapply(ii, function(j) sum(abs(E[, j] - median(E[, j])) > 3 * IQR(E[, j]) / 1.349), 1), sep = ":", collapse = "  "), "\n")
}
