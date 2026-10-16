#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Contrasta si la variabilidad por SELECCION de la base (que el IC trata como fija) explica el SE corto en ciertas
# dosis. Para cada dosis a0 se agrupan las replicas por el numero de nudos seleccionados en [a0 - w, a0 + w] y se
# compara el error cuadratico medio real (MSE) con el SE^2 medio del IC en cada grupo.
# Requiere .rds de mc_cf.R con sel y knots (base k = 1: columna 1 = constante, 2 = lineal, 2 + j = bisagra en el nudo j).
# Uso: Rscript experiments/diag_seleccion.R <outdir> [dosis_indices = 5,9,13,21,25] [ventana w = 0.4]
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "mc_out_cf"
idx <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else c(5, 9, 13, 21, 25)
w <- if (length(args) >= 3) as.numeric(args[3]) else 0.4
files <- sort(list.files(outdir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) == 0) stop("No hay .rds en ", outdir)
for (f in files) {
  o <- readRDS(f); res <- if (!is.null(o$results)) o$results else o$res; tr <- o$truth
  ok <- Filter(function(x) isTRUE(x$ok) && !is.null(x$sel) && !is.null(x$knots), res)
  ok <- Filter(function(x) is.null(o$B) || x$b <= o$B, ok)
  ok <- ok[!duplicated(vapply(ok, function(x) x$b, 1L))]
  if (length(ok) < 10) { cat("\n", basename(f), ": sin sel/knots en los .rds\n", sep = ""); next }
  H <- do.call(rbind, lapply(ok, `[[`, "psi_hat")); SE <- do.call(rbind, lapply(ok, `[[`, "se"))
  # posiciones de los nudos seleccionados por replica (solo bisagras: columnas >= 3)
  kpos <- lapply(ok, function(x) { j <- x$sel[x$sel >= 3] - 2; x$knots[j] })
  cat("\n", sub("\\.rds$", "", basename(f)), " (B = ", length(ok), ", ventana +-", w, ")\n", sep = "")
  for (j in idx[idx <= ncol(H)]) {
    a0 <- o$grid[j]
    nin <- vapply(kpos, function(k) sum(abs(k - a0) <= w), 1)
    grp <- cut(nin, c(-Inf, 0, 1, 2, Inf), labels = c("0", "1", "2", "3+"))
    e2 <- (H[, j] - tr[j])^2; s2 <- SE[, j]^2
    tab <- do.call(rbind, lapply(levels(grp), function(g) {
      i <- which(grp == g); if (length(i) == 0) return(NULL)
      data.frame(nudos_en_ventana = g, n = length(i), MSE = mean(e2[i]), SE2_medio = mean(s2[i], na.rm = TRUE),
                 MSE_sobre_SE2 = mean(e2[i]) / mean(s2[i], na.rm = TRUE))
    }))
    cat("\ndosis ", round(a0, 2), ": MSE total ", round(mean(e2), 3), ", SE2 medio ", round(mean(s2, na.rm = TRUE), 3), "\n", sep = "")
    print(tab, digits = 3, row.names = FALSE)
  }
}
cat("\nLectura: si MSE_sobre_SE2 > 1 y cambia mucho entre grupos, la seleccion de la base (no recogida en el IC) aporta varianza.\n")
