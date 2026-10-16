#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Contrasta si la varianza Monte Carlo de psi_hat sigue la pendiente de la verdad al cuadrado, psi0'(a)^2, como
# ocurriria si el estimador es una funcion escalonada o lineal a trozos con nudos que cambian de una muestra a otra
# (variabilidad horizontal x pendiente). Compara con el SE^2 medio del IC.
# Uso: Rscript experiments/diag_pendiente.R <outdir>
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "mc_out_cf"
files <- sort(list.files(outdir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) == 0) stop("No hay .rds en ", outdir)
for (f in files) {
  o <- readRDS(f); res <- if (!is.null(o$results)) o$results else o$res; tr <- o$truth; g <- o$grid
  ok <- Filter(function(x) isTRUE(x$ok), res)
  ok <- Filter(function(x) is.null(o$B) || x$b <= o$B, ok)
  ok <- ok[!duplicated(vapply(ok, function(x) x$b, 1L))]
  H <- do.call(rbind, lapply(ok, `[[`, "psi_hat")); SE <- do.call(rbind, lapply(ok, `[[`, "se"))
  m <- length(g); h <- g[2] - g[1]
  sl2 <- c(NA, ((tr[3:m] - tr[1:(m - 2)]) / (2 * h))^2, NA)        # psi0'(a)^2 (diferencia central)
  vmc <- apply(H, 2, var); se2 <- colMeans(SE^2, na.rm = TRUE); ii <- 2:(m - 1)
  cat("\n", sub("\\.rds$", "", basename(f)), " (B = ", nrow(H), ")\n", sep = "")
  cat(sprintf("Spearman con psi0'(a)^2 sobre %d dosis:  var_MC = %.2f   SE2_medio = %.2f\n",
              length(ii), cor(vmc[ii], sl2[ii], method = "spearman"), cor(se2[ii], sl2[ii], method = "spearman")))
  jj <- seq(1, m, by = 3)
  tab <- rbind(dosis = g[jj], pendiente2 = sl2[jj], var_MC = vmc[jj], SE2_medio = se2[jj])
  colnames(tab) <- NULL; print(round(tab, 2))
}
