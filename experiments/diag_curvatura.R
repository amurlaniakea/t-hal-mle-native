#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Relaciona la varianza Monte Carlo de psi_hat y el SE^2 medio del IC con la curvatura de la verdad |psi0''(a)|
# (segunda diferencia sobre la rejilla). Si la varianza real sigue la curvatura y el SE^2 no, el IC no recoge la
# variabilidad asociada a aproximar zonas curvas con una base lineal a trozos de nudos/seleccion variables.
# Uso: Rscript experiments/diag_curvatura.R <outdir>
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
  curv <- c(NA, abs(diff(tr, differences = 2)) / h^2, NA)           # |psi0''| en los puntos interiores
  vmc <- apply(H, 2, var); se2 <- colMeans(SE^2, na.rm = TRUE)
  ii <- 2:(m - 1)
  cat("\n", sub("\\.rds$", "", basename(f)), " (B = ", nrow(H), ")\n", sep = "")
  cat(sprintf("correlacion de Spearman con |psi0''|:  var_MC = %.2f   SE2_medio = %.2f   cociente var_MC/SE2 = %.2f\n",
              cor(vmc[ii], curv[ii], method = "spearman"), cor(se2[ii], curv[ii], method = "spearman"),
              cor(vmc[ii] / se2[ii], curv[ii], method = "spearman")))
  jj <- seq(1, m, by = 3)
  tab <- rbind(dosis = g[jj], curvatura = curv[jj], var_MC = vmc[jj], SE2_medio = se2[jj], cociente = vmc[jj] / se2[jj])
  colnames(tab) <- NULL; print(round(tab, 2))
}
