#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Comprobacion POST-HOC (sin nuevas simulaciones): suaviza el SE de cada replica a lo largo de la rejilla de dosis
# (spline de suavizado sobre la varianza SE^2) y recalcula la cobertura de Wald. Si sube, parte del defecto es ruido
# del SE entre dosis dentro de cada replica; si no sube, el defecto es sistematico (el SE medio esta mal en esa dosis).
# Uso: Rscript experiments/diag_se_suavizado.R <outdir> [grados de libertad del spline = 6]
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "mc_out"
dfs <- if (length(args) >= 2) as.numeric(args[2]) else 6
idx <- c(1, 5, 9, 13, 17, 21, 25)
files <- sort(list.files(outdir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) == 0) stop("No hay .rds en ", outdir)
out <- list()
for (f in files) {
  o <- readRDS(f); res <- if (!is.null(o$results)) o$results else o$res; tr <- o$truth
  ok <- Filter(function(x) isTRUE(x$ok), res)
  ok <- Filter(function(x) is.null(o$B) || x$b <= o$B, ok)
  ok <- ok[!duplicated(vapply(ok, function(x) x$b, 1L))]
  M <- function(nm) do.call(rbind, lapply(ok, `[[`, nm))
  H <- M("psi_hat"); SE <- M("se"); B <- nrow(H); g <- o$grid
  SS <- t(apply(SE, 1, function(s) {
    okk <- is.finite(s) & s > 0
    if (sum(okk) < 8) return(s)
    sqrt(pmax(predict(smooth.spline(g[okk], s[okk]^2, df = dfs), g)$y, 1e-12))   # se suaviza la VARIANZA (no log SE: evita sesgo a la baja)
  }))
  E <- abs(sweep(H, 2, tr)); sdmc <- apply(H, 2, sd)
  cov_o <- colMeans(E <= 1.96 * SE, na.rm = TRUE); cov_s <- colMeans(E <= 1.96 * SS, na.rm = TRUE)
  cat("\n", sub("\\.rds$", "", basename(f)), " (B = ", B, ", df = ", dfs, ")\n", sep = "")
  cat(sprintf("media sobre las dosis: cob_Wald = %.3f, cob_Wald_suavizado = %.3f, cob_SE_ideal = %.3f\n",
              mean(cov_o), mean(cov_s), mean(sweep(E, 2, 1.96 * sdmc, "<="))))
  tab <- rbind(dosis = g[idx], cob_Wald = cov_o[idx], cob_suavizado = cov_s[idx])
  colnames(tab) <- NULL; print(round(tab, 2))
}
