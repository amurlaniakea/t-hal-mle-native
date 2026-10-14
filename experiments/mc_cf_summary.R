#!/usr/bin/env Rscript
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Resumen PAREADO de los .rds de mc_cf.R: mismo ajuste, dos errores estandar (residuos en muestra vs
# cross-fitting K = 2). Uso: Rscript experiments/mc_cf_summary.R <outdir>
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "mc_out_cf"
files <- sort(list.files(outdir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) == 0) stop("No hay .rds en ", outdir)

summ <- function(f) {
  o <- readRDS(f); res <- if (!is.null(o$results)) o$results else o$res
  ok <- Filter(function(x) isTRUE(x$ok) && !is.null(x$se_insample), res)
  ok <- Filter(function(x) x$b <= o$B, ok)
  ok <- ok[!duplicated(vapply(ok, function(x) x$b, 1L))]
  B <- length(ok); if (B < 2) return(NULL)
  M <- function(nm) do.call(rbind, lapply(ok, `[[`, nm))
  H <- M("psi_hat"); S1 <- M("se_insample"); S2 <- M("se"); tr <- o$truth
  err <- abs(sweep(H, 2, tr)); sdmc <- apply(H, 2, sd)
  c1 <- err <= 1.96 * S1; c2 <- err <= 1.96 * S2
  d <- rowMeans(c2, na.rm = TRUE) - rowMeans(c1, na.rm = TRUE)      # diferencia pareada por replica
  data.frame(
    escenario = sub("\\.rds$", "", basename(f)), B_validas = B,
    cob_Wald_insample = mean(c1, na.rm = TRUE), cob_Wald_cf = mean(c2, na.rm = TRUE),
    dif_cob_cf = mean(d), EE_dif = sd(d) / sqrt(B),
    ancho_insample = mean(2 * 1.96 * S1, na.rm = TRUE), ancho_cf = mean(2 * 1.96 * S2, na.rm = TRUE),
    ancho_MCSE = mean(2 * 1.96 * sdmc),
    SEsd_insample = mean(colMeans(S1, na.rm = TRUE) / sdmc), SEsd_cf = mean(colMeans(S2, na.rm = TRUE) / sdmc),
    sesgo_sd_abs = mean(abs(colMeans(H) - tr) / sdmc),
    cob_MCSE = mean(err <= 1.96 * rep(sdmc, each = B)))
}
out <- do.call(rbind, lapply(files, summ))
print(out, digits = 3, row.names = FALSE)
write.csv(out, file.path(outdir, "summary_cf.csv"), row.names = FALSE)
cat("\nGuardado:", file.path(outdir, "summary_cf.csv"), "\n")
