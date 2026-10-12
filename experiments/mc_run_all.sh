#!/bin/sh
# Copyright (C) 2026 Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
# Lanza el piloto (n = 500) de los 3 DGP x 2 tratamientos. Reanudable.
# Uso: sh experiments/mc_run_all.sh <B> <outdir> <cores>
B=${1:-50}; OUT=${2:-mc_out}; CORES=${3:-1}; N=${4:-500}
mkdir -p "$OUT"
for DGP in 1 2 3; do
  for TRT in uniform normal; do
    Rscript experiments/mc_paper.R "$DGP" "$TRT" "$N" "$B" "$OUT" "$CORES" >> "$OUT/log_dgp${DGP}_${TRT}_n${N}.txt" 2>&1
  done
done
echo "FIN $(date '+%F %T')" >> "$OUT/FIN.txt"
