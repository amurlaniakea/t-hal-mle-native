# Piloto Monte Carlo (n = 500, B = 50)

- Script: `experiments/mc_paper.R` (DGP exactos del paper, seccion 7.1), lanzado con
  `sh experiments/mc_run_all.sh 50 ~/mc_out 6 500`; semilla base 20261012.
- Codigo: paquete en `main` (52fcaca) + scripts de `feat/montecarlo` (b408b00).
- Resumen: `piloto_n500_B50.csv` (generado con `experiments/mc_summary.R`).
- Limitaciones: solo n = 500 y B = 50 (el paper usa B = 1000 y n = 500, 1000, 2000).

# Topes de peso (experiments/mc_wcap.R)

- `wcap_exploratorio_n500_B50.csv`: DGP 1 y 2, n = 500, B = 50, tratamiento normal, semilla base 20261012;
  topes: paper, `obs_max`, 10, 5; nuisances estimadas y oraculo. Exploratorio (topes comparados en la misma simulacion).
- `wcap_confirmacion_n500_B60_s777001.csv`: DGP 1-3, n = 500, B = 60, tratamiento normal, semilla base 777001;
  topes: paper, 3, 5, 8, `obs_max`; solo nuisances estimadas. Lanzado con
  `Rscript experiments/mc_wcap.R <dgp> 500 60 ~/mc_out_conf 6 777001 "defecto,3,5,8,obs_max" est`.
- Resumenes generados con `experiments/mc_wcap_summary.R`.

# Cobertura con topes bajos

- `cobertura_topes_n500_B60_s888001.csv`: DGP 1-3, n = 500, B = 60, tratamiento normal, semilla base 888001;
  topes 3 y 2. Lanzado con `Rscript experiments/mc_paper.R <dgp> normal 500 60 ~/mc_out_cov 6 888001 <tope>`
  y resumido con `experiments/mc_summary.R`.
- `experiments/diag_cobertura.R <outdir>`: cobertura de Wald, SE medio / sd Monte Carlo y sesgo / sd por dosis.
- Nota: el RMSE de `mc_summary.R` es la media de los RMSE por dosis; el de `mc_wcap_summary.R` es el RMSE
  agregado sobre dosis y replicas. No son comparables en nivel (las diferencias frente al HAL si).

# Cross-fitting K = 2 para la varianza (pareado)

- `cf_pareado_n500_B60_s222001.csv`: DGP 1-3, n = 500, B = 60, tratamiento normal, tope 3, semilla base 222001;
  nuisances `est` (g y omega estimadas) y `orac_gw` (g y omega oraculo, Q ajustada). Lanzado con
  `Rscript experiments/mc_cf.R <dgp> normal 500 60 ~/mc_out_cfp 6 222001 3 <est|orac_gw>` y resumido con
  `experiments/mc_cf_summary.R` (comparacion pareada de los dos SE del mismo ajuste).
- Aviso: los resultados del DGP 3 frente al plug-in HAL anteriores al commit `f4d86a5` (piloto, topes y cobertura)
  usan un `Q_fun` mal escalado para `k_Q = 1` y no son validos; los DGP 1 y 2 no se ven afectados.

# Tratamiento uniforme, topes bajos del DGP 3 y tamano del sieve

- `uniforme_cf_n500_B60_s333001.csv` / `uniforme_rmse_n500_B60_s333001.csv`: DGP 1-3, tratamiento uniforme, ventana por
  defecto. `mc_cf_summary.R` y `mc_summary.R` sobre `mc_cf.R <dgp> uniform 500 60 ~/mc_out_cfu 6 333001 3 est`.
- `dgp3_topes_bajos_cf_n500_B60_s444001.csv` / `dgp3_topes_bajos_rmse_n500_B60_s444001.csv`: DGP 3, tratamiento normal,
  topes 2, 1.5 y 1 (`mc_cf.R 3 normal 500 60 ~/mc_out_d3caps 6 444001 <tope> est`).
- `sieve_cf_n500_B60_s555001_c9-14.csv` / `sieve_rmse_n500_B60_s555001_c9-14.csv`: DGP 1-3, uniforme, ventana
  `c1 = 9, c2 = 14` (`mc_cf.R <dgp> uniform 500 60 ~/mc_out_sieve 6 555001 3 est 9 14`).
- El RMSE de los `*_rmse_*.csv` es la media de los RMSE por dosis (como `piloto_n500_B50.csv`).

# Resolucion del HAL de Q (n_knots_Q y k_Q)

- `nkq20_*_dgp1_uniforme_n500_B30_s888001.csv`: DGP 1 uniforme, ventana 9-14, `n_knots_Q = 20`
  (`mc_cf.R 1 uniform 500 30 ~/mc_out_nkq 6 888001 3 est 9 14 20`).
- `kq1_*_dgp1_uniforme_n500_B60_s999001.csv`: DGP 1 uniforme, ventana 9-14, `k_Q = 1` con 10 nudos
  (`mc_cf.R 1 uniform 500 60 ~/mc_out_kq1 6 999001 3 est 9 14 10 1`).
- `nkq20_*_dgp123_n500_B30_s123001.csv`: DGP 3 normal, DGP 2 uniforme (ventana 9-14) y DGP 1 normal, `n_knots_Q = 20`
  (`mc_cf.R <dgp> <trt> 500 30 ~/mc_out_nkq2 6 123001 3 est <c1> <c2> 20`).
- `*_cf_*.csv` es la salida de `mc_cf_summary.R` y `*_rmse_*.csv` la de `mc_summary.R`.
- Diagnosticos: `experiments/diag_pendiente.R` (relacion de la varianza con la pendiente^2),
  `diag_modos.R`, `diag_seleccion.R`, `diag_componentes.R`, `diag_se_ruido.R`.

