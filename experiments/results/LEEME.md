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
