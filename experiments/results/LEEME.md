# Piloto Monte Carlo (n = 500, B = 50)

- Script: `experiments/mc_paper.R` (DGP exactos del paper, seccion 7.1), lanzado con
  `sh experiments/mc_run_all.sh 50 ~/mc_out 6 500`; semilla base 20261012.
- Codigo: paquete en `main` (52fcaca) + scripts de `feat/montecarlo` (b408b00).
- Resumen: `piloto_n500_B50.csv` (generado con `experiments/mc_summary.R`).
- Limitaciones: solo n = 500 y B = 50 (el paper usa B = 1000 y n = 500, 1000, 2000).
