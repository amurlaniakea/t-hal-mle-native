# validation/

Enfoque actual (Opción A, repo 100% nativo):
- La validación numérica se hace con tests unitarios puros en `tests/testthat/`.
- No se usan dependencias externas (hal9001, glmnet, haldensify) ni siquiera como oráculos.
- La fidelidad matemática se verifica contra van der Laan (2023), arXiv:2301.13354.
- Para el modelo de outcome Q_n se usan oráculos (funciones Q_fun) en los tests.

Si en el futuro se desea validar contra hal9001/glmnet, se pueden añadir scripts opcionales aquí,
pero no son necesarios para el funcionamiento del paquete.
