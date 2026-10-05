#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

PYTHON_BIN="${PYTHON_BIN:-python}"

echo "[0/4] Verificando entradas"
bash scripts/00_verify_inputs.sh

echo "[1/4] Preparando centróides"
Rscript --vanilla scripts/01_prepare_centroids.R

echo "[2/4] Rodando Random Forest"
Rscript --vanilla scripts/02_run_random_forest.R

echo "[3/4] Rodando classificador por centróides"
"$PYTHON_BIN" scripts/03_run_centroid_classifier.py

echo "[4/4] Comparando métodos e validando resultados"
"$PYTHON_BIN" scripts/04_compare_methods.py

echo
echo "Pipeline finalizado com sucesso."
