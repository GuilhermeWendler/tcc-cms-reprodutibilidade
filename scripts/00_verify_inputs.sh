#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DATA="data/TCGACRC_expression-merged.tsv"
MODEL="CMSclassifier/data/model.rda"
CENTROIDS="CMSclassifier/data/centroids.RData"
DESC="CMSclassifier/DESCRIPTION"

for f in "$DATA" "$MODEL" "$CENTROIDS" "$DESC"; do
  if [[ ! -f "$f" ]]; then
    echo "ERRO: arquivo ausente: $f" >&2
    exit 1
  fi
done

echo "== Verificando entradas =="

DATA_MD5="$(md5sum "$DATA" | awk '{print $1}')"
MODEL_SHA="$(sha256sum "$MODEL" | awk '{print $1}')"
CENTROIDS_SHA="$(sha256sum "$CENTROIDS" | awk '{print $1}')"
CMS_VER="$(awk -F': ' '$1=="Version"{print $2}' "$DESC" | tr -d '\r')"

[[ "$DATA_MD5" == "48335693339d54eceef4b57274cf8b75" ]] || {
  echo "ERRO: MD5 da matriz não corresponde ao arquivo Synapse syn2325328.1." >&2
  echo "Obtido: $DATA_MD5" >&2
  exit 1
}

[[ "$MODEL_SHA" == "d0a43aea4f682719d10f4f80054f28134f6715af6b1f2b8a6a081af0b130b756" ]] || {
  echo "ERRO: SHA-256 de model.rda inesperado." >&2
  echo "Obtido: $MODEL_SHA" >&2
  exit 1
}

[[ "$CENTROIDS_SHA" == "b03b9663cbf23a989aed822c75a8fe26915b629f26c7d57a33fb90fe4ef50fb3" ]] || {
  echo "ERRO: SHA-256 de centroids.RData inesperado." >&2
  echo "Obtido: $CENTROIDS_SHA" >&2
  exit 1
}

[[ "$CMS_VER" == "1.0.0" ]] || {
  echo "ERRO: versão do CMSclassifier diferente de 1.0.0: $CMS_VER" >&2
  exit 1
}

echo "OK: matriz = Synapse syn2325328.1 (MD5 confirmado)"
echo "OK: CMSclassifier = 1.0.0"
echo "OK: model.rda e centroids.RData com hashes esperados"
