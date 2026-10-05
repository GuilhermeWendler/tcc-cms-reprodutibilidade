#!/usr/bin/env python3

"""
Classificador por centróides inspirado no SSP.

IMPORTANTE:
- não é uma reprodução integral de classifyCMS.SSP;
- usa quatro centróides TCGA HiSeq (CMS1-CMS4);
- usa z-score e produto interno médio, preservando a lógica do script histórico;
- gera os cenários sem limiar, 0.15 e 0.20 em uma única execução.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
X_PATH = ROOT / "data" / "TCGACRC_expression-merged.tsv"
C_PATH = ROOT / "outputs" / "centroids_symbols_tcga_hiseq.csv"
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)


def norm_gene(idx):
    s = pd.Series(idx).astype(str).str.strip().str.upper()
    s = s.str.split("|", n=1).str[0]
    s = s.str.replace(r"\.\d+$", "", regex=True)
    return pd.Index(s)


if not X_PATH.exists():
    raise FileNotFoundError(X_PATH)
if not C_PATH.exists():
    raise FileNotFoundError(
        f"{C_PATH} não existe. Execute scripts/01_prepare_centroids.R primeiro."
    )

X_raw = pd.read_csv(X_PATH, sep="\t", index_col=0)
X = X_raw.T if X_raw.shape[0] > X_raw.shape[1] else X_raw
X.index = X.index.astype(str)

cent = pd.read_csv(C_PATH)
if "gene_id" in cent.columns:
    cent = cent.set_index("gene_id")

X.columns = norm_gene(X.columns)
cent.index = norm_gene(cent.index)

common_genes = X.columns.intersection(cent.index)
print("Genes em comum:", len(common_genes))

if len(common_genes) == 0:
    raise ValueError("0 genes em comum entre matriz de expressão e centróides.")

# Valor observado no ambiente documentado.
if len(common_genes) != 632:
    raise RuntimeError(
        f"Reprodução divergente: esperavam-se 632 genes em comum; obtidos {len(common_genes)}."
    )

Xg = X[common_genes].copy()
Cg = cent.loc[common_genes].copy()

# Preserva exatamente a transformação histórica.
Xz = (Xg - Xg.mean(axis=0)) / (Xg.std(axis=0, ddof=0) + 1e-8)
Cz = (Cg - Cg.mean(axis=0)) / (Cg.std(axis=0, ddof=0) + 1e-8)

scores = (Xz.values @ Cz.values) / Xz.shape[1]
scores_df = pd.DataFrame(scores, index=Xz.index, columns=Cz.columns)

pred = scores_df.idxmax(axis=1)
score_max = scores_df.max(axis=1)


def save_scenario(threshold, filename):
    pred_thr = pred.where(score_max >= threshold, other="UNCLASSIFIED")
    out = pd.DataFrame(
        {
            "sample_id": Xz.index,
            "CMS_SSP": pred.astype(str),          # nome histórico mantido por compatibilidade
            "score_max": score_max.values,
            "CMS_SSP_thr": pred_thr.astype(str), # nome histórico mantido por compatibilidade
            "threshold": threshold,
        }
    )
    out.to_csv(OUT / filename, index=False)
    print(
        f"{filename}: UNCLASSIFIED="
        f"{int((out['CMS_SSP_thr'] == 'UNCLASSIFIED').sum())}"
    )
    return out


save_scenario(0.0, "cms_ssp_pred.csv")
save_scenario(0.15, "cms_ssp_pred_thr015.csv")
save_scenario(0.20, "cms_ssp_pred_thr020.csv")

print("Classificador por centróides concluído.")
