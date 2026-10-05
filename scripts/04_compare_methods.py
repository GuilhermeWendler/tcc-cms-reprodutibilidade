#!/usr/bin/env python3

"""
Compara o RF do CMSclassifier com o classificador por centróides.

Reproduz:
- matrizes de concordância;
- concordância no subconjunto RF-classificado (n=442);
- kappa de Cohen no cenário sem limiar;
- cobertura dos limiares 0.15 e 0.20;
- sensibilidade do RF a minPosterior 0.5, 0.4 e 0.3;
- verificações estritas contra os resultados reportados no TCC.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"


def norm_sample(s):
    return (
        s.astype(str)
        .str.upper()
        .str.replace(".", "-", regex=False)
        .str.strip()
    )


def read_rf(filename):
    df = pd.read_csv(OUT / filename)
    df["sample_key"] = norm_sample(df["sample_id"])
    return df


def read_cent(filename):
    df = pd.read_csv(OUT / filename)
    df["sample_key"] = norm_sample(df["sample_id"])
    return df


rf05 = read_rf("cms_rf_pred_reproduzido.csv")
rf04 = read_rf("cms_rf_pred_reproduzido_mp040.csv")
rf03 = read_rf("cms_rf_pred_reproduzido_mp030.csv")

c0 = read_cent("cms_ssp_pred.csv")
c15 = read_cent("cms_ssp_pred_thr015.csv")
c20 = read_cent("cms_ssp_pred_thr020.csv")

if rf05["sample_key"].duplicated().any():
    raise RuntimeError("IDs duplicados no RF.")
if c0["sample_key"].duplicated().any():
    raise RuntimeError("IDs duplicados no classificador por centróides.")


def merge_pair(rf, cent):
    m = rf.merge(
        cent,
        on="sample_key",
        how="inner",
        suffixes=("_RF", "_CENT")
    )
    if len(m) != 577:
        raise RuntimeError(f"Esperavam-se 577 amostras após merge; obtidas {len(m)}.")
    return m


m0 = merge_pair(rf05, c0)
m15 = merge_pair(rf05, c15)
m20 = merge_pair(rf05, c20)

classes = ["CMS1", "CMS2", "CMS3", "CMS4"]
rows = classes + ["NA"]


def rf_table_label(series):
    return series.fillna("NA").astype(str)


def counts_table(m, cent_col, cols):
    temp = pd.DataFrame(
        {
            "CMS_RF": rf_table_label(m["predictedCMS"]),
            "CMS_CENT": m[cent_col].astype(str),
        }
    )
    tab = pd.crosstab(temp["CMS_RF"], temp["CMS_CENT"])
    return tab.reindex(index=rows, columns=cols, fill_value=0)


tab0 = counts_table(m0, "CMS_SSP", classes)
tab15 = counts_table(m15, "CMS_SSP_thr", classes + ["UNCLASSIFIED"])
tab20 = counts_table(m20, "CMS_SSP_thr", classes + ["UNCLASSIFIED"])

tab0.to_csv(OUT / "concordancia_RF_oficial_vs_SSP_noThr_counts.csv")
tab15.to_csv(OUT / "concordancia_RF_oficial_vs_SSP_thr015_counts.csv")
tab20.to_csv(OUT / "concordancia_RF_oficial_vs_SSP_thr020_counts.csv")

for name, tab in [
    ("noThr", tab0),
    ("thr015", tab15),
    ("thr020", tab20),
]:
    rowpct = tab.div(tab.sum(axis=1).replace(0, np.nan), axis=0).round(4)
    rowpct.to_csv(OUT / f"concordancia_RF_oficial_vs_SSP_{name}_rowpct.csv")


def classified_subset(m):
    return m[m["predictedCMS"].notna()].copy()


s0 = classified_subset(m0)
s15 = classified_subset(m15)
s20 = classified_subset(m20)

n_rf = len(s0)
agree0 = int((s0["predictedCMS"] == s0["CMS_SSP"]).sum())
agree15 = int((s15["predictedCMS"] == s15["CMS_SSP_thr"]).sum())
agree20 = int((s20["predictedCMS"] == s20["CMS_SSP_thr"]).sum())

rate0 = agree0 / n_rf
rate15 = agree15 / n_rf
rate20 = agree20 / n_rf

unclass15_all = int((m15["CMS_SSP_thr"] == "UNCLASSIFIED").sum())
unclass20_all = int((m20["CMS_SSP_thr"] == "UNCLASSIFIED").sum())
unclass15_rf = int((s15["CMS_SSP_thr"] == "UNCLASSIFIED").sum())
unclass20_rf = int((s20["CMS_SSP_thr"] == "UNCLASSIFIED").sum())


def cohen_kappa(y1, y2, labels):
    cm = pd.crosstab(y1, y2).reindex(
        index=labels, columns=labels, fill_value=0
    ).to_numpy(dtype=float)

    n = cm.sum()
    po = np.trace(cm) / n
    pe = (cm.sum(axis=1) * cm.sum(axis=0)).sum() / (n * n)
    return (po - pe) / (1.0 - pe)


kappa = cohen_kappa(s0["predictedCMS"], s0["CMS_SSP"], classes)

# Resumos compatíveis com os arquivos históricos.
def write_summary(path, rate, unclass_rate=None):
    lines = [
        f"n_merge={len(m0)}",
        f"rf_na_rate={rf05['predictedCMS'].isna().mean():.4f}",
        f"n_rf_classified={n_rf}",
        f"concordancia_global_rfClassified={rate:.4f}",
    ]
    if unclass_rate is not None:
        lines.append(f"ssp_unclassified_rate_on_rfClassified={unclass_rate:.4f}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


write_summary(
    OUT / "concordancia_RF_oficial_vs_SSP_noThr_summary.txt",
    rate0
)
write_summary(
    OUT / "concordancia_RF_oficial_vs_SSP_thr015_summary.txt",
    rate15,
    unclass15_rf / n_rf
)
write_summary(
    OUT / "concordancia_RF_oficial_vs_SSP_thr020_summary.txt",
    rate20,
    unclass20_rf / n_rf
)

# Sensibilidade do RF.

def rf_stats(df):
    na_n = int(df["predictedCMS"].isna().sum())

    ties_predicted = int(
        df["predictedCMS"]
        .fillna("")
        .astype(str)
        .str.contains(",", regex=False)
        .sum()
    )

    return na_n, ties_predicted

na05, ties05 = rf_stats(rf05)
na04, ties04 = rf_stats(rf04)
na03, ties03 = rf_stats(rf03)

base = rf05.loc[
    rf05["predictedCMS"].notna(),
    ["sample_key", "predictedCMS"]
].copy()

def stability(other):
    x = base.merge(
        other[["sample_key", "predictedCMS"]],
        on="sample_key",
        suffixes=("_base", "_other")
    )
    return float(
        (x["predictedCMS_base"] == x["predictedCMS_other"]).mean()
    )

stab04 = stability(rf04)
stab03 = stability(rf03)

def stability(other):
    x = base.merge(
        other[["sample_key", "predictedCMS"]],
        on="sample_key",
        suffixes=("_base", "_other")
    )
    return float((x["predictedCMS_base"] == x["predictedCMS_other"]).mean())


stab04 = stability(rf04)
stab03 = stability(rf03)

summary = pd.DataFrame(
    [
        ["n_total", 577],
        ["rf_classified_mp050", n_rf],
        ["rf_na_mp050", na05],
        ["agreement_no_threshold_n", agree0],
        ["agreement_no_threshold", rate0],
        ["cohen_kappa_no_threshold", kappa],
        ["centroid_unclassified_thr015_all", unclass15_all],
        ["centroid_coverage_thr015", 1 - unclass15_all / 577],
        ["agreement_thr015_fixed_n442", rate15],
        ["centroid_unclassified_thr015_on_rf442", unclass15_rf],
        ["centroid_unclassified_thr020_all", unclass20_all],
        ["centroid_coverage_thr020", 1 - unclass20_all / 577],
        ["agreement_thr020_fixed_n442", rate20],
        ["centroid_unclassified_thr020_on_rf442", unclass20_rf],
        ["rf_na_mp040", na04],
        ["rf_ties_mp040", ties04],
        ["rf_stability_original442_mp040", stab04],
        ["rf_na_mp030", na03],
        ["rf_ties_mp030", ties03],
        ["rf_stability_original442_mp030", stab03],
    ],
    columns=["metric", "value"]
)
summary.to_csv(OUT / "analysis_summary.csv", index=False)

# Matrizes esperadas, derivadas dos resultados finais preservados no TCC.
expected0 = np.array([
    [64, 0, 2, 0],
    [0, 132, 5, 4],
    [3, 7, 75, 0],
    [4, 3, 0, 143],
    [20, 51, 30, 34],
])

expected15 = np.array([
    [64, 0, 2, 0, 0],
    [0, 126, 3, 3, 9],
    [3, 7, 75, 0, 0],
    [4, 3, 0, 142, 1],
    [16, 34, 22, 21, 42],
])

expected20 = np.array([
    [64, 0, 2, 0, 0],
    [0, 119, 3, 2, 17],
    [3, 7, 74, 0, 1],
    [4, 2, 0, 133, 11],
    [14, 26, 14, 16, 65],
])

checks = [
    (np.array_equal(tab0.to_numpy(), expected0), "matriz sem limiar"),
    (np.array_equal(tab15.to_numpy(), expected15), "matriz thr=0.15"),
    (np.array_equal(tab20.to_numpy(), expected20), "matriz thr=0.20"),
    (n_rf == 442, "n RF classificado = 442"),
    (na05 == 135, "NA RF 0.5 = 135"),
    (agree0 == 414, "414 concordâncias sem limiar"),
    (abs(rate0 - 414/442) < 1e-12, "concordância 93.67%"),
    (abs(kappa - 0.9126210850348782) < 1e-12, "kappa = 0.912621..."),
    (unclass15_all == 52, "52 UNCLASSIFIED em 0.15"),
    (unclass20_all == 94, "94 UNCLASSIFIED em 0.20"),
    (unclass15_rf == 10, "10 UNCLASSIFIED no n=442 em 0.15"),
    (unclass20_rf == 29, "29 UNCLASSIFIED no n=442 em 0.20"),
    (agree15 == 407, "407 concordâncias em 0.15"),
    (agree20 == 390, "390 concordâncias em 0.20"),
    (na04 == 33, "NA RF 0.4 = 33"),
    (na03 == 2, "NA RF 0.3 = 2"),
    (ties04 == 4, "4 empates em 0.4"),
    (ties03 == 6, "6 empates em 0.3"),
    (stab04 == 1.0, "estabilidade 100% do n=442 em 0.4"),
    (stab03 == 1.0, "estabilidade 100% do n=442 em 0.3"),
]

failed = [name for ok, name in checks if not ok]
if failed:
    print("\nREPRODUÇÃO DIVERGENTE:")
    for item in failed:
        print(" -", item)
    raise SystemExit(1)

print("\n=== REPRODUÇÃO OK ===")
print(f"Amostras: {len(m0)}")
print(f"RF classificado (minPosterior=0.5): {n_rf}")
print(f"Concordância sem limiar: {agree0}/{n_rf} = {rate0:.4%}")
print(f"Kappa de Cohen: {kappa:.6f}")
print(f"thr=0.15: UNCLASSIFIED={unclass15_all}/577; concordância fixa={rate15:.4%}")
print(f"thr=0.20: UNCLASSIFIED={unclass20_all}/577; concordância fixa={rate20:.4%}")
print(f"RF 0.4: NA={na04}, empates={ties04}, estabilidade={stab04:.1%}")
print(f"RF 0.3: NA={na03}, empates={ties03}, estabilidade={stab03:.1%}")
