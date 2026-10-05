# Reprodutibilidade da classificação CMS em câncer colorretal

Pipeline organizado para reproduzir a análise final do TCC: aplicação da floresta aleatória
pré-treinada do **CMSclassifier 1.0.0** e comparação com um **classificador por centróides
inspirado no SSP** na coorte TCGA COREAD.

## O que este repositório reproduz

1. preparação dos quatro centróides TCGA HiSeq;
2. harmonização SYMBOL → Entrez e aplicação do Random Forest;
3. análise de sensibilidade de `minPosterior` = 0.5, 0.4 e 0.3;
4. classificação por centróides sem limiar e com limiares 0.15 e 0.20;
5. matrizes de concordância, cobertura e kappa de Cohen;
6. verificações automáticas contra os resultados finais reportados.

**Nota metodológica:** o classificador por centróides usado aqui não é uma reprodução
integral de `classifyCMS.SSP`. O nome histórico `CMS_SSP` é mantido apenas em alguns
arquivos CSV para compatibilidade com as saídas antigas.

## Estrutura

```text
.
├── README.md
├── EXPECTED_RESULTS.md
├── run_all.sh
├── scripts/
│   ├── 00_verify_inputs.sh
│   ├── 01_prepare_centroids.R
│   ├── 02_run_random_forest.R
│   ├── 03_run_centroid_classifier.py
│   └── 04_compare_methods.py
├── environment/
│   ├── sessionInfo_analysis.txt
│   ├── python_environment.txt
│   └── requirements-minimal.txt
├── checksums/
│   └── KNOWN_HASHES.txt
├── data/
│   └── README.md
└── outputs/
```

## Entradas necessárias

### 1. Matriz de expressão

Coloque em:

`data/TCGACRC_expression-merged.tsv`

Origem: Synapse `syn2325328.1`, DOI `10.7303/syn2325328.1`.

O MD5 esperado é:

`48335693339d54eceef4b57274cf8b75`

### 2. CMSclassifier 1.0.0

Coloque a pasta `CMSclassifier/` na raiz do projeto. O pipeline utiliza:

- `CMSclassifier/data/model.rda`
- `CMSclassifier/data/centroids.RData`
- os arquivos R do pacote em `CMSclassifier/R/`

Antes da análise, `scripts/00_verify_inputs.sh` valida a versão e os hashes dos arquivos.

## Ambiente

R:

- R 4.3.3
- randomForest 4.7-1.2
- AnnotationDbi 1.64.1
- org.Hs.eg.db 3.18.0
- CMSclassifier 1.0.0

Python:

- Python 3.12.3
- NumPy 2.3.0
- pandas 2.3.0

## Execução

Com o ambiente Python correto ativo:

```bash
chmod +x run_all.sh scripts/00_verify_inputs.sh
./run_all.sh
```

Ou etapa por etapa:

```bash
bash scripts/00_verify_inputs.sh
Rscript --vanilla scripts/01_prepare_centroids.R
Rscript --vanilla scripts/02_run_random_forest.R
python scripts/03_run_centroid_classifier.py
python scripts/04_compare_methods.py
```

A última etapa imprime `=== REPRODUÇÃO OK ===` somente se os resultados coincidirem
com os valores finais do TCC.

## Principais saídas

- `outputs/cms_rf_pred_reproduzido.csv`
- `outputs/cms_rf_pred_reproduzido_mp040.csv`
- `outputs/cms_rf_pred_reproduzido_mp030.csv`
- `outputs/cms_ssp_pred.csv`
- `outputs/cms_ssp_pred_thr015.csv`
- `outputs/cms_ssp_pred_thr020.csv`
- `outputs/concordancia_RF_oficial_vs_SSP_*`
- `outputs/analysis_summary.csv`

Consulte `EXPECTED_RESULTS.md` para os valores esperados.

## Proveniência

A proveniência da matriz, do modelo e dos centróides foi verificada criptograficamente.
Os hashes conhecidos estão em `checksums/KNOWN_HASHES.txt`.
