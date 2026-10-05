# Dados de entrada

O arquivo grande de expressão gênica não é incluído neste repositório.

## Matriz TCGA COREAD

Arquivo esperado:

`data/TCGACRC_expression-merged.tsv`

Origem documentada:

- Synapse ID: `syn2325328`
- versão: `1`
- identificador versionado: `syn2325328.1`
- DOI: `10.7303/syn2325328.1`
- MD5 oficial/local confirmado: `48335693339d54eceef4b57274cf8b75`

A cópia utilizada na análise apresentou SHA-256:

`e9fc48ae1dcd2fe15e23eec6693919441bf207a175fcd1d67bc24de462088624`

## CMSclassifier

O pipeline espera uma árvore local:

`CMSclassifier/`

correspondente ao pacote CMSclassifier versão 1.0.0, contendo pelo menos:

- `CMSclassifier/R/*.R`
- `CMSclassifier/data/model.rda`
- `CMSclassifier/data/centroids.RData`
- `CMSclassifier/DESCRIPTION`

Os hashes esperados estão em `checksums/KNOWN_HASHES.txt`.
