# Resultados esperados

A etapa `scripts/04_compare_methods.py` interrompe a execução caso os resultados
não reproduzam os valores finais preservados no TCC.

## Random Forest — minPosterior = 0.5

- total: 577
- CMS1: 66
- CMS2: 141
- CMS3: 85
- CMS4: 150
- NA: 135
- classificados: 442

## Concordância RF × centróides, sem limiar

- concordantes: 414/442
- concordância bruta: 93,6652% (93,67%)
- kappa de Cohen: 0,912621

Matriz de contagens:

| RF \ centróide | CMS1 | CMS2 | CMS3 | CMS4 |
|---|---:|---:|---:|---:|
| CMS1 | 64 | 0 | 2 | 0 |
| CMS2 | 0 | 132 | 5 | 4 |
| CMS3 | 3 | 7 | 75 | 0 |
| CMS4 | 4 | 3 | 0 | 143 |
| NA | 20 | 51 | 30 | 34 |

## Limiar 0.15

- UNCLASSIFIED no total: 52/577
- cobertura: 91,0%
- UNCLASSIFIED no subconjunto RF n=442: 10
- concordâncias usando denominador fixo 442: 407
- concordância: 92,08%

## Limiar 0.20

- UNCLASSIFIED no total: 94/577
- cobertura: 83,7%
- UNCLASSIFIED no subconjunto RF n=442: 29
- concordâncias usando denominador fixo 442: 390
- concordância: 88,24%

## Sensibilidade do RF

- minPosterior 0.4: 33 NA, 4 empates; 100% de estabilidade no n=442 original
- minPosterior 0.3: 2 NA, 6 empates; 100% de estabilidade no n=442 original
