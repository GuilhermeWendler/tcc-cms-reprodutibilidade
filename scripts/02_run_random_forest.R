#!/usr/bin/env Rscript

# Reproduz o Random Forest do CMSclassifier para minPosterior 0.5, 0.4 e 0.3.
# A lógica de harmonização é a mesma dos scripts finais históricos.

required <- c("randomForest", "AnnotationDbi", "org.Hs.eg.db")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Pacotes R ausentes: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

dir.create("outputs", showWarnings = FALSE)

rfiles <- list.files("CMSclassifier/R", pattern = "\\.R$", full.names = TRUE)
if (length(rfiles) == 0) stop("Não foram encontrados arquivos R em CMSclassifier/R.")
for (f in rfiles) source(f)

load("CMSclassifier/data/model.rda")  # finalModel + mgenes

# Patch histórico necessário para evitar a falha de dispatch observada no projeto.
# Preserva center=TRUE, probabilidades, nearestCMS e minPosterior.
classifyCMS.RF.repro <- function(Exp, center = TRUE, minPosterior = 0.5) {
  if (center) Exp <- Exp - rowMeans(Exp)

  res <- as.data.frame(
    randomForest:::predict.randomForest(finalModel, t(Exp), type = "prob")
  )

  res$predictedCMS <- res$nearestCMS <- apply(res, 1, function(z) {
    y <- which(z == max(z))
    paste(c("CMS1", "CMS2", "CMS3", "CMS4")[y], collapse = ",")
  })

  res$predictedCMS[
    which(apply(res[, 1:4, drop = FALSE], 1, max) < minPosterior)
  ] <- NA

  res
}

expression_path <- "data/TCGACRC_expression-merged.tsv"
if (!file.exists(expression_path)) stop("Arquivo ausente: ", expression_path)

# Matriz original: genes x amostras.
Exp <- read.table(
  expression_path,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# Mantém a heurística do script histórico.
if (nrow(Exp) < ncol(Exp)) Exp <- t(Exp)

# IDs de amostras no mesmo formato usado pelo pipeline histórico do RF.
colnames(Exp) <- gsub("-", ".", toupper(colnames(Exp)))
rownames(Exp) <- toupper(rownames(Exp))

sym <- rownames(Exp)
sym <- sub("\\|.*$", "", sym)
sym <- sub("\\.[0-9]+$", "", sym)
sym <- toupper(sym)

map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(sym),
  keytype = "SYMBOL",
  columns = "ENTREZID"
)

map <- map[!is.na(map$ENTREZID), ]
map <- map[!duplicated(map$SYMBOL), ]  # primeira correspondência por símbolo

keep <- intersect(map$SYMBOL, sym)
Exp2 <- Exp[match(keep, sym), , drop = FALSE]
rownames(Exp2) <- keep

map2 <- map[match(rownames(Exp2), map$SYMBOL), ]
stopifnot(all(map2$SYMBOL == rownames(Exp2)))
rownames(Exp2) <- as.character(map2$ENTREZID)

# Agrega múltiplos símbolos que convergem para o mesmo Entrez ID pela média.
u <- unique(rownames(Exp2))
Exp_entrez <- sapply(u, function(g) {
  idx <- which(rownames(Exp2) == g)
  if (length(idx) == 1) return(Exp2[idx, ])
  colMeans(Exp2[idx, , drop = FALSE])
})
Exp_entrez <- t(Exp_entrez)
Exp_entrez <- as.matrix(Exp_entrez)
mode(Exp_entrez) <- "numeric"

mgenes_chr <- as.character(mgenes)
present <- intersect(mgenes_chr, rownames(Exp_entrez))
missing <- setdiff(mgenes_chr, present)

cat("Matriz original:", nrow(Exp), "genes x", ncol(Exp), "amostras\n")
cat("mgenes no modelo:", length(mgenes_chr), "\n")
cat("mgenes presentes:", length(present), "\n")
cat("mgenes faltantes:", length(missing), "\n")

Exp_model <- Exp_entrez[present, , drop = FALSE]

if (length(missing) > 0) {
  add <- matrix(
    0,
    nrow = length(missing),
    ncol = ncol(Exp_model),
    dimnames = list(missing, colnames(Exp_model))
  )
  Exp_model <- rbind(Exp_model, add)
}

Exp_model <- Exp_model[mgenes_chr, , drop = FALSE]

stopifnot(nrow(Exp_model) == length(mgenes_chr))
stopifnot(all(rownames(Exp_model) == mgenes_chr))
stopifnot(is.numeric(Exp_model[1, 1]))

# Registra o mapeamento efetivamente usado para aumentar a rastreabilidade.
write.csv(
  map2[, c("SYMBOL", "ENTREZID")],
  "outputs/symbol_to_entrez_used.csv",
  row.names = FALSE
)

run_one <- function(threshold, outfile) {
  res <- classifyCMS.RF.repro(
    Exp_model,
    center = TRUE,
    minPosterior = threshold
  )

  out <- res
  out$sample_id <- rownames(out)
  out <- out[, c("sample_id", setdiff(colnames(out), "sample_id"))]
  write.csv(out, outfile, row.names = FALSE)

  cat("\nminPosterior =", threshold, "\n")
  print(table(out$predictedCMS, useNA = "ifany"))

  invisible(out)
}

run_one(0.5, "outputs/cms_rf_pred_reproduzido.csv")
run_one(0.4, "outputs/cms_rf_pred_reproduzido_mp040.csv")
run_one(0.3, "outputs/cms_rf_pred_reproduzido_mp030.csv")

cat("\nRF concluído.\n")
