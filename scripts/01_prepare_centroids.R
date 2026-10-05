#!/usr/bin/env Rscript

# Prepara os quatro centróides TCGA HiSeq usados na análise comparativa.
# Preserva a lógica do script histórico make_centroids_symbols_hiseq.R.

required <- c("AnnotationDbi", "org.Hs.eg.db")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Pacotes R ausentes: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

dir.create("outputs", showWarnings = FALSE)

CENTROIDS_RDATA <- "CMSclassifier/data/centroids.RData"
if (!file.exists(CENTROIDS_RDATA)) stop("Arquivo ausente: ", CENTROIDS_RDATA)

e <- new.env(parent = emptyenv())
load(CENTROIDS_RDATA, envir = e)

pick_centroids <- function(env) {
  candidates <- ls(env)[vapply(ls(env), function(o) {
    x <- get(o, envir = env)
    is.matrix(x) || is.data.frame(x)
  }, logical(1))]
  if (length(candidates) == 0) stop("Nenhuma matriz/data.frame encontrada em centroids.RData.")
  if (length(candidates) > 1) {
    message("Aviso: múltiplos objetos tabulares encontrados; usando o primeiro em ordem de ls(): ",
            candidates[1])
  }
  as.data.frame(get(candidates[1], envir = env))
}

Cfull <- pick_centroids(e)

if ("gene_id" %in% colnames(Cfull)) {
  Cfull$ENTREZID <- as.character(Cfull$gene_id)
} else {
  Cfull$ENTREZID <- as.character(rownames(Cfull))
}

cols <- c(
  "TCGA.HiSeq.centered.CMS1",
  "TCGA.HiSeq.centered.CMS2",
  "TCGA.HiSeq.centered.CMS3",
  "TCGA.HiSeq.centered.CMS4"
)

missing_cols <- setdiff(cols, colnames(Cfull))
if (length(missing_cols) > 0) {
  stop("Colunas de centróides ausentes: ", paste(missing_cols, collapse = ", "))
}

C <- Cfull[, c("ENTREZID", cols)]
for (cc in cols) C[[cc]] <- as.numeric(C[[cc]])

map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(C$ENTREZID),
  keytype = "ENTREZID",
  columns = "SYMBOL"
)

map <- map[!is.na(map$SYMBOL), ]
map <- map[!duplicated(map$ENTREZID), ]

C2 <- merge(map, C, by = "ENTREZID", all = FALSE)
Cagg <- aggregate(C2[, cols], by = list(gene_id = C2$SYMBOL), FUN = mean)
colnames(Cagg) <- c("gene_id", "CMS1", "CMS2", "CMS3", "CMS4")

write.csv(
  Cagg,
  "outputs/centroids_symbols_tcga_hiseq.csv",
  row.names = FALSE
)

cat("Centróides preparados.\n")
cat("Genes em símbolo:", nrow(Cagg), "\n")

if (nrow(Cagg) != 683) {
  warning("Esperavam-se 683 genes após o mapeamento no ambiente documentado; obtidos: ", nrow(Cagg))
}
