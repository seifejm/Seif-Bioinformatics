# ============================================================================
# SEIF BIOINFORMATICS | AULA 03
# GENE ONTOLOGY (GO): ENRIQUECIMENTO FUNCIONAL EM DADOS DE scRNA-seq
# Exemplo didático com células humanas (PBMC) e clusterProfiler
# ============================================================================
# Objetivo: identificar termos GO super-representados nos genes marcadores
# de um tipo celular previamente identificado com Seurat.
#
# ENTRADAS (geradas na Aula 2):
#   data/pbmc_seurat_final.rds  -> objeto Seurat processado
#   data/markers_celltypes.csv -> tabela de marcadores por tipo celular
#
# SAÍDAS:
#   results/                   -> tabelas CSV, figuras PNG e objetos RDS
#
# ATENÇÃO: enriquecimento de GO NÃO demonstra, por si só, que uma via esteja
# ativada ou inibida. Ele identifica termos mais representados do que o
# esperado para um universo de genes definido.
# ============================================================================

# 1) INSTALAÇÃO (execute apenas uma vez, se necessário) -----------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot",
#                        "AnnotationDbi", "GO.db"))
# install.packages(c("ggplot2", "dplyr"))

# 2) PACOTES E PASTAS ---------------------------------------------------------
library(Seurat)          # leitura do objeto gerado na aula anterior
library(clusterProfiler) # teste de enriquecimento de GO
library(org.Hs.eg.db)   # anotações de genes humanos
library(enrichplot)     # gráficos de enriquecimento
library(AnnotationDbi) # consulta às bases de anotação
library(GO.db)          # consulta opcional a termos e definições do GO
library(dplyr)         # manipulação da tabela de marcadores
library(ggplot2)       # aparência dos gráficos e exportação

# Execute este script a partir da raiz do projeto RStudio.
arquivo_pbmc    <- "data/pbmc_seurat_final.rds"
arquivo_markers <- "data/markers_celltypes.csv"
pasta_saida    <- "results"
dir.create(pasta_saida, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(arquivo_pbmc) || !file.exists(arquivo_markers)) {
  stop("Arquivos de entrada não encontrados. Confira a pasta data/ e o README.")
}

# 3) CARREGAR OS RESULTADOS DA AULA 2 -----------------------------------------
pbmc <- readRDS(arquivo_pbmc)
markers_celltype <- read.csv(arquivo_markers, check.names = FALSE)

# A tabela de FindAllMarkers() deve conter estas colunas.
colunas_necessarias <- c("cluster", "gene", "p_val_adj", "avg_log2FC")
if (!all(colunas_necessarias %in% names(markers_celltype))) {
  stop("A tabela de marcadores não contém: ",
       paste(setdiff(colunas_necessarias, names(markers_celltype)), collapse = ", "))
}

print(unique(markers_celltype$cluster))  # tipos celulares disponíveis

# 4) ESCOLHER TIPO CELULAR E FILTRAR GENES ------------------------------------
celltype_alvo <- "B cells"  # altere conforme os nomes mostrados acima
padj_max      <- 0.05       # significância após correção múltipla
log2fc_min    <- 0.25       # expressão aumentada em relação ao grupo referência

if (!celltype_alvo %in% markers_celltype$cluster) {
  stop("Tipo celular não encontrado na tabela: ", celltype_alvo)
}

# Selecionamos apenas genes com log2FC positivo e p-valor ajustado significativo.
# A análise é de genes marcadores; não compara tratamentos experimentais.
genes_selecionados <- markers_celltype %>%
  filter(cluster == celltype_alvo,
         !is.na(p_val_adj), p_val_adj < padj_max,
         !is.na(avg_log2FC), avg_log2FC > log2fc_min) %>%
  pull(gene) %>%
  unique()

cat("Genes marcadores filtrados:", length(genes_selecionados), "\n")

# 5) DEFINIR O UNIVERSO / BACKGROUND ------------------------------------------
# O universo deve representar os genes que PODERIAM ter sido selecionados
# no experimento. Aqui usamos os genes presentes no assay RNA da Aula 2
# como aproximação didática. Em análises formais, prefira os genes que foram
# efetivamente testados na comparação diferencial.
DefaultAssay(pbmc) <- "RNA"
universo_experimento <- rownames(pbmc[["RNA"]])

# A base org.Hs.eg.db usa símbolos de genes humanos (ex.: CD3D, MS4A1).
# Intersectamos os símbolos do experimento com os reconhecidos na anotação.
simbolos_humanos <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
universo <- intersect(universo_experimento, simbolos_humanos)
genes_go <- intersect(genes_selecionados, universo)

cat("Genes no universo anotado:", length(universo), "\n")
cat("Genes selecionados presentes no universo:", length(genes_go), "\n")
if (length(genes_go) < 5L) {
  stop("Poucos genes disponíveis para GO. Confira os filtros e identificadores.")
}

# 6) EXECUTAR ENRIQUECIMENTO NAS TRÊS ONTOLOGIAS -----------------------------
# BP = Biological Process (processo biológico)
# MF = Molecular Function (função molecular)
# CC = Cellular Component (componente celular)
#
# O teste avalia super-representação dos genes selecionados dentro de cada
# termo GO, em relação ao universo. BH controla a taxa de falsas descobertas.
# O qvalueCutoff é um filtro adicional; ajuste-o conforme o objetivo.
rodar_go <- function(ontologia) {
  clusterProfiler::enrichGO(
    gene          = genes_go,
    universe      = universo,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ontologia,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    minGSSize     = 10,
    maxGSSize     = 500,
    readable      = TRUE
  )
}

resultados_go <- setNames(lapply(c("BP", "MF", "CC"), rodar_go),
                          c("BP", "MF", "CC"))

# 7) TABELAS E SALVAMENTO -----------------------------------------------------
# GeneRatio = genes selecionados anotados no termo / genes selecionados
#             reconhecidos no universo de anotação utilizado pela análise.
# BgRatio   = genes do universo anotados no termo / genes do universo
#             reconhecidos no universo de anotação utilizado pela análise.
# Count     = número de genes selecionados associados ao termo.
# p.adjust  = p-valor após correção de múltiplos testes (BH).
# geneID    = genes da lista associados ao termo.
# IMPORTANTE: as razões são relativas ao universo efetivamente analisado.

tabelas <- lapply(resultados_go, as.data.frame)
for (ont in names(tabelas)) {
  tab <- tabelas[[ont]]
  cat(ont, "- termos enriquecidos:", nrow(tab), "\n")
  write.csv(tab, file.path(pasta_saida, paste0("GO_", ont, ".csv")),
            row.names = FALSE)
}

tabela_completa <- dplyr::bind_rows(tabelas, .id = "Ontology")
write.csv(tabela_completa,
          file.path(pasta_saida, "GO_todas_ontologias.csv"), row.names = FALSE)
saveRDS(resultados_go, file.path(pasta_saida, "GO_objetos.rds"))

# Mostra resultados mais relevantes (se existirem).
if (nrow(tabela_completa) > 0L) {
  print(head(tabela_completa[, intersect(c("Ontology", "ID", "Description",
                                         "GeneRatio", "BgRatio", "p.adjust", "Count"),
                                       names(tabela_completa))], 10))
} else {
  message("Nenhum termo passou pelos filtros de enriquecimento.")
}

# 8) GRÁFICOS DE ENRIQUECIMENTO -----------------------------------------------
# Dotplot: normalmente combina significância, proporção de genes e contagem.
# Barplot: mostra as categorias enriquecidas em barras.
# Evitamos plotar resultados vazios.

salvar_grafico <- function(obj, ont, tipo = c("dotplot", "barplot")) {
  tipo <- match.arg(tipo)
  if (nrow(as.data.frame(obj)) == 0L) return(invisible(NULL))

  p <- if (tipo == "dotplot") {
    enrichplot::dotplot(obj, showCategory = 15)
  } else {
    graphics::barplot(obj, showCategory = 15)
  }

  p <- p + ggplot2::labs(title = paste("Gene Ontology -", ont)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   axis.text.y = ggplot2::element_text(size = 9))

  print(p)
  ggplot2::ggsave(filename = file.path(pasta_saida,
                                       paste0("GO_", ont, "_", tipo, ".png")),
                  plot = p, width = 9, height = 6, dpi = 300)
  invisible(p)
}

for (ont in names(resultados_go)) {
  salvar_grafico(resultados_go[[ont]], ont, "dotplot")
  salvar_grafico(resultados_go[[ont]], ont, "barplot")
}

# 9) REDES DE ENRIQUECIMENTO (OPCIONAIS) --------------------------------------
# Redes são úteis para explorar relações entre genes e termos, mas podem
# ficar difíceis de interpretar quando há muitas categorias.
# Execute se houver resultados em BP.
if (nrow(tabelas[["BP"]]) > 0L) {
  p_rede_genes <- enrichplot::cnetplot(resultados_go[["BP"]], showCategory = 3)
  print(p_rede_genes)
  ggplot2::ggsave(file.path(pasta_saida, "GO_BP_cnetplot.png"),
                  plot = p_rede_genes, width = 10, height = 7, dpi = 300)
}

if (nrow(tabelas[["BP"]]) >= 2L) {
  go_bp_sim <- enrichplot::pairwise_termsim(resultados_go[["BP"]])
  p_rede_termos <- enrichplot::emapplot(go_bp_sim, showCategory = 15)
  print(p_rede_termos)
  ggplot2::ggsave(file.path(pasta_saida, "GO_BP_emapplot.png"),
                  plot = p_rede_termos, width = 10, height = 7, dpi = 300)
}

# 10) EXPLORAR O BANCO GENE ONTOLOGY (CONSULTA OPCIONAL) ----------------------
# Essas consultas são ilustrativas e não são necessárias para rodar enrichGO.
# Mostra identificador, nome, tipo de ontologia e definição de um termo.
termo_exemplo <- AnnotationDbi::select(
  GO.db,
  keys = "GO:0006955", keytype = "GOID",
  columns = c("TERM", "ONTOLOGY", "DEFINITION")
)
print(termo_exemplo)

# Exemplo: consultar anotações de genes humanos no org.Hs.eg.db.
exemplo_genes <- intersect(c("CD3D", "CD4", "MS4A1", "CD14"), simbolos_humanos)
if (length(exemplo_genes) > 0L) {
  print(head(AnnotationDbi::select(
    org.Hs.eg.db, keys = exemplo_genes, keytype = "SYMBOL",
    columns = c("ENTREZID", "GENENAME", "GO", "ONTOLOGY", "EVIDENCE")
  ), 15))
}

# Registrar versões para reprodução e interpretação futura.
writeLines(capture.output(sessionInfo()),
           con = file.path(pasta_saida, "sessionInfo.txt"))

# FIM DA AULA 3 ---------------------------------------------------------------
