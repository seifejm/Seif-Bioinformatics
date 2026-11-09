# ============================================================================
# SEIF BIOINFORMATICS | AULA 02
# Análise de scRNA-seq com Seurat: dos dados à anotação celular
# ============================================================================
# Material didático em português. Executar os blocos na ordem.
# Exemplo: dados PBMC no formato 10x Genomics (genes humanos).
# Consulte o README.md antes de iniciar.

# 1. PACOTES ---------------------------------------------------------------
# Instalar previamente, se necessário: install.packages(c("Seurat", "dplyr", "ggplot2"))
library(Seurat)
library(dplyr)
library(ggplot2)

# Parâmetros usados na aula: ajuste conforme seu conjunto de dados.
diretorio_10x <- "data/filtered_feature_bc_matrix"
n_genes_variaveis <- 2000
n_pcs_calculados <- 50
pcs_usados <- 1:10   # Exemplo didático: escolha após inspecionar o ElbowPlot.
resolucao <- 0.5

# 2. LEITURA DOS DADOS -----------------------------------------------------
# Read10X() lê os arquivos matrix.mtx(.gz), features.tsv(.gz) e
# barcodes.tsv(.gz). Linhas da matriz = genes; colunas = barcodes/células.
if (!dir.exists(diretorio_10x)) {
  stop("Pasta 10x não encontrada: ", diretorio_10x,
       ". Confira a estrutura indicada no README.md.")
}
dados <- Read10X(data.dir = diretorio_10x)
if (is.list(dados)) {
  # Alguns dados 10x incluem modalidades adicionais, como anticorpos.
  # Aqui usamos somente a matriz de expressão gênica.
  dados <- dados[["Gene Expression"]]
}
print(dim(dados))              # Número de genes x número de barcodes.
print(head(rownames(dados)))   # Nomes dos genes.
print(head(colnames(dados)))   # Identificadores das células.

# 3. OBJETO SEURAT ---------------------------------------------------------
# min.cells: retém genes detectados em pelo menos 3 células.
# min.features: retém células com pelo menos 200 genes detectados.
# NÃO é um limite máximo de genes por célula.
pbmc <- CreateSeuratObject(
  counts = dados,
  project = "PBMC_Aula02",
  min.cells = 3,
  min.features = 200
)
print(pbmc)
# RNA/counts = contagens brutas. Após as próximas etapas, outras camadas
# poderão conter dados normalizados (data) e escalonados (scale.data).
# meta.data contém métricas como nCount_RNA e nFeature_RNA.

# 4. CONTROLE DE QUALIDADE -------------------------------------------------
# Para genes humanos, o prefixo mitocondrial normalmente é "MT-".
# Para camundongos, costuma ser "mt-". Confirme os nomes no seu dataset.
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# nFeature_RNA: genes detectados em cada célula.
# nCount_RNA: soma das contagens em cada célula.
# percent.mt: proporção das contagens associadas a genes mitocondriais.
# Valores altos podem sugerir estresse/dano; valores extremos de nFeature
# e nCount podem sugerir doublets, mas nenhuma métrica comprova isso sozinha.

# Na aula original, o filtro refinado estava comentado. Mantemos opcional
# para NÃO eliminar células sem antes avaliar as distribuições de QC.
aplicar_filtro_qc <- FALSE
if (aplicar_filtro_qc) {
  # Limiares APENAS ILUSTRATIVOS, não universais. Ajuste após avaliar os dados.
  pbmc <- subset(pbmc, subset =
                   nFeature_RNA > 200 &
                   nFeature_RNA < 2500 &
                   percent.mt < 5)
  print(pbmc)
}

# 5. NORMALIZAÇÃO ----------------------------------------------------------
# Corrige diferenças de profundidade total de contagens entre células.
# O padrão LogNormalize escala cada célula e aplica log1p aos valores.
# No Seurat v5, NormalizeData() precisa ser executado explicitamente aqui.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize",
                      scale.factor = 10000)
print(Layers(pbmc[["RNA"]]))
print(GetAssayData(pbmc, assay = "RNA", layer = "data")[1:5, 1:5])

# 6. GENES ALTAMENTE VARIÁVEIS E ESCALONAMENTO -----------------------------
# Seleciona genes cuja variação entre células é informativa para o PCA.
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst",
                             nfeatures = n_genes_variaveis)
VariableFeaturePlot(pbmc)
print(head(VariableFeatures(pbmc)))
# ScaleData centraliza e padroniza, por gene, os valores usados no PCA.
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))

# 7. PCA -------------------------------------------------------------------
# Resume a variação dos genes selecionados em componentes principais (PCs).
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc),
               npcs = n_pcs_calculados)
ElbowPlot(pbmc, ndims = n_pcs_calculados)
DimPlot(pbmc, reduction = "pca")
# Confira o ElbowPlot e as contribuições biológicas antes de manter 1:10.
# A visualização nos primeiros 2 PCs não substitui a análise dos demais.

# 8. VIZINHANÇA E CLUSTERIZAÇÃO --------------------------------------------
# FindNeighbors: constrói grafos de vizinhança com os PCs selecionados.
# FindClusters: detecta comunidades no grafo; resolution altera a granularidade.
pbmc <- FindNeighbors(pbmc, dims = pcs_usados)
pbmc <- FindClusters(pbmc, resolution = resolucao)
print(table(Idents(pbmc)))

# 9. UMAP ------------------------------------------------------------------
# Projeta as células em duas dimensões PARA VISUALIZAR a estrutura dos dados.
# UMAP não é o algoritmo que atribui os clusters calculados acima.
pbmc <- RunUMAP(pbmc, dims = pcs_usados, n.neighbors = 30,
                min.dist = 0.3, spread = 1, seed.use = 123)
DimPlot(pbmc, reduction = "umap", label = TRUE)
# Em uma única amostra, este exemplo não executa integração/batch correction.
# Para múltiplas amostras, avalie integração (RPCA/CCA/Harmony etc.) separadamente.

# 10. MARCADORES DIFERENCIAIS ENTRE CLUSTERS -------------------------------
# Cada cluster é comparado às demais células (one versus rest).
# only.pos: mantém genes mais expressos no cluster de interesse.
# min.pct: fração mínima de células expressando o gene em pelo menos um grupo.
# logfc.threshold: limiar mínimo de efeito usado na seleção para teste.
# p_val_adj: p-valor ajustado; avalie junto de efeito e % de expressão.
Idents(pbmc) <- "seurat_clusters"
markers_clusters <- FindAllMarkers(
  pbmc, assay = "RNA", test.use = "wilcox",
  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25
)
print(head(markers_clusters))
if (nrow(markers_clusters) > 0) {
  top10 <- markers_clusters %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
    ungroup()
  print(top10)
}

# 11. MARCADORES CANÔNICOS E ANOTAÇÃO --------------------------------------
# Avalie o conjunto de marcadores: um gene isolado raramente define um tipo.
# A presença e a numeração dos clusters mudam com QC, PCs e resolution.
genes_canonicos <- c("CD3D", "CD4", "CD8A", "MS4A1",
                     "CD14", "FCGR3A", "GNLY", "NKG7")
genes_canonicos <- intersect(genes_canonicos, rownames(pbmc))
FeaturePlot(pbmc, features = genes_canonicos)
DotPlot(pbmc, features = genes_canonicos) + RotatedAxis()
print(levels(pbmc))

# MAPEAMENTO ILUSTRATIVO da gravação original (clusters 0 a 9).
# Não aplique automaticamente: confira DotPlot, FeaturePlot, FindAllMarkers
# e a identidade biológica de cada cluster no SEU resultado.
mapa_exemplo <- c(
  "0" = "Mono (CD14)", "1" = "CD4 T", "2" = "CD4 T",
  "3" = "CD8 T", "4" = "B cells", "5" = "Mono (FCGR3A)",
  "6" = "Mono (CD14)", "7" = "B cells", "8" = "NK cells",
  "9" = "Unknown"
)
# Só ative depois de validar TODAS as correspondências acima.
anotacao_validada <- FALSE
if (anotacao_validada) {
  if (!setequal(levels(pbmc), names(mapa_exemplo))) {
    stop("Os IDs de clusters diferem do exemplo. Atualize mapa_exemplo.")
  }
  pbmc <- RenameIdents(pbmc, mapa_exemplo)
  pbmc$celltype <- as.character(Idents(pbmc))
  DimPlot(pbmc, reduction = "umap", group.by = "celltype", label = TRUE)
} else {
  message("Anotação não aplicada: valide mapa_exemplo antes de ativá-la.")
}

# 12. SALVAR RESULTADOS ----------------------------------------------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)
write.csv(markers_clusters, "results/markers_clusters.csv", row.names = FALSE)
if (exists("top10")) {
  write.csv(top10, "results/top10_markers_clusters.csv", row.names = FALSE)
}
# Exportar marcadores por TIPO CELULAR somente quando a anotação estiver validada.
if (anotacao_validada) {
  Idents(pbmc) <- "celltype"
  markers_celltype <- FindAllMarkers(
    pbmc, assay = "RNA", test.use = "wilcox",
    only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25
  )
  write.csv(markers_celltype, "results/markers_celltypes.csv", row.names = FALSE)
}
saveRDS(pbmc, "results/pbmc_seurat_final.rds")
# Para reutilizar: pbmc <- readRDS("results/pbmc_seurat_final.rds")

# FIM: análise introdutória, não substitui QC/DE por desenho experimental.
