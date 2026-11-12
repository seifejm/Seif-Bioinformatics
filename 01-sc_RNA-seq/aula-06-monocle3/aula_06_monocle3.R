# ================================================================
# SEIF BIOINFORMATICS | AULA 06
# Pseudotempo e trajetórias celulares com Monocle 3
# PBMC: trajetórias independentes por grupo celular
# ================================================================
# Execução: abra este arquivo no RStudio a partir da pasta da aula.
# Entrada: data/pbmc_seurat_final.rds (objeto Seurat da Aula 2).
# ATENÇÃO: execute primeiro com ROOT_NODE <- NULL para inspecionar
# o grafo. Depois defina um nó biologicamente justificável e rode
# novamente para calcular pseudotempo, módulos e figuras.
# ================================================================

# 0. INSTALAÇÃO (executar apenas uma vez, se necessário) -----------
# install.packages(c("Seurat", "remotes", "dplyr", "ggplot2", "pheatmap"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SingleCellExperiment")
# remotes::install_github("cole-trapnell-lab/monocle3")
# remotes::install_github("satijalab/seurat-wrappers")

# 1. PACOTES E CONFIGURAÇÕES --------------------------------------
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(pheatmap)

set.seed(12345)

ARQUIVO_SEURAT <- file.path("data", "pbmc_seurat_final.rds")
PASTA_RESULTADOS <- file.path("results", "monocle3")
GRUPO_ESCOLHIDO <- "Monocytes" # Monocytes, B_cells, CD4_T ou CD8_NK
ROOT_NODE <- NULL              # Exemplo: "Y_1" (NÃO escolher sem olhar o grafo)
N_BINS <- 30
QVALUE_CUTOFF <- 0.05
MORANS_CUTOFF <- 0.01
MODULE_RESOLUTION <- 1e-3
N_TOP_GENES <- 5
CORES <- 4

# Evita salvar arquivos de grande porte no diretório principal.
dir.create(PASTA_RESULTADOS, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(ARQUIVO_SEURAT)) {
  stop("Arquivo não encontrado: ", ARQUIVO_SEURAT,
       ". Consulte o README e coloque o objeto na pasta data/.")
}
if (!GRUPO_ESCOLHIDO %in% c("Monocytes", "B_cells", "CD4_T", "CD8_NK")) {
  stop("GRUPO_ESCOLHIDO inválido. Consulte as opções no início do script.")
}

salvar_figura <- function(grafico, nome, largura = 9, altura = 6) {
  ggplot2::ggsave(
    filename = file.path(PASTA_RESULTADOS, nome),
    plot = grafico, width = largura, height = altura, dpi = 300
  )
}

# 2. LER O OBJETO E DEFINIR GRUPOS CELULARES -----------------------
pbmc <- readRDS(ARQUIVO_SEURAT)
stopifnot(inherits(pbmc, "Seurat"))
if (!"celltype" %in% colnames(pbmc@meta.data)) {
  stop("O objeto Seurat precisa conter a coluna 'celltype'.")
}
if (!"umap" %in% names(pbmc@reductions)) {
  stop("O objeto Seurat precisa conter uma redução UMAP chamada 'umap'.")
}
DefaultAssay(pbmc) <- "RNA"

print(table(pbmc$celltype, useNA = "ifany"))

# Altere este mapeamento se as anotações do seu objeto forem diferentes.
mapa_grupos <- list(
  Monocytes = c("Mono (CD14)", "Mono (FCGR3A)"),
  B_cells   = "B cells",
  CD4_T     = "CD4 T",
  CD8_NK    = c("CD8 T", "NK cells")
)

pbmc$trajectory_group <- NA_character_
for (grupo in names(mapa_grupos)) {
  idx <- !is.na(pbmc$celltype) & pbmc$celltype %in% mapa_grupos[[grupo]]
  pbmc$trajectory_group[idx] <- grupo
}
print(table(pbmc$trajectory_group, useNA = "ifany"))

p_grupos <- DimPlot(pbmc, reduction = "umap", group.by = "trajectory_group") +
  ggtitle("Grupos celulares para análise de pseudotempo")
salvar_figura(p_grupos, "01_grupos_celulares.png")

# 3. CONVERTER O GRUPO ESCOLHIDO PARA MONOCLE 3 --------------------
# Cada grupo tem sua própria trajetória. Evitamos impor uma única
# trajetória a linhagens biologicamente distintas.
celulas <- colnames(pbmc)[pbmc$trajectory_group %in% GRUPO_ESCOLHIDO]
if (length(celulas) < 30) {
  stop("Poucas células no grupo ", GRUPO_ESCOLHIDO,
       " (", length(celulas), "). Verifique as anotações.")
}

pbmc_grupo <- subset(pbmc, cells = celulas)
cds <- SeuratWrappers::as.cell_data_set(pbmc_grupo, assay = "RNA")
SingleCellExperiment::rowData(cds)$gene_short_name <- rownames(cds)

# as.cell_data_set() normalmente importa o UMAP existente do Seurat.
# Se não houver UMAP no cds, será necessário recalculá-lo antes.
if (!"UMAP" %in% SingleCellExperiment::reducedDimNames(cds)) {
  stop("Conversão sem redução UMAP. Confira SeuratWrappers e o objeto original.")
}

# 4. CLUSTERIZAÇÃO E GRAFO PRINCIPAL -------------------------------
cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP")
cds <- monocle3::learn_graph(cds, use_partition = TRUE)

# Os nós Y_1, Y_2 etc. podem mudar entre execuções/datasets.
p_nos <- monocle3::plot_cells(
  cds, color_cells_by = "celltype", show_trajectory_graph = TRUE,
  label_principal_points = TRUE, label_branch_points = TRUE,
  label_leaves = TRUE, label_cell_groups = FALSE
) + ggtitle(paste(GRUPO_ESCOLHIDO, "— nós do grafo"))
salvar_figura(p_nos, "02_nos_grafo.png")

p_particoes <- monocle3::plot_cells(
  cds, color_cells_by = "partition", show_trajectory_graph = FALSE,
  label_cell_groups = FALSE
) + ggtitle(paste(GRUPO_ESCOLHIDO, "— partições"))
salvar_figura(p_particoes, "03_particoes.png")

if (is.null(ROOT_NODE) || identical(ROOT_NODE, "")) {
  saveRDS(cds, file.path(PASTA_RESULTADOS, "cds_antes_da_raiz.rds"))
  message("Confira results/monocle3/02_nos_grafo.png e 03_particoes.png.")
  stop("Selecione ROOT_NODE (ex.: 'Y_1') após inspecionar o grafo e rode novamente.",
       call. = FALSE)
}

nos_validos <- igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name
if (!ROOT_NODE %in% nos_validos) {
  stop("ROOT_NODE = '", ROOT_NODE, "' não existe neste grafo. Nós válidos: ",
       paste(nos_validos, collapse = ", "))
}

# IMPORTANTE: o nó raiz deve representar um estado inicial com
# justificativa biológica; não escolher apenas por estar na ponta.
# Se existirem partições desconectadas sem raiz, algumas células
# terão pseudotempo infinito (serão removidas dos gráficos a seguir).
cds <- monocle3::order_cells(cds, root_pr_nodes = ROOT_NODE)

# 5. VISUALIZAR E EXTRAIR O PSEUDOTEMPO ----------------------------
p_pseudotempo <- monocle3::plot_cells(
  cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE,
  label_cell_groups = FALSE, label_leaves = FALSE,
  label_branch_points = FALSE
) + ggtitle(paste(GRUPO_ESCOLHIDO, "— pseudotempo"))
salvar_figura(p_pseudotempo, "04_pseudotempo_umap.png")

pseudo <- monocle3::pseudotime(cds)
pseudo_df <- data.frame(
  cell = names(pseudo), pseudotime = as.numeric(pseudo),
  celltype = as.character(SummarizedExperiment::colData(cds)[names(pseudo), "celltype"]),
  stringsAsFactors = FALSE
) |>
  dplyr::filter(is.finite(pseudotime))

if (nrow(pseudo_df) < 10) {
  stop("Menos de 10 células com pseudotempo finito. Revise raiz e partições.")
}

write.csv(pseudo_df, file.path(PASTA_RESULTADOS, "pseudotempo_celulas.csv"),
          row.names = FALSE)

p_violino <- ggplot(pseudo_df, aes(x = celltype, y = pseudotime)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  theme_classic() +
  labs(x = "Tipo celular", y = "Pseudotempo",
       title = paste("Pseudotempo —", GRUPO_ESCOLHIDO))
salvar_figura(p_violino, "05_pseudotempo_por_tipo_celular.png")

# 6. GENES ASSOCIADOS À ESTRUTURA DO GRAFO -------------------------
# graph_test() testa autocorrelação de expressão sobre o grafo por
# Moran's I; associação ao grafo NÃO demonstra regulação causal.
graph_res <- monocle3::graph_test(
  cds, neighbor_graph = "principal_graph", cores = CORES
)
graph_df <- as.data.frame(graph_res) |>
  tibble::rownames_to_column("gene_id")
write.csv(graph_df, file.path(PASTA_RESULTADOS, "graph_test_todos_genes.csv"),
          row.names = FALSE)

genes_pseudotime <- graph_df |>
  dplyr::filter(!is.na(q_value), !is.na(morans_I),
                q_value < QVALUE_CUTOFF, morans_I > MORANS_CUTOFF) |>
  dplyr::arrange(dplyr::desc(morans_I))
if (nrow(genes_pseudotime) < 2) {
  stop("Não há genes suficientes após os filtros de Moran's I e q-value.")
}
write.csv(genes_pseudotime,
          file.path(PASTA_RESULTADOS, "genes_associados_ao_grafo.csv"),
          row.names = FALSE)

# 7. MÓDULOS GÊNICOS COM MONOCLE 3 --------------------------------
# Módulos do Monocle 3 não são equivalentes aos módulos WGCNA.
# O pré-processamento nesta cópia prepara os componentes necessários
# à identificação de módulos; preservamos o cds original e seu grafo.
cds_modulos <- monocle3::preprocess_cds(cds, method = "PCA", num_dim = 50)
ids_genes <- intersect(genes_pseudotime$gene_id, rownames(cds_modulos))

gene_module_df <- monocle3::find_gene_modules(
  cds_modulos[ids_genes, ], resolution = MODULE_RESOLUTION
)
gene_module_df <- as.data.frame(gene_module_df)
gene_module_df$gene_name <- as.character(
  SummarizedExperiment::rowData(cds)[gene_module_df$id, "gene_short_name"]
)
write.csv(gene_module_df, file.path(PASTA_RESULTADOS, "modulos_genicos.csv"),
          row.names = FALSE)
print(table(gene_module_df$module))

# 8. EXPRESSÃO AGREGADA DOS MÓDULOS AO LONGO DO PSEUDOTEMPO -------
# Bins por quantis: aproximadamente o mesmo número de células por bin.
# São posições relativas na trajetória, não intervalos de tempo real.
n_bins_usados <- min(N_BINS, nrow(pseudo_df))
pseudo_df$pseudotime_bin <- paste0(
  "PT_", sprintf("%02d", dplyr::ntile(pseudo_df$pseudotime, n_bins_usados))
)
bin_levels <- sort(unique(pseudo_df$pseudotime_bin))
cell_group_df <- data.frame(
  cell = pseudo_df$cell,
  cell_group = pseudo_df$pseudotime_bin
)

matriz_modulos <- monocle3::aggregate_gene_expression(
  cds, gene_module_df, cell_group_df
)
matriz_modulos <- as.matrix(matriz_modulos)
matriz_modulos <- matriz_modulos[, intersect(bin_levels, colnames(matriz_modulos)),
                                 drop = FALSE]
rownames(matriz_modulos) <- paste0("Module_", rownames(matriz_modulos))
write.csv(matriz_modulos,
          file.path(PASTA_RESULTADOS, "matriz_modulos_pseudotempo.csv"))

# Escala z-score por linha, exclusivamente para visualizar padrões.
# Módulos sem variação são retirados do heatmap (não dos resultados).
linhas_variaveis <- apply(matriz_modulos, 1, sd, na.rm = TRUE) > 0
if (sum(linhas_variaveis) >= 2 && ncol(matriz_modulos) >= 2) {
  pheatmap::pheatmap(
    matriz_modulos[linhas_variaveis, , drop = FALSE],
    scale = "row", cluster_rows = TRUE, cluster_cols = FALSE,
    main = paste("Módulos ao longo do pseudotempo —", GRUPO_ESCOLHIDO),
    filename = file.path(PASTA_RESULTADOS, "06_heatmap_modulos.png"),
    width = 10, height = 7
  )
}

# 9. TOP GENES DE CADA MÓDULO -------------------------------------
# Ordenação por Moran's I dentro do módulo. Não são necessariamente
# os genes de maior expressão ou hub genes da rede WGCNA.
top_genes <- gene_module_df |>
  dplyr::left_join(genes_pseudotime, by = c("id" = "gene_id")) |>
  dplyr::filter(!is.na(morans_I)) |>
  dplyr::group_by(module) |>
  dplyr::slice_max(morans_I, n = N_TOP_GENES, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::arrange(module, dplyr::desc(morans_I))
top_genes$gene_name[is.na(top_genes$gene_name) |
                      top_genes$gene_name == ""] <-
  top_genes$id[is.na(top_genes$gene_name) | top_genes$gene_name == ""]
write.csv(top_genes, file.path(PASTA_RESULTADOS, "top_genes_por_modulo.csv"),
          row.names = FALSE)

# Matriz de expressão média normalizada por bin de pseudotempo.
# normalized_counts(log) é uma transformação de visualização;
# suas médias não representam contagens brutas nem efeito causal.
expr <- monocle3::normalized_counts(cds, norm_method = "log")
expr <- expr[top_genes$id, pseudo_df$cell, drop = FALSE]
matriz_genes <- vapply(bin_levels, function(bin) {
  celulas_bin <- pseudo_df$cell[pseudo_df$pseudotime_bin == bin]
  Matrix::rowMeans(expr[, celulas_bin, drop = FALSE])
}, numeric(nrow(expr)))
if (is.null(dim(matriz_genes))) {
  matriz_genes <- matrix(matriz_genes, nrow = nrow(expr),
                         dimnames = list(rownames(expr), bin_levels))
}
rownames(matriz_genes) <- make.unique(top_genes$gene_name)
colnames(matriz_genes) <- bin_levels

anotacao_linhas <- data.frame(Module = paste0("Module_", top_genes$module))
rownames(anotacao_linhas) <- rownames(matriz_genes)
tamanhos_modulos <- as.numeric(table(factor(top_genes$module,
                                           levels = unique(top_genes$module))))
gaps_row <- head(cumsum(tamanhos_modulos), -1)
write.csv(matriz_genes, file.path(PASTA_RESULTADOS, "matriz_top_genes.csv"))

linhas_genes_variaveis <- apply(matriz_genes, 1, sd, na.rm = TRUE) > 0
if (sum(linhas_genes_variaveis) >= 2 && ncol(matriz_genes) >= 2) {
  matriz_plot <- matriz_genes[linhas_genes_variaveis, , drop = FALSE]
  anotacao_plot <- anotacao_linhas[rownames(matriz_plot), , drop = FALSE]
  tamanhos_plot <- as.numeric(table(factor(anotacao_plot$Module,
                                          levels = unique(anotacao_plot$Module))))
  gaps_plot <- head(cumsum(tamanhos_plot), -1)
  pheatmap::pheatmap(
    matriz_plot, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
    annotation_row = anotacao_plot, gaps_row = gaps_plot,
    fontsize_row = 8, fontsize_col = 7,
    main = paste("Top genes por módulo —", GRUPO_ESCOLHIDO),
    filename = file.path(PASTA_RESULTADOS, "07_heatmap_top_genes.png"),
    width = 12, height = max(6, nrow(matriz_plot) * 0.23)
  )
}

# 10. EXEMPLOS DE GENES AO LONGO DA TRAJETÓRIA ---------------------
# Seleciona os 4 genes com maior Moran's I do primeiro módulo.
modulo_exemplo <- unique(top_genes$module)[1]
genes_exemplo <- head(top_genes$id[top_genes$module == modulo_exemplo], 4)
if (length(genes_exemplo) > 0) {
  p_genes <- monocle3::plot_genes_in_pseudotime(
    cds[genes_exemplo, ], min_expr = 0.5
  )
  salvar_figura(p_genes, "08_genes_no_pseudotempo.png", 11, 7)
}

# 11. SALVAR OBJETO E FINALIZAR -----------------------------------
saveRDS(cds,
        file.path(PASTA_RESULTADOS, paste0("monocle3_", GRUPO_ESCOLHIDO, ".rds")))
message("Análise concluída. Resultados em: ", PASTA_RESULTADOS)

# PIPELINE: Seurat -> grupo celular -> Monocle 3 -> grafo -> raiz ->
# pseudotempo -> Moran's I -> módulos Monocle 3 -> heatmaps.
