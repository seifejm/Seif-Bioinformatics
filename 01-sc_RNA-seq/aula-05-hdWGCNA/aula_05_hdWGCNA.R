# ============================================================
# SEIF BIOINFORMATICS | AULA 05
# hdWGCNA: redes de coexpressão em dados scRNA-seq (PBMC)
# ============================================================
# Este é um roteiro didático para execução sequencial no RStudio.
# Entrada: objeto Seurat da Aula 2, anotado em pbmc$celltype.
# Saída: resultados em results/wgcna/.
# Execute a partir da raiz do projeto RStudio.
# IMPORTANTE: resultados dependem dos dados e parâmetros; não são fixos.

# 0) INSTALAÇÃO (executar somente quando necessário) ------------
# install.packages(c("remotes", "dplyr", "ggplot2", "patchwork", "cowplot"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("WGCNA", "UCell", "GenomicRanges", "GeneOverlap"))
# remotes::install_github("smorabit/hdWGCNA")

# 1) PACOTES, CAMINHOS E OBJETO -------------------------------
library(Seurat)
library(hdWGCNA)
library(WGCNA)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

set.seed(12345)
theme_set(cowplot::theme_cowplot())
# Ajuste ao número de núcleos disponíveis; pode ser omitido em alguns sistemas.
WGCNA::enableWGCNAThreads(nThreads = 4)

arquivo_entrada <- "data/pbmc_seurat_final.rds"
pasta_resultados <- "results/wgcna"
pasta_tom <- file.path(pasta_resultados, "TOM")
dir.create(pasta_tom, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(arquivo_entrada)) {
  stop("Arquivo não encontrado: ", arquivo_entrada,
       "\nCopie o RDS da Aula 2 para a pasta data/.")
}
pbmc <- readRDS(arquivo_entrada)
if (!"celltype" %in% colnames(pbmc@meta.data)) {
  stop("A coluna 'celltype' não existe. Execute a anotação celular da Aula 2.")
}
# orig.ident distingue amostras na criação das metacélulas.
if (!"orig.ident" %in% colnames(pbmc@meta.data)) {
  stop("É necessária uma coluna orig.ident com os identificadores das amostras.")
}
if (!"pca" %in% Reductions(pbmc)) {
  stop("O objeto deve conter uma redução 'pca' calculada previamente.")
}
DefaultAssay(pbmc) <- "RNA"
Idents(pbmc) <- "celltype"
print(table(pbmc$celltype))

# 2) SELEÇÃO DE GENES -----------------------------------------
# fraction = 0.05: gene detectado em >= 5% das células.
# A seleção de genes altera os módulos obtidos.
pbmc <- SetupForWGCNA(
  pbmc, gene_select = "fraction", fraction = 0.05,
  wgcna_name = "PBMC_WGCNA"
)

# 3) METACÉLULAS ----------------------------------------------
# Metacélulas agregam células similares para reduzir a esparsidade.
# group.by evita misturar tipos celulares e amostras diferentes.
# k = vizinhos; min_cells = mínimo para construir o grupo.
pbmc <- MetacellsByGroups(
  seurat_obj = pbmc,
  group.by = c("celltype", "orig.ident"),
  ident.group = "celltype", reduction = "pca",
  k = 20, max_shared = 10, min_cells = 25,
  assay = "RNA", layer = "counts"
)
pbmc <- NormalizeMetacells(pbmc)
metacelulas <- GetMetacellObject(pbmc)
print(table(metacelulas$celltype))

# Visualização opcional: PCA + UMAP das metacélulas (não constrói a rede).
pbmc <- ScaleMetacells(pbmc, features = VariableFeatures(pbmc))
pbmc <- RunPCAMetacells(pbmc, features = VariableFeatures(pbmc))
pbmc <- RunUMAPMetacells(pbmc, reduction = "pca", dims = 1:10)
p_metacelulas <- DimPlotMetacells(pbmc, group.by = "celltype") +
  ggtitle("Metacélulas por tipo celular")
print(p_metacelulas)
ggsave(file.path(pasta_resultados, "01_umap_metacelulas.png"),
       p_metacelulas, width = 8, height = 6, dpi = 300)

# 4) MATRIZ DE EXPRESSÃO PARA WGCNA ---------------------------
# Linhas = metacélulas; colunas = genes.
# Usamos os tipos presentes nas metacélulas, em uma rede geral.
celltypes_wgcna <- sort(unique(as.character(metacelulas$celltype)))
celltypes_wgcna <- setdiff(celltypes_wgcna, NA_character_)
if (length(celltypes_wgcna) == 0L) stop("Nenhum tipo celular com metacélulas.")
pbmc <- SetDatExpr(
  pbmc, group_name = celltypes_wgcna, group.by = "celltype",
  assay = "RNA", layer = "data", use_metacells = TRUE
)

# 5) SOFT POWER -----------------------------------------------
# A potência controla quanto a rede enfatiza correlações mais fortes.
# O ajuste scale-free é um critério diagnóstico, não garantia biológica.
pbmc <- TestSoftPowers(pbmc, networkType = "signed")
p_power <- wrap_plots(PlotSoftPowers(pbmc), ncol = 2)
print(p_power)
ggsave(file.path(pasta_resultados, "02_soft_power.png"),
       p_power, width = 10, height = 7, dpi = 300)
power_table <- GetPowerTable(pbmc)
write.csv(power_table, file.path(pasta_resultados, "tabela_soft_power.csv"),
          row.names = FALSE)

# Seleção didática automática: menor potência >=3 com SFT.R.sq >= 0.80.
# Se nenhuma cumprir, selecionar o melhor ajuste e REVISAR o diagnóstico.
possiveis <- power_table %>%
  filter(!is.na(SFT.R.sq), Power >= 3, SFT.R.sq >= 0.80) %>%
  arrange(Power)
if (nrow(possiveis) > 0) {
  soft_power_use <- possiveis$Power[1]
} else {
  validos <- power_table %>% filter(!is.na(SFT.R.sq))
  if (nrow(validos) == 0) stop("Sem valores válidos de soft power.")
  soft_power_use <- validos$Power[which.max(validos$SFT.R.sq)]
  warning("Nenhum power atingiu 0.80; examine a tabela antes de interpretar.")
}
message("Soft power selecionado: ", soft_power_use)

# 6) REDE, TOM E MÓDULOS -------------------------------------
# TOM: sobreposição topológica baseada nas conexões dos genes.
# minModuleSize controla tamanho mínimo do módulo; mergeCutHeight,
# a fusão de módulos semelhantes. Ajuste conforme os dados.
pbmc <- ConstructNetwork(
  pbmc, soft_power = soft_power_use,
  tom_name = "PBMC_all", tom_outdir = pasta_tom,
  networkType = "signed", TOMType = "signed",
  minModuleSize = 30, mergeCutHeight = 0.50,
  overwrite_tom = FALSE
)
# O módulo grey representa genes não atribuídos a módulos definidos.
p_dendro <- PlotDendrogram(pbmc, main = "PBMC - módulos hdWGCNA")
print(p_dendro)
# Evitamos exportar toda a TOM como CSV: pode ser muito grande.
# Para inspecionar uma pequena parte: GetTOM(pbmc)[1:5, 1:5]

# 7) EIGENGENES E CONECTIVIDADE -------------------------------
# ME: primeiro componente principal da expressão do módulo.
# kME: correlação gene-ME, útil para classificar candidatos a hubs.
pbmc <- ModuleEigengenes(pbmc, assay = "RNA")
pbmc <- ModuleConnectivity(
  pbmc, group.by = "celltype", group_name = celltypes_wgcna,
  assay = "RNA", slot = "data", harmonized = FALSE
)
# Renomeia os módulos para um prefixo padronizado (exceto grey).
pbmc <- ResetModuleNames(pbmc, new_name = "PBMC-M")
modules <- GetModules(pbmc)
modules_no_grey <- modules %>% filter(as.character(module) != "grey")
hub_df <- GetHubGenes(pbmc, n_hubs = 10)
write.csv(modules_no_grey,
          file.path(pasta_resultados, "pbmc_hdWGCNA_modules.csv"),
          row.names = FALSE)
write.csv(hub_df,
          file.path(pasta_resultados, "pbmc_hdWGCNA_hub_genes.csv"),
          row.names = FALSE)
print(head(hub_df))

# 8) VISUALIZAÇÃO DA ATIVIDADE DOS MÓDULOS --------------------
# ModuleFeaturePlot: valores de ME projetados no UMAP celular.
p_me <- ModuleFeaturePlot(pbmc, features = "MEs", reduction = "umap",
                          order_points = TRUE)
print(wrap_plots(p_me, ncol = 3))

# Comparação entre tipos celulares usando os MEs, sem confundir com
# expressão percentual de genes. Cada ponto resume um ME médio.
MEs <- GetMEs(pbmc, harmonized = FALSE)
# Em geral, GetMEs retorna uma linha por célula; alinhamos por barcode.
if (!all(rownames(pbmc@meta.data) %in% rownames(MEs))) {
  warning("Índices de MEs e metadados não coincidem; ignorando gráfico por tipo.")
} else {
  me_cols <- colnames(MEs)
  me_cols <- me_cols[!grepl("grey", me_cols, ignore.case = TRUE)]
  if (length(me_cols) > 0L) {
    tabela_me <- cbind(
      celltype = as.character(pbmc$celltype),
      as.data.frame(MEs[rownames(pbmc@meta.data), me_cols, drop = FALSE])
    )
    resumo_me <- tabela_me %>%
      group_by(celltype) %>%
      summarise(across(all_of(me_cols), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
      tidyr::pivot_longer(-celltype, names_to = "modulo", values_to = "ME_medio")
    p_me_tipo <- ggplot(resumo_me, aes(modulo, celltype, fill = ME_medio)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "#3569B0", mid = "white", high = "#CE4A44") +
      labs(x = "Módulo", y = "Tipo celular", fill = "ME médio",
           title = "Módulos por tipo celular") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p_me_tipo)
    ggsave(file.path(pasta_resultados, "03_modulos_por_tipo.png"),
           p_me_tipo, width = 9, height = 5, dpi = 300)
  }
}

# Correlação entre os eigengenes dos módulos.
print(ModuleCorrelogram(pbmc))

# 9) REDES E HUB GENES (OPCIONAL, PODE SER DEMORADO) ----------
# UMAP abaixo representa GENES, não células. Não substitui o UMAP do Seurat.
# Execute os comandos opcionais apenas se o objeto tiver módulos válidos.
mods <- unique(as.character(modules_no_grey$module))
if (length(mods) > 0) {
  pbmc <- RunModuleUMAP(
    pbmc, n_hubs = 10, n_neighbors = 15,
    min_dist = 0.1, exclude_grey = TRUE
  )
  print(ModuleUMAPPlot(
    pbmc, edge.alpha = 0.25, sample_edges = TRUE,
    edge_prop = 0.10, label_hubs = 2, keep_grey_edges = FALSE
  ))

  # Rede de um módulo escolhido; ajuste ao módulo de interesse.
  module_use <- mods[1]
  print(hub_df %>% filter(as.character(module) == module_use))
  print(HubGeneNetworkPlot(
    pbmc, n_hubs = 10, n_other = 15,
    edge_prop = 0.75, mods = module_use
  ))
}
# ModuleNetworkPlot() pode gerar uma figura por módulo e custar tempo.
# Consulte a documentação da versão instalada antes de usá-lo.

# 10) SALVAR O OBJETO PROCESSADO -------------------------------
# O RDS pode ser grande; não enviar objetos volumosos ao GitHub.
saveRDS(pbmc, file.path(pasta_resultados, "pbmc_hdWGCNA_final.rds"))
message("Análise concluída. Saídas em: ", pasta_resultados)

# Fluxo: Seurat -> seleção genes -> metacélulas -> soft power ->
# rede/TOM -> módulos -> eigengenes -> kME -> hub genes -> figuras.
