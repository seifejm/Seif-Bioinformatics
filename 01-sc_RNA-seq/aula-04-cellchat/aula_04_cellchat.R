# ================================================================
# SEIF BIOINFORMATICS | AULA 04
# COMUNICAÇÃO CÉLULA-CÉLULA COM CELLCHAT (PBMC HUMANAS)
# ================================================================
# Objetivo: inferir interações ligante-receptor entre tipos celulares,
# analisar vias de sinalização e visualizar redes de comunicação.
# Este script NÃO inclui LIANA.
#
# Pré-requisito: objeto Seurat da Aula 2, anotado em pbmc$celltype.
# Execute as seções na ordem. Os resultados dependem dos seus dados.
# ================================================================

# 1) PACOTES -------------------------------------------------------
# Instale previamente, quando necessário:
# install.packages(c("Seurat", "dplyr", "ggplot2", "remotes"))
# remotes::install_github("jinworks/CellChat")
# Atenção: consulte a documentação da versão instalada do CellChat.

library(Seurat)
library(CellChat)
library(dplyr)
library(ggplot2)

# 2) DADOS DE ENTRADA ----------------------------------------------
# Ajuste o caminho caso o objeto tenha sido salvo em outra pasta.
arquivo_pbmc <- "data/pbmc_seurat_final.rds"
if (!file.exists(arquivo_pbmc)) {
  stop("Arquivo não encontrado: ", arquivo_pbmc,
       ". Copie o RDS da Aula 2 para a pasta data/.")
}
pbmc <- readRDS(arquivo_pbmc)

# 'celltype' deve conter a anotação dos tipos celulares feita na Aula 2.
if (!"celltype" %in% colnames(pbmc[[]])) {
  stop("A coluna 'celltype' não existe nos metadados do objeto Seurat.")
}
if (!"RNA" %in% Assays(pbmc)) {
  stop("O objeto precisa conter o assay RNA.")
}
DefaultAssay(pbmc) <- "RNA"
Idents(pbmc) <- "celltype"
print(table(Idents(pbmc)))  # Número de células por tipo celular

# CellChat utiliza expressão normalizada, não valores de scale.data.
# Execute NormalizeData na Aula 2 antes de salvar o objeto.
# Se houver múltiplas layers de amostras (Seurat v5), harmonize/junte
# as layers antes de criar o objeto CellChat, conforme sua estratégia.

# Pasta para salvar resultados, sem poluir a pasta principal.
outdir <- "results/cellchat"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 3) CRIAR OBJETO CELLCHAT E EXPLORAR O BANCO ----------------------
cellchat <- createCellChat(
  object = pbmc,
  group.by = "celltype",
  assay = "RNA"
)

# Banco humano de interações ligante-receptor.
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# As três categorias principais de comunicação no banco:
# Secreted Signaling = moléculas secretadas;
# Cell-Cell Contact  = contato direto entre células;
# ECM-Receptor       = matriz extracelular e receptores.
table(CellChatDB$interaction$annotation)

# Inspecione uma categoria e as vias que ela contém.
categoria_exemplo <- "Secreted Signaling"
interacoes_banco <- CellChatDB$interaction %>%
  filter(annotation == categoria_exemplo)
head(interacoes_banco)
sort(unique(interacoes_banco$pathway_name))
# Se desejar explorar todas as interações: View(CellChatDB$interaction)

# 4) INFERÊNCIA DE COMUNICAÇÃO -------------------------------------
# Seleciona genes de sinalização presentes na expressão do experimento.
cellchat <- subsetData(cellchat)

# Identifica genes e pares ligante-receptor com expressão relevante.
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Estima a probabilidade de comunicação entre grupos celulares.
# O cálculo pode demorar: utiliza testes/permutação de acordo com o método.
cellchat <- computeCommunProb(cellchat)

# Descarta interações de grupos com menos de 10 células.
# Ajuste min.cells de acordo com os tamanhos dos grupos e o desenho.
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Agrega pares ligante-receptor em vias de sinalização.
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Centralidade: identifica papéis relativos de emissão/recepção na rede.
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Salve após a etapa custosa para não precisar recalcular tudo.
saveRDS(cellchat, file.path(outdir, "pbmc_cellchat.rds"))

# 5) EXTRAIR E SALVAR RESULTADOS -----------------------------------
# net: informação no nível de interações ligante-receptor.
# netP: comunicação agregada por vias de sinalização.
print(names(cellchat@net))
print(names(cellchat@netP))

interacoes <- subsetCommunication(cellchat)
vias_resultado <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(interacoes, file.path(outdir, "interacoes_ligante_receptor.csv"), row.names = FALSE)
write.csv(vias_resultado, file.path(outdir, "comunicacao_por_via.csv"), row.names = FALSE)

# Examine source (emissor), target (receptor), ligand, receptor,
# prob (força/probabilidade inferida), pval e annotation.
head(interacoes)

# Exemplo de filtro por categoria de sinalização.
interacoes_secreted <- interacoes %>% filter(annotation == "Secreted Signaling")
head(interacoes_secreted)

# 6) REDE GLOBAL: NÚMERO E FORÇA DAS INTERAÇÕES --------------------
# O tamanho dos nós representa o número de células em cada grupo.
# count = quantidade de interações inferidas entre grupos;
# weight = soma das forças/probabilidades inferidas pelo método.
groupSize <- as.numeric(table(cellchat@idents))

netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Número de interações"
)

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Força das interações"
)

# 7) INVESTIGAR UMA VIA ESPECÍFICA ---------------------------------
# As vias disponíveis dependem dos seus dados: NÃO presuma uma via.
vias_disponiveis <- cellchat@netP$pathways
print(vias_disponiveis)

if (length(vias_disponiveis) > 0) {
  # Para a aula, começa pela primeira via; altere após examinar a lista.
  pathway.use <- vias_disponiveis[1]
  message("Via selecionada: ", pathway.use)

  # Compare o que existe no banco com o que foi inferido no experimento.
  componentes_banco <- CellChatDB$interaction %>%
    filter(pathway_name == pathway.use)
  interacoes_via <- subsetCommunication(cellchat, signaling = pathway.use)
  head(componentes_banco)
  head(interacoes_via)

  # Quais pares de ligante/receptor mais contribuem para essa via?
  pares_LR <- extractEnrichedLR(
    cellchat, signaling = pathway.use, geneLR.return = FALSE
  )
  print(pares_LR)
  netAnalysis_contribution(cellchat, signaling = pathway.use)

  # Gráfico circular, hierárquico e heatmap da via selecionada.
  netVisual_aggregate(cellchat, signaling = pathway.use, layout = "circle")

  # Selecione receptores que existam nos nomes reais do seu objeto.
  vertex.receiver <- which(levels(cellchat@idents) %in%
                             c("CD4 T", "CD8 T", "B cells"))
  if (length(vertex.receiver) > 0 &&
      length(vertex.receiver) < length(levels(cellchat@idents))) {
    netVisual_aggregate(
      cellchat, signaling = pathway.use,
      layout = "hierarchy", vertex.receiver = vertex.receiver
    )
  } else {
    message("Gráfico hierárquico ignorado: ajuste vertex.receiver.")
  }

  netVisual_heatmap(cellchat, signaling = pathway.use, color.heatmap = "Reds")
  plotGeneExpression(cellchat, signaling = pathway.use)

  # Visualizar UM par ligante-receptor (primeira linha, se houver).
  # Evita índices fixos como [11, ], que podem não existir.
  if (nrow(pares_LR) > 0) {
    LR.use <- pares_LR[1, , drop = FALSE]
    netVisual_individual(
      cellchat, signaling = pathway.use,
      pairLR.use = LR.use, layout = "circle"
    )
  }

  # Papel dos grupos celulares na via escolhida.
  netAnalysis_signalingRole_network(
    cellchat, signaling = pathway.use,
    width = 8, height = 2.5, font.size = 10
  )

  write.csv(interacoes_via, file.path(outdir, "interacoes_via_exemplo.csv"), row.names = FALSE)
} else {
  message("Nenhuma via disponível após os filtros; revise grupos e parâmetros.")
}

# 8) CONCLUSÃO -----------------------------------------------------
# CellChat infere interações a partir da expressão de ligantes/receptores
# e do banco de conhecimento; NÃO prova comunicação experimental direta.
# Analise o desenho experimental, o número de amostras independentes,
# os tamanhos dos grupos e o contexto biológico antes de interpretar.
message("Análise concluída. Resultados em: ", outdir)
