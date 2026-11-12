# Aula 6 — Análise de Pseudotempo com Monocle 3 no R

**Seif Bioinformatics — Da Biologia aos Dados**  
**Série:** Transcriptoma de Célula Única (scRNA-seq)

Nesta aula, usamos **Monocle 3** para estudar a heterogeneidade transcriptômica dentro de grupos de células PBMC, inferir trajetórias, definir um ponto inicial (*root node*), calcular **pseudotempo**, identificar genes associados à estrutura do grafo usando **Moran's I** e explorar **módulos de genes** e suas dinâmicas de expressão.

> **Pseudotempo não é tempo cronológico.** É uma ordenação relativa das células ao longo de uma trajetória inferida. A escolha da raiz e a interpretação das transições exigem justificativa biológica.

## Materiais

- **Código:** [`aula_06_monocle3.R`](aula_06_monocle3.R)
- **Vídeo da Aula 6:** [adicionar link do YouTube]
- **Playlist — Transcriptoma de Célula Única:** [adicionar link da playlist]

## Estrutura sugerida

```text
01-sc_RNA-seq/
└── aula-06-monocle3/
    ├── README.md
    ├── aula_06_monocle3.R
    ├── .gitignore
    ├── data/
    │   └── pbmc_seurat_final.rds   # fornecer localmente; não incluído
    └── results/
        └── monocle3/             # criado pelo script
```

## Pré-requisitos

- R e RStudio instalados.
- Um objeto **Seurat** com ensaio `RNA`, redução `umap` e anotação de tipos celulares na coluna `celltype`.
- O objeto `pbmc_seurat_final.rds` produzido na **Aula 2** (ou um objeto equivalente adaptando o mapeamento de tipos celulares).
- Pacotes `Seurat`, `SeuratWrappers`, `monocle3`, `SingleCellExperiment`, `dplyr`, `tibble`, `ggplot2`, `pheatmap`, `Matrix` e `igraph` (e dependências dos pacotes).

As instalações podem variar conforme versão do R e sistema operacional. Um exemplo inicial:

```r
install.packages(c("Seurat", "remotes", "dplyr", "tibble", "ggplot2", "pheatmap", "Matrix", "igraph"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
remotes::install_github("cole-trapnell-lab/monocle3")
remotes::install_github("satijalab/seurat-wrappers")
```

Confira as instruções atualizadas nos projetos oficiais: [Monocle 3](https://cole-trapnell-lab.github.io/monocle3/) e [SeuratWrappers](https://github.com/satijalab/seurat-wrappers).

## Como executar

1. Abra o projeto ou defina como diretório de trabalho a pasta `aula-06-monocle3/`.
2. Coloque `pbmc_seurat_final.rds` em `data/`.
3. Abra `aula_06_monocle3.R` no RStudio.
4. Configure as variáveis no início do arquivo, principalmente `GRUPO_ESCOLHIDO`.
5. **Primeira execução:** mantenha `ROOT_NODE <- NULL`. O script constrói o grafo, salva as figuras dos nós e das partições e interrompe a execução intencionalmente.
6. Abra `results/monocle3/02_nos_grafo.png` e `03_particoes.png`. Escolha um nó `Y_...` que represente biologicamente o estado inicial; não escolha a raiz somente por sua posição geométrica.
7. **Segunda execução:** substitua `ROOT_NODE <- NULL` pelo nó escolhido (por exemplo, `ROOT_NODE <- "Y_1"`, somente se esse nó existir e fizer sentido). Execute novamente.
8. Verifique os arquivos gerados em `results/monocle3/`.

### Parâmetros principais

| Parâmetro | Padrão | Função |
|---|---|---|
| `GRUPO_ESCOLHIDO` | `"Monocytes"` | Grupo celular analisado |
| `ROOT_NODE` | `NULL` | Nó inicial a ser definido após inspecionar o grafo |
| `N_BINS` | `30` | Número de grupos ordenados de células ao longo do pseudotempo |
| `QVALUE_CUTOFF` | `0.05` | Limite de significância ajustada no `graph_test()` |
| `MORANS_CUTOFF` | `0.01` | Limite mínimo de Moran's I |
| `MODULE_RESOLUTION` | `1e-3` | Resolução de `find_gene_modules()` |
| `N_TOP_GENES` | `5` | Quantidade de genes selecionados por módulo para o heatmap |

### Grupos celulares utilizados

O código mapeia os nomes da anotação original:

| Grupo de trajetória | Tipos celulares do Seurat |
|---|---|
| `Monocytes` | `Mono (CD14)`, `Mono (FCGR3A)` |
| `B_cells` | `B cells` |
| `CD4_T` | `CD4 T` |
| `CD8_NK` | `CD8 T`, `NK cells` |

Células não mapeadas são excluídas dessa análise. A divisão é **didática**: agrupar `CD8 T` e `NK cells` não implica que essas células formem uma única linhagem de diferenciação. Avalie a coerência biológica de cada trajetória e considere analisar populações separadamente.

## Etapas da análise

1. **Leitura e anotação:** carregar o Seurat e definir grupos celulares.
2. **Conversão:** criar um `cell_data_set` com `as.cell_data_set()`.
3. **Grafo:** executar `cluster_cells()` e `learn_graph()`.
4. **Raiz:** selecionar o nó inicial e ordenar com `order_cells()`.
5. **Pseudotempo:** visualizar a distribuição entre células e tipos celulares.
6. **Genes associados ao grafo:** usar `graph_test(neighbor_graph = "principal_graph")` com Moran's I.
7. **Módulos Monocle 3:** agrupar genes usando `find_gene_modules()`.
8. **Dinâmica de expressão:** agregar módulos e calcular médias de genes em intervalos de pseudotempo.
9. **Visualização:** gerar gráficos de trajetória, violino e heatmaps.

## Resultados

Em `results/monocle3/`, o script pode gerar:

- `01_grupos_celulares.png`, `02_nos_grafo.png`, `03_particoes.png` — exploração inicial.
- `04_pseudotempo_umap.png`, `05_pseudotempo_por_tipo_celular.png` — pseudotempo.
- `graph_test_todos_genes.csv`, `genes_associados_ao_grafo.csv` — estatísticas gênicas.
- `modulos_genicos.csv`, `top_genes_por_modulo.csv` — genes por módulo.
- `matriz_modulos_pseudotempo.csv`, `matriz_top_genes.csv` — dados dos heatmaps.
- `06_heatmap_modulos.png`, `07_heatmap_top_genes.png`, `08_genes_no_pseudotempo.png` — visualizações, se houver genes/módulos suficientes.
- `monocle3_<grupo>.rds` — objeto Monocle 3 processado.

## Cuidados de interpretação

- Os módulos identificados por `find_gene_modules()` **não equivalem automaticamente** aos módulos de coexpressão identificados pelo **hdWGCNA** na Aula 5.
- `graph_test()` avalia associação da expressão com a estrutura do grafo; isso não demonstra causalidade nem garante uma alteração monotônica com o pseudotempo.
- Se houver **partições desconectadas**, células sem uma raiz conectada podem ficar com pseudotempo infinito e são excluídas das agregações.
- Os heatmaps com `scale = "row"` exibem **expressão relativa padronizada por gene ou módulo**, não valores absolutos diretamente comparáveis entre genes.
- Os intervalos `PT_01`, `PT_02`, etc. agrupam células por **quantis de pseudotempo** e não representam unidades cronológicas equivalentes.
- A validade biológica depende de QC, amostragem, escolha de populações, estrutura do grafo e raiz.

## Reprodutibilidade

O script define uma semente aleatória (`set.seed(12345)`). Para um projeto científico, registre também as versões utilizadas (`sessionInfo()`), documente a origem do objeto de entrada e valide biologicamente as trajetórias.

**Nota:** este material foi reorganizado a partir do roteiro da aula. Não foi executado com o arquivo RDS de entrada; teste e revise os resultados no RStudio antes de utilizá-los em análises científicas.

---

**Seif Bioinformatics — Da Biologia aos Dados**
