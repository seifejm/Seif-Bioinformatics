# Aula 5 — Redes de coexpressão gênica com hdWGCNA no R

Material de apoio à série **Transcriptoma de Célula Única**, do canal **Seif Bioinformatics — Da Biologia aos Dados**.

**Vídeo da aula:** [INSERIR LINK DO YOUTUBE]  
**Playlist:** [INSERIR LINK DA PLAYLIST]

## Objetivo

Construir uma rede de coexpressão gênica a partir de dados scRNA-seq de PBMC, identificar módulos de genes correlacionados, calcular *module eigengenes* (ME), conectividade kME e explorar *hub genes* usando **hdWGCNA**.

## Conteúdo

1. Preparação do objeto Seurat e seleção de genes.
2. Criação e normalização de metacélulas.
3. Seleção do *soft-thresholding power*.
4. Construção da rede, TOM e detecção de módulos.
5. Cálculo de eigengenes e conectividade kME.
6. Identificação de *hub genes* e visualização dos módulos.

## Arquivos

- `aula_05_hdWGCNA.R`: script didático com comentários em português.
- `data/pbmc_seurat_final.rds`: entrada necessária **não incluída no repositório**; produzida na Aula 2.
- `results/wgcna/`: pasta criada automaticamente com figuras, tabelas e o objeto processado.

## Como executar

Abra um projeto RStudio na pasta desta aula (ou adapte os caminhos ao projeto principal) e utilize a estrutura:

```text
aula-05-hdWGCNA/
├── README.md
├── aula_05_hdWGCNA.R
├── data/
│   └── pbmc_seurat_final.rds     # copiar da Aula 2
└── results/
    └── wgcna/                 # gerado pelo script
```

O objeto de entrada precisa conter `pbmc$celltype` (anotação celular), `pbmc$orig.ident` (amostra), o ensaio `RNA` normalizado e a redução `pca`; para alguns gráficos, também é necessária a redução `umap`.

Instale os pacotes indicados no início do script, ajuste os parâmetros e execute as seções em ordem. Para bases maiores, pode ser necessário mais memória e processamento.

Pacotes principais: `Seurat`, `hdWGCNA`, `WGCNA`, `dplyr`, `tidyr`, `ggplot2`, `patchwork` e `cowplot`.

Documentação: [hdWGCNA](https://smorabit.github.io/hdWGCNA/) e [repositório do projeto](https://github.com/smorabit/hdWGCNA).

## Parâmetros didáticos

| Parâmetro | Valor usado | Interpretação |
|---|---:|---|
| `fraction` | 0.05 | Genes detectados em pelo menos 5% das células |
| `k` | 20 | Número de vizinhos para formar metacélulas |
| `min_cells` | 25 | Mínimo de células por grupo para formar metacélulas |
| `SFT.R.sq` | 0.80 | Critério inicial de seleção de potência (revisar diagnóstico) |
| `minModuleSize` | 30 | Número mínimo de genes por módulo |
| `mergeCutHeight` | 0.50 | Limiar utilizado para unir módulos próximos |

Esses valores **não são universais**; revise-os de acordo com o desenho experimental, a estrutura amostral e o tamanho do conjunto de dados. Uma rede geral de tipos celulares pode capturar diferenças entre populações, além de programas de coexpressão dentro delas.

## Interpretação e cuidados

- Coexpressão **não demonstra regulação direta ou causalidade**.
- O módulo `grey` representa genes não atribuídos a módulos bem definidos e é excluído de algumas tabelas.
- Valores altos de kME identificam candidatos a genes centrais no módulo, mas não comprovam função reguladora.
- A escolha do *soft power* depende da inspeção do ajuste scale-free e da conectividade; o critério automatizado serve apenas como ponto de partida.
- A construção de metacélulas deve respeitar a estrutura das amostras e não substitui replicação biológica.
- Os resultados e a compatibilidade de funções podem variar conforme a versão do `hdWGCNA` e do `Seurat`.

## Resultados gerados

Os arquivos serão armazenados em `results/wgcna/`, incluindo tabela de *soft powers*, módulos, *hub genes*, figuras e `pbmc_hdWGCNA_final.rds`. Os arquivos de dados grandes e a TOM não devem ser versionados no GitHub.

> **Aviso:** código reorganizado para fins educacionais a partir do material da aula. É necessário executar e validar com os dados reais antes de considerar os resultados reproduzidos.
