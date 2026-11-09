# Aula 3 — Enriquecimento Funcional com Gene Ontology (GO) no R

Material didático da série **Transcriptoma de Célula Única**, do canal **Seif Bioinformatics**.

Nesta aula, utilizamos genes marcadores de tipos celulares identificados com **Seurat** para investigar termos de **Gene Ontology** super-representados com o pacote **clusterProfiler**.

🎥 **Vídeo da Aula 3 no YouTube:** [INSERIR LINK]

## Conteúdo

- Seleção de genes marcadores a partir de `FindAllMarkers()`.
- Definição do universo (*background*) de genes.
- Enriquecimento GO para **BP** (processo biológico), **MF** (função molecular) e **CC** (componente celular).
- Interpretação de `GeneRatio`, `BgRatio`, `Count` e `p.adjust`.
- Visualização por dotplot, barplot e redes exploratórias.
- Consulta opcional às anotações de `GO.db` e `org.Hs.eg.db`.

## Arquivos

```text
aula-03-gene-ontology/
├── README.md
├── aula_03_gene_ontology.R
├── data/                         # criar localmente
│   ├── pbmc_seurat_final.rds
│   └── markers_celltypes.csv
└── results/                      # gerada automaticamente
```

Os dois arquivos da pasta `data/` são produzidos na **Aula 2**. Não é necessário enviar dados grandes ao GitHub.

## Como executar

1. Abra a pasta da aula como um projeto do RStudio (ou ajuste o diretório de trabalho para ela).
2. Crie a pasta `data/` e coloque nela os arquivos `pbmc_seurat_final.rds` e `markers_celltypes.csv` gerados na Aula 2.
3. Instale os pacotes indicados no início do script, se necessário.
4. Abra `aula_03_gene_ontology.R` e execute os blocos em sequência.
5. Altere `celltype_alvo <- "B cells"` para outro tipo celular presente na sua tabela, caso deseje.
6. Consulte tabelas e figuras salvas automaticamente em `results/`.

**Pacotes principais:** `Seurat`, `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `AnnotationDbi`, `GO.db`, `dplyr` e `ggplot2`.

## Parâmetros do exemplo

| Parâmetro | Valor | Significado |
|---|---|---|
| Tipo celular | `B cells` | Grupo usado na seleção dos marcadores |
| `p_val_adj` | `< 0.05` | Filtro de significância para marcadores |
| `avg_log2FC` | `> 0.25` | Genes com maior expressão no grupo escolhido |
| `pAdjustMethod` | `BH` | Ajuste de múltiplos testes |
| `pvalueCutoff` | `0.05` | Filtro de enriquecimento |
| `qvalueCutoff` | `0.20` | Filtro adicional de q-value |
| Universo | genes do assay RNA reconhecidos na base humana | Aproximação didática dos genes disponíveis para seleção |

## Cuidados na interpretação

- **Enriquecimento não é sinônimo de ativação de via.** O teste identifica termos GO mais representados entre os genes selecionados do que o esperado pelo universo considerado.
- O **universo ideal** corresponde aos genes elegíveis ao teste diferencial que originou a lista. Aqui, os genes do assay RNA são usados como aproximação para fins de ensino.
- As anotações são para **Homo sapiens**. Para outras espécies, troque a base `OrgDb` e confira os identificadores.
- Genes e termos podem mudar conforme filtros, versões de pacotes e bases de anotação.
- Este script não faz inferência de diferenças entre condições ou replicatas biológicas; ele interpreta uma lista de genes marcadores de tipos celulares.

---

**Seif Bioinformatics — Da Biologia aos Dados.**
