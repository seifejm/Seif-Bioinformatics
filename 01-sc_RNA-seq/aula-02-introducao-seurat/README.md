# Aula 02 — Análise de scRNA-seq com Seurat

Material complementar da série **Transcriptômica de Célula Única**, do canal **Seif Bioinformatics — Da Biologia aos Dados**.

Nesta aula, realizamos um fluxo introdutório de análise de dados de RNA de célula única (scRNA-seq) no **R**, usando o pacote **Seurat**: leitura de dados 10x Genomics, criação do objeto, controle de qualidade, normalização, genes variáveis, PCA, clusterização, UMAP, marcadores e anotação de tipos celulares.

📺 **Vídeo no YouTube:** [INSERIR LINK DA AULA]

## Arquivo da aula

- [`aula_02_seurat_pipeline.R`](aula_02_seurat_pipeline.R): script completo, comentado em português e organizado por etapas.

## Requisitos

- R e RStudio instalados.
- Pacotes R: `Seurat`, `dplyr` e `ggplot2`.
- Dados de expressão gênica no formato **10x Genomics**, preferencialmente PBMC 1k (genes humanos), conforme a demonstração da aula.

Para instalar os pacotes no R:

```r
install.packages(c("Seurat", "dplyr", "ggplot2"))
```

## Como preparar os dados

Crie um projeto no RStudio com esta estrutura:

```text
meu-projeto/
├── aula_02_seurat_pipeline.R
├── data/
│   └── filtered_feature_bc_matrix/
│       ├── matrix.mtx.gz
│       ├── features.tsv.gz
│       └── barcodes.tsv.gz
└── results/                     # criado automaticamente pelo script
```

Os arquivos podem estar descompactados (`.mtx`/`.tsv`). Caso você use outro conjunto de dados, altere `diretorio_10x` no início do script. Os dados não acompanham este repositório; baixe-os da fonte original (por exemplo, [10x Genomics Datasets](https://www.10xgenomics.com/datasets)).

## Etapas do pipeline

1. **`Read10X()`**: importa a matriz de contagens (genes × células).
2. **`CreateSeuratObject()`**: organiza os dados no Seurat e aplica filtros *mínimos* iniciais (`min.cells` e `min.features`).
3. **`PercentageFeatureSet()` / `VlnPlot()`**: avaliam contagens, genes detectados e fração mitocondrial.
4. **`NormalizeData()`**: normaliza as contagens, reduzindo diferenças de profundidade de sequenciamento entre células.
5. **`FindVariableFeatures()` / `ScaleData()`**: selecionam genes informativos e padronizam os valores usados no PCA.
6. **`RunPCA()` / `ElbowPlot()`**: resumem a variação dos dados e ajudam a escolher os componentes principais.
7. **`FindNeighbors()` / `FindClusters()`**: constroem o grafo de vizinhança e identificam grupos de células.
8. **`RunUMAP()`**: visualiza as células em duas dimensões; **não atribui os clusters**.
9. **`FindAllMarkers()`**: busca genes com maior expressão em cada cluster em comparação aos demais.
10. **`FeaturePlot()` / `DotPlot()`**: ajudam a avaliar marcadores canônicos para a anotação biológica.
11. **Exportação**: salva tabelas `.csv` e o objeto final `.rds` na pasta `results/`.

## Atenção antes de executar

- **Filtragem de QC:** no script, `aplicar_filtro_qc <- FALSE`. Os valores de 200, 2500 e 5% são **ilustrativos**, não recomendações universais. Avalie o seu próprio conjunto de dados antes de alterar para `TRUE`.
- **Marcadores mitocondriais:** `^MT-` é apropriado para a nomenclatura usual de genes humanos; para outras espécies, confira o prefixo.
- **PCA:** o exemplo calcula 50 PCs, mas usa os PCs 1 a 10 para clustering e UMAP. Essa escolha deve ser revista conforme o `ElbowPlot` e a biologia das amostras.
- **Anotação celular:** os nomes dos clusters no exemplo original foram preservados somente como **referência didática** em `mapa_exemplo`. A opção `anotacao_validada <- FALSE` impede que rótulos possivelmente incorretos sejam atribuídos automaticamente. Valide os marcadores e atualize os IDs antes de mudar para `TRUE`.
- **Expressão diferencial:** os testes célula a célula de `FindAllMarkers()` servem para identificar marcadores exploratórios. Comparações entre condições biológicas com réplicas exigem atenção ao desenho experimental, por exemplo com análises pseudobulk.
- **Integração de amostras:** não foi incluída nesta aula introdutória de uma amostra.

## Saídas

A pasta `results/` conterá:

- `markers_clusters.csv` — marcadores por cluster;
- `top10_markers_clusters.csv` — até 10 marcadores significativos por cluster, quando disponíveis;
- `pbmc_seurat_final.rds` — objeto Seurat para continuar as análises;
- `markers_celltypes.csv` — gerado **apenas após validação da anotação**.

## Sobre o material

O script foi simplificado a partir do código utilizado em aula: removemos a segunda execução redundante de normalização, PCA, clusterização, UMAP e a repetição de `FindAllMarkers()` para os mesmos grupos. Comentários explicativos foram mantidos ou ampliados para facilitar o acompanhamento dos vídeos.

**Seif Bioinformatics — Da Biologia aos Dados.**
