# Aula 8 — Análise de scRNA-seq com Scanpy no Python

Material da série **Seif Bioinformatics — Da Biologia aos Dados**. Este tutorial apresenta um pipeline introdutório de RNA-seq de célula única (scRNA-seq) usando **Scanpy**, do carregamento de dados 10X à identificação de marcadores e à anotação celular.

> **Atenção:** os dados 10X não estão incluídos neste repositório. A atribuição de tipos celulares no script é **um exemplo didático** e precisa ser revisada conforme os clusters e marcadores obtidos no seu próprio conjunto de dados.

## Conteúdo

1. Leitura de matriz 10X com `sc.read_10x_mtx()`.
2. Exploração do objeto `AnnData` (`.obs`, `.var`, `.obsm`, `.layers`, etc.).
3. Filtragem inicial e métricas de controle de qualidade.
4. Normalização, transformação logarítmica e seleção de genes altamente variáveis.
5. Escalonamento, PCA e avaliação da variância explicada.
6. Construção do grafo de vizinhos, Leiden e UMAP.
7. Identificação de marcadores por cluster com teste de Wilcoxon.
8. Exploração de marcadores conhecidos e anotação de tipos celulares.
9. Marcadores por tipo celular e exportação de figuras, tabelas e `AnnData`.

## Estrutura sugerida

```text
05-r-python/
└── aula-08-scanpy/
    ├── README.md
    ├── aula_08_scanpy.py
    ├── requirements.txt
    ├── .gitignore
    ├── data/
    │   └── filtered_feature_bc_matrix/
    │       ├── matrix.mtx.gz
    │       ├── features.tsv.gz
    │       └── barcodes.tsv.gz
    └── results/
        └── python/
            └── open/
```

O conteúdo de `data/` e os resultados não precisam ser publicados no GitHub. O script cria automaticamente a pasta de resultados.

## Instalação no Windows

Abra o terminal na pasta da aula e crie um ambiente virtual:

```powershell
py -3.12 -m venv scRNA_env
```

Ative no **PowerShell**:

```powershell
.\scRNA_env\Scripts\Activate.ps1
```

Ou no **Prompt de Comando (cmd)**:

```bat
scRNA_env\Scripts\activate.bat
```

Instale as dependências:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

> O ambiente Python precisa ter os pacotes `python-igraph` e `leidenalg` para executar o clustering Leiden configurado com `flavor="igraph"`.

## Dados de entrada

Use uma saída **10X Genomics / Cell Ranger** com a matriz filtrada, contendo os três arquivos na pasta `data/filtered_feature_bc_matrix/`. O exemplo foi escrito para **genes humanos**, pois identifica genes mitocondriais pelo prefixo `MT-`. Para camundongos ou outras espécies, adapte esse padrão conforme a anotação gênica.

## Execução

```bash
python aula_08_scanpy.py
```

O script também pode ser executado por seções em um editor Python apropriado.

## Resultados

Arquivos criados em `results/python/open/`:

- `QC_violin.png` — distribuição das métricas de QC.
- `highly_variable_genes.png` — genes altamente variáveis.
- `pca_variance_ratio.png` e `pca_plot.png` — PCA.
- `UMAP_clusters_leiden.png` — clusters.
- `UMAP_marker_genes.png` e `DotPlot_marker_genes.png` — genes marcadores conhecidos.
- `DotPlot_top5_markers_clusters.png` — genes marcadores por cluster.
- `UMAP_celltypes.png`, `DotPlot_top10_celltypes.png` e `DotPlot_top5_markers_celltypes.png` — anotação e marcadores.
- `markers_clusters.csv`, `markers_celltypes.csv`, `top10_markers_celltypes.csv` — tabelas.
- `pbmc_scanpy_final.h5ad` — objeto final.

## Observações metodológicas

- Os limiares de QC mostrados são didáticos: ajuste-os com base na qualidade e na origem dos seus dados.
- A seleção de genes contendo `CD` serve apenas como exemplo de categoria gênica; **não é uma métrica de QC convencional**.
- O script mantém uma cópia da expressão normalizada em `adata.raw` antes de restringir a matriz aos 2.000 genes variáveis, permitindo buscar marcadores em todos os genes retidos após QC.
- O `logfoldchanges` calculado pelo Scanpy em `rank_genes_groups` é uma estimativa de fold-change baseada nas médias de expressão log-transformada e não deve ser confundido com o resultado de modelos específicos de contagens.
- Os grupos numéricos produzidos por Leiden dependem dos dados e dos parâmetros. **Nunca reutilize a tabela `celltype_map` sem conferir os marcadores e os agrupamentos resultantes.**
- O código não executa correção de *batch*: essa etapa é mencionada como possibilidade para análises com múltiplas amostras.

## Materiais

- 🎥 **Videoaula:** [adicionar link]
- 📚 **Playlist:** [adicionar link]
- 📖 Documentação oficial: <https://scanpy.readthedocs.io/>
- 📖 AnnData: <https://anndata.readthedocs.io/>

**Seif Bioinformatics — Da Biologia aos Dados**
