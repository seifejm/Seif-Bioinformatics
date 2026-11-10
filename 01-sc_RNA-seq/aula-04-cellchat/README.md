# Aula 4 — Comunicação celular com CellChat no R

Material da série **Transcriptômica de Célula Única** do canal **Seif Bioinformatics**.

Nesta aula, utilizamos um objeto Seurat com dados humanos de PBMC para explorar a **inferência de comunicação célula–célula** por meio de interações ligante–receptor e vias de sinalização com o pacote **CellChat**.

> **Escopo:** este repositório contém somente o conteúdo de **CellChat**. A parte de LIANA foi retirada para corresponder ao conteúdo publicado nesta aula.

## Conteúdo

1. Preparação do objeto Seurat e anotação `celltype`.
2. Criação do objeto CellChat e exploração do banco `CellChatDB.human`.
3. Categorias: *Secreted Signaling*, *Cell-Cell Contact* e *ECM-Receptor*.
4. Identificação de genes e interações, inferência de comunicação e agregação por vias.
5. Exportação de tabelas de interações e de vias.
6. Visualizações globais (número e força das interações).
7. Exploração de uma via, gráficos circular/hierárquico, heatmap e pares ligante–receptor.

## Arquivos

```text
aula-04-cellchat/
├── README.md
└── aula_04_cellchat.R
```

O script cria, na execução, a pasta `results/cellchat/` com o objeto processado (`pbmc_cellchat.rds`) e tabelas CSV.

## Pré-requisitos

- R e RStudio.
- Pacotes `Seurat`, `CellChat`, `dplyr` e `ggplot2`.
- Um objeto Seurat humano salvo em `data/pbmc_seurat_final.rds`, com assay `RNA`, dados normalizados e coluna `celltype` nos metadados.
- Para usar o mesmo exemplo, conclua antes a **Aula 2 — Análise de scRNA-seq com Seurat**.

Instalação de pacotes, se necessária:

```r
install.packages(c("Seurat", "dplyr", "ggplot2", "remotes"))
remotes::install_github("jinworks/CellChat")
```

Confira a documentação da versão instalada do CellChat. Dependências e requisitos podem variar.

## Como executar

1. Abra um projeto no RStudio com a pasta `data/` contendo o arquivo `pbmc_seurat_final.rds`.
2. Ajuste o caminho `arquivo_pbmc` no início do script, se necessário.
3. Confira os nomes e a quantidade de células por tipo em `pbmc$celltype`.
4. Execute `aula_04_cellchat.R` **na ordem**. A inferência pode exigir tempo de processamento.
5. Consulte as tabelas geradas em `results/cellchat/` e explore as vias presentes em `cellchat@netP$pathways`.

**Nota:** os arquivos `.rds` e as tabelas de resultados não estão incluídos aqui, porque dependem da execução com os dados utilizados pelo aluno. Se a sua versão do Seurat utiliza múltiplas camadas (*layers*) de RNA, revise essa configuração antes de construir o CellChat.

## Cuidados na interpretação

O CellChat calcula **interações inferidas** a partir de perfis de expressão e de um banco de ligantes e receptores. Esses resultados **não comprovam, isoladamente, comunicação celular direta**. O contexto biológico, a composição celular e a existência de amostras independentes devem orientar a interpretação.

## Videoaula

🎬 **Aula 4 — Comunicação Celular com CellChat no R**  
🔗 YouTube: **[INSERIR LINK DO VÍDEO]**

📚 Playlist — Transcriptômica de Célula Única: **[INSERIR LINK]**

---

**Seif Bioinformatics — Da Biologia aos Dados.**
