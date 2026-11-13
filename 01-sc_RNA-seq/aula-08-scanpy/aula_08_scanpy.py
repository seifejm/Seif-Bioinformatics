# ============================================================
# PIPELINE BÁSICO DE SINGLE-CELL RNA-SEQ COM SCANPY
# ============================================================

# ============================================================
# 1. CARREGAR PACOTES
# ============================================================

import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

print("Scanpy:", sc.__version__)
print("NumPy:", np.__version__)

# Criar pasta onde todas as figuras serão salvas
RESULTS_DIR = Path("results/python/open")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# 2. ABRIR DADOS 10X
# ============================================================

# Equivalente ao Read10X() + CreateSeuratObject() do Seurat
# Coloque matrix.mtx.gz, features.tsv.gz e barcodes.tsv.gz em data/filtered_feature_bc_matrix/
adata = sc.read_10x_mtx(
    "data/filtered_feature_bc_matrix/",
    var_names="gene_symbols",
    make_unique=True
)

adata
adata.shape
adata.n_obs
adata.n_vars

# Genes
adata.var_names[:50]

# Barcodes
adata.obs_names[:10]



# ============================================================
# ESTRUTURA DO OBJETO AnnData
# ============================================================
#
# SEURAT                     AnnData
# ------------------------------------------------------------
# objeto Seurat              adata
# pbmc@meta.data             adata.obs
# assay RNA                  adata.X / adata.layers
# VariableFeatures()         adata.var["highly_variable"]
# Embeddings(PCA)            adata.obsm["X_pca"]
# Embeddings(UMAP)           adata.obsm["X_umap"]
# Loadings(PCA)              adata.varm["PCs"]
# Graphs()                   adata.obsp
# misc / resultados          adata.uns
#
# Principais partes:
#
# adata.X       -> matriz células x genes
# adata.obs     -> metadados das células
# adata.var     -> metadados dos genes
# adata.layers  -> matrizes adicionais
# adata.obsm    -> PCA, UMAP etc.
# adata.varm    -> loadings dos genes
# adata.obsp    -> grafos célula-célula
# adata.uns     -> resultados e parâmetros
# adata.raw     -> cópia opcional dos dados

# Explorar objeto
adata.obs.head()
adata.var.head()
adata.obs.columns
adata.var.columns
list(adata.layers.keys())
list(adata.obsm.keys())
list(adata.varm.keys())
list(adata.obsp.keys())
list(adata.uns.keys())

# ============================================================
# 3. FILTRAGEM INICIAL
# ============================================================

# Equivalente a min.cells = 3
sc.pp.filter_genes(
    adata,
    min_cells=3
)

# Equivalente a min.features = 200
sc.pp.filter_cells(
    adata,
    min_genes=200
)

adata.shape

# Guardar contagens brutas
adata.layers["counts"] = adata.X.copy()

# ============================================================
# 4. CONTROLE DE QUALIDADE
# ============================================================

# Genes mitocondriais
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# Genes que contêm "CD" no símbolo (apenas exemplo didático; não é uma métrica de QC padrão)
adata.var["CD"] = adata.var_names.str.contains("CD")

# Calcular métricas
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "CD"],
    percent_top=None,
    inplace=True
)

adata.obs.head()

# ============================================================
# VIOLIN PLOT DO QC
# ============================================================

sc.pl.violin(
    adata,["n_genes_by_counts","total_counts",
           "pct_counts_mt","pct_counts_CD"],
    jitter=0.4,
    multi_panel=True,
    show=False
    )

fig = plt.gcf()
fig.set_size_inches(5, 5)

plt.subplots_adjust(wspace=0.8, hspace=0.4)

plt.savefig("results/python/open/QC_violin.png",
    dpi=300,bbox_inches="tight")

plt.show()
plt.close()

# Filtragem opcional
# adata = adata[
#     (adata.obs["n_genes_by_counts"] > 200) &
#     (adata.obs["n_genes_by_counts"] < 2500) &
#     (adata.obs["pct_counts_mt"] < 5)
# ].copy()

# ============================================================
# 5. NORMALIZAÇÃO
# ============================================================

# Normalizar profundidade por célula
sc.pp.normalize_total(
    adata,
    target_sum=1e4
)

# Aplicar log(1+x)
sc.pp.log1p(
    adata
)

# Guardar dados normalizados de todos os genes
adata.raw = adata.copy()

# ============================================================
# 6. GENES ALTAMENTE VARIÁVEIS
# ============================================================

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat"
)

# Ver genes variáveis
adata.var[[ "highly_variable","means",
        "dispersions","dispersions_norm"]].head()

# Quantidade selecionada
adata.var["highly_variable"].sum()

# ============================================================
# PLOT DOS GENES ALTAMENTE VARIÁVEIS
# ============================================================

sc.pl.highly_variable_genes(
    adata,
    show=False
)

fig = plt.gcf()
fig.set_size_inches(8, 5)

plt.savefig(
    "results/python/open/highly_variable_genes.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
plt.close()

# Manter genes variáveis para redução dimensional
adata = adata[:,adata.var["highly_variable"]].copy()

# ============================================================
# 7. SCALE DATA
# ============================================================

sc.pp.scale(adata,max_value=10)

# ============================================================
# 8. PCA
# ============================================================

N_PCS = min(50, adata.n_obs - 1, adata.n_vars - 1)
if N_PCS < 2:
    raise ValueError("Dados insuficientes para executar PCA. Confira a filtragem.")
sc.tl.pca(adata, n_comps=N_PCS)

# Coordenadas das células
adata.obsm["X_pca"][:10, :5]

# Loadings dos genes
adata.varm["PCs"][:10, :5]

# ============================================================
# ELBOW PLOT / PCA VARIANCE RATIO
# ============================================================

sc.pl.pca_variance_ratio(
    adata,
    n_pcs=N_PCS,
    log=False,
    show=False
)

fig = plt.gcf()
fig.set_size_inches(7, 5)

plt.savefig(
    "results/python/open/pca_variance_ratio.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
plt.close()

# ============================================================
# PCA PLOT
# ============================================================

sc.pl.pca(
    adata,
    show=False
)

fig = plt.gcf()
fig.set_size_inches(7, 5)

plt.savefig(
    "results/python/open/pca_plot.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
plt.close()

# ============================================================
# 9. GRAFO DE VIZINHANÇA
# ============================================================

sc.pp.neighbors(
    adata,
    n_neighbors=15,
    n_pcs=min(10, N_PCS)
)

# Ver grafos armazenados
adata.obsp

# ============================================================
# 10. CLUSTERIZAÇÃO
# ============================================================

# Para usar flavor="igraph", instale python-igraph e leidenalg.
sc.tl.leiden(
    adata,
    resolution=0.5,
    flavor="igraph",
    n_iterations=2,
    key_added="leiden"
)

# Quantidade de células por cluster
adata.obs["leiden"].value_counts().sort_index()

# ============================================================
# 11. UMAP
# ============================================================

sc.tl.umap(
    adata,
    min_dist=0.3,
    spread=1.0,
    random_state=123
)

# ============================================================
# UMAP DOS CLUSTERS
# ============================================================

sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="on data",
    legend_fontweight="bold",
    frameon=False,
    show=False
)

fig = plt.gcf()
fig.set_size_inches(7, 6)

plt.savefig(
    "results/python/open/UMAP_clusters_leiden.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
plt.close()

# ============================================================
# 12. INTEGRAÇÃO / CORREÇÃO DE BATCH
# ============================================================

# Se houver apenas uma amostra esta etapa não é necessária.
#
# Para múltiplas amostras podemos utilizar:
#
# Harmony
# Scanorama
# scVI
# BBKNN

# ============================================================
# 13. MARCADORES DE CADA CLUSTER
# ============================================================

# Equivalente ao FindAllMarkers()
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon",
    use_raw=True,
    key_added="markers_clusters"
)

# Transformar em DataFrame
markers = sc.get.rank_genes_groups_df(
    adata,
    group=None,
    key="markers_clusters"
)

# Filtrar marcadores
markers = markers[
    (markers["logfoldchanges"] > 0.25) &
    (markers["pvals_adj"] < 0.05)
].copy()

markers.head()

# ============================================================
# 14. AVALIAR MARCADORES CONHECIDOS
# ============================================================

marker_genes = ["CD3D", "CD4","CD8A",  "MS4A1",
                "CD14","FCGR3A","GNLY","NKG7"]

marker_genes = [gene for gene in marker_genes if gene in adata.raw.var_names]
if not marker_genes:
    raise ValueError("Nenhum gene marcador encontrado no objeto. Confira gene_symbols.")

# ============================================================
# FEATURE PLOT DOS MARCADORES
# ============================================================

sc.pl.umap(
    adata,
    color=marker_genes,
    use_raw=True,
    ncols=4,
    frameon=False,
    show=False
)

fig = plt.gcf()
fig.set_size_inches(14, 7)

plt.savefig(
    "results/python/open/UMAP_marker_genes.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
plt.close()

# ============================================================
# DOTPLOT DOS MARCADORES
# ============================================================

plot_obj = sc.pl.dotplot(
    adata,
    var_names=marker_genes,
    groupby="leiden",
    use_raw=True,
    standard_scale="var",
    show=False,
    return_fig=True
)

plot_obj.savefig("results/python/open/DotPlot_marker_genes.png")
plt.close("all")

# ============================================================
# TOP MARCADORES POR CLUSTER
# ============================================================

plot_obj = sc.pl.rank_genes_groups_dotplot(
    adata,
    key="markers_clusters",
    n_genes=5,
    groupby="leiden",
    show=False,
    return_fig=True
)

plot_obj.savefig("results/python/open/DotPlot_top5_markers_clusters.png")
plt.close("all")

# ============================================================
# 15. ANOTAÇÃO DOS TIPOS CELULARES
# ============================================================

# Ver clusters antes da anotação
adata.obs["leiden"].value_counts().sort_index()

# IMPORTANTE:
# o mapeamento abaixo é APENAS ilustrativo!
# ajuste os números conforme o resultado do Leiden,
# FeaturePlot, DotPlot e genes diferenciais.

celltype_map = {
    "0": "Mono (CD14)",
    "1": "Unknown",
    "2": "B cells",
    "3": "Unknown",
    "4": "CD8 T",
    "5": "CD4 T",
    "6": "CD4 T",
    "7": "Mono (FCGR3A)",
    "8": "NK cells",
    "9": "Mono (CD14)",
    "10": "Unknown"
}

# Criar coluna celltype
adata.obs["celltype"] = (
    adata.obs["leiden"]
    .astype(str)
    .map(celltype_map)
    .fillna("Unknown")
)

adata.obs[["leiden","celltype"]].head()

adata.obs["celltype"].value_counts()

# ============================================================
# UMAP ANOTADO
# ============================================================

sc.pl.umap(
    adata,
    color="celltype",
    legend_loc="on data",
    legend_fontweight="bold",
    frameon=False,
    show=False
)

fig = plt.gcf()
fig.set_size_inches(8, 6)

plt.savefig(
    "results/python/open/UMAP_celltypes.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
plt.close()

# ============================================================
# 16. MARCADORES POR TIPO CELULAR
# ============================================================

sc.tl.rank_genes_groups(
    adata,
    groupby="celltype",
    method="wilcoxon",
    use_raw=True,
    key_added="markers_celltype"
)

markers_celltype = sc.get.rank_genes_groups_df(
    adata,
    group=None,
    key="markers_celltype"
)

# Filtrar
markers_celltype = markers_celltype[
    (markers_celltype["logfoldchanges"] > 0.25) &
    (markers_celltype["pvals_adj"] < 0.05)
].copy()

markers_celltype.head()

# ============================================================
# 17. TOP 10 GENES POR TIPO CELULAR
# ============================================================

top10 = (markers_celltype
    .sort_values(
        ["group", "logfoldchanges"],
        ascending=[True, False]
    )
    .groupby(
        "group",
        observed=True
    )
    .head(10)
)

top10

# Extrair genes
top_genes = (
    top10["names"]
    .drop_duplicates()
    .tolist()
)

# ============================================================
# DOTPLOT DOS TOP 10 GENES
# ============================================================

plot_obj = sc.pl.dotplot(
    adata,
    var_names=top_genes,
    groupby="celltype",
    use_raw=True,
    standard_scale="var",
    show=False,
    return_fig=True
)

plot_obj.savefig("results/python/open/DotPlot_top10_celltypes.png")
plt.close("all")

# ============================================================
# TOP MARCADORES POR TIPO CELULAR
# ============================================================

plot_obj = sc.pl.rank_genes_groups_dotplot(
    adata,
    key="markers_celltype",
    n_genes=5,
    groupby="celltype",
    show=False,
    return_fig=True
)

plot_obj.savefig("results/python/open/DotPlot_top5_markers_celltypes.png")
plt.close("all")

# ============================================================
# 18. SALVAR RESULTADOS
# ============================================================

markers.to_csv(
    "results/python/open/markers_clusters.csv",
    index=False
)

markers_celltype.to_csv(
    "results/python/open/markers_celltypes.csv",
    index=False
)

top10.to_csv(
    "results/python/open/top10_markers_celltypes.csv",
    index=False
)

# Salvar objeto AnnData
adata.write_h5ad(
    "results/python/open/pbmc_scanpy_final.h5ad"
)

# Para abrir novamente:
#
# adata = sc.read_h5ad(
#     "results/python/open/pbmc_scanpy_final.h5ad"
# )

# ============================================================
# FIM DO PIPELINE
# ============================================================
