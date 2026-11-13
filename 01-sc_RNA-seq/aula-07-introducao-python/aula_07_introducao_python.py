# ============================================================
# AULA 7 - INTRODUÇÃO AO PYTHON
# SEIF BIOINFORMATICS | DA BIOLOGIA AOS DADOS
# ============================================================
#
# Execute este arquivo por partes no VS Code, Spyder ou em um
# notebook Jupyter para acompanhar os resultados de cada etapa.
# Se executar como script completo, utilize print() para mostrar
# os resultados no terminal.
#
# ============================================================
# 0. INSTALAÇÃO E AMBIENTE VIRTUAL
# ============================================================
#
# Os comandos abaixo devem ser executados no Prompt de Comando
# ou PowerShell do Windows, e NÃO no interpretador Python.
#
# Verificar a versão instalada:
# python --version
#
# Instalar o Python pelo winget (Windows, se disponível):
# winget install Python.Python.3.12
#
# Criar um ambiente virtual na pasta do projeto:
# python -m venv scRNA_env
#
# Ativar no Prompt de Comando (cmd):
# scRNA_env\Scripts\activate.bat
#
# Ativar no PowerShell:
# .\scRNA_env\Scripts\Activate.ps1
#
# Quando o ambiente estiver ativo, aparecerá (scRNA_env)
# no início da linha do terminal.
#
# Atualizar o pip:
# python -m pip install --upgrade pip
#
# Instalar os pacotes básicos:
# python -m pip install numpy pandas matplotlib scipy scikit-learn jupyter
#
# Pacotes que utilizaremos nas próximas aulas de single-cell:
# python -m pip install scanpy anndata
#
# Desativar o ambiente:
# deactivate
#
# Um ambiente virtual reúne seus próprios pacotes e dependências,
# ajudando a evitar conflitos entre projetos.


# ============================================================
# 1. ATRIBUIÇÃO DE OBJETOS
# ============================================================
#
# Em R:      x <- 10
# Em Python: x = 10

x = 10
print(x)

nome = "Elias"
print(nome)


# ============================================================
# 2. TIPOS DE OBJETOS
# ============================================================

idade = 30            # int: número inteiro
altura = 1.75         # float: número decimal
cidade = "São Paulo"  # str: texto
ativo = True          # bool: verdadeiro/falso

print(type(idade))
print(type(altura))
print(type(cidade))
print(type(ativo))


# ============================================================
# 3. LISTAS
# ============================================================
#
# Em R:      genes <- c("CD3D", "MS4A1", "LYZ")
# Em Python: usamos colchetes.

genes = ["CD3D", "MS4A1", "LYZ"]

# ATENÇÃO: em Python a indexação começa em zero.
# Índices:      0       1       2
# Valores:   CD3D    MS4A1     LYZ

print(genes)
print(genes[0])  # CD3D
print(genes[1])  # MS4A1
print(genes[2])  # LYZ


# ============================================================
# 4. DICIONÁRIOS
# ============================================================
#
# Dicionários armazenam pares chave -> valor.

celltypes = {
    "cluster_0": "CD4 T",
    "cluster_1": "B cells",
    "cluster_2": "Monocytes"
}

print(celltypes)
print(celltypes["cluster_1"])


# ============================================================
# 5. CARREGAR PACOTES
# ============================================================
#
# Em R:      library(dplyr)
# Em Python: import pacote

import numpy
import pandas

print(numpy.mean([1, 2, 3, 4, 5]))


# ============================================================
# 6. UTILIZAR FUNÇÕES DOS PACOTES
# ============================================================
#
# É comum escrever pacote.funcao().

print(numpy.mean([1, 2, 3, 4, 5]))


# ============================================================
# 7. APELIDOS PARA PACOTES
# ============================================================
#
# import pacote as apelido

import numpy as np

print(np.mean([1, 2, 3, 4, 5]))
print(np.median([1, 2, 3, 4, 5]))
print(np.std([1, 2, 3, 4, 5]))


# ============================================================
# 8. PANDAS: CRIAR UMA TABELA
# ============================================================
#
# pandas é um dos principais pacotes para trabalhar com tabelas.

import pandas as pd

dados = pd.DataFrame({
    "gene": ["CD3D", "MS4A1", "LYZ"],
    "expression": [5.2, 8.1, 6.4]
})

print(dados)


# ============================================================
# 9. EXPLORAR UMA TABELA
# ============================================================

print(dados.head())           # Primeiras linhas
print(dados.columns)          # Nomes das colunas
print(dados.shape)            # (número de linhas, número de colunas)
print(dados["gene"])          # Selecionar coluna de genes
print(dados["expression"])    # Selecionar coluna de expressão


# ============================================================
# 10. FUNÇÕES BÁSICAS DO PYTHON
# ============================================================

numeros = [1, 2, 3, 4, 5]

print(len(numeros))  # Número de elementos
print(max(numeros))  # Maior valor
print(min(numeros))  # Menor valor
print(sum(numeros))  # Soma


# ============================================================
# 11. PRINT
# ============================================================
#
# print() mostra informações no terminal.

print("Hello single-cell!")
print("Número de genes:", len(genes))


# ============================================================
# 12. COMPARAÇÃO R x PYTHON
# ============================================================
#
# ATRIBUIÇÃO
# R:      x <- 10
# Python: x = 10
#
# CARREGAR PACOTE
# R:      library(dplyr)
# Python: import pandas as pd
#
# USAR FUNÇÕES DE PACOTES
# R:      dplyr::filter()
# Python: pd.DataFrame()   (exemplo de função do pandas)
#
# ACESSAR UMA COLUNA
# R:      dados$gene
# Python: dados["gene"]
#
# CRIAR VETOR / LISTA
# R:      c(1, 2, 3)
# Python: [1, 2, 3]
#
# BOOLEANOS
# R:      TRUE / FALSE
# Python: True / False


# ============================================================
# 13. PREPARAÇÃO PARA SINGLE-CELL
# ============================================================
#
# Nas próximas aulas usaremos o Scanpy:
#
# import scanpy as sc
#
# Principais grupos de funções:
# sc.pp -> pré-processamento
# sc.tl -> ferramentas de análise
# sc.pl -> visualizações
#
# Exemplos:
# sc.pp.filter_cells()
# sc.pp.normalize_total()
# sc.pp.log1p()
# sc.pp.highly_variable_genes()
# sc.tl.pca()
# sc.pp.neighbors()
# sc.tl.umap()
# sc.tl.leiden()
# sc.pl.umap()


# ============================================================
# RESUMO
# ============================================================
#
# objeto = valor
# import pacote as apelido
# pacote.funcao()
# lista = [1, 2, 3]
# tabela = pd.DataFrame(...)
# tabela["coluna"]
#
# ============================================================
# FIM
# ============================================================
