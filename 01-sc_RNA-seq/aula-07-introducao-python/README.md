# Aula 7 — Introdução ao Python para Bioinformática

**Seif Bioinformatics — Da Biologia aos Dados**  
**Série:** Transcriptoma de Célula Única

Nesta aula, apresentamos os fundamentos de Python para quem está começando a programar, especialmente estudantes e pesquisadores que já conhecem um pouco de **R** e pretendem utilizar **Scanpy** em análises de dados de célula única (*scRNA-seq*).

O material traz exemplos simples e comentados para acompanhar a aula prática.

## Conteúdo da aula

1. Verificação da instalação do Python e criação de ambiente virtual (`venv`).
2. Instalação de bibliotecas com `pip`.
3. Atribuição de objetos e tipos de dados (`int`, `float`, `str`, `bool`).
4. Listas, indexação a partir de zero e dicionários.
5. Importação de pacotes e uso de apelidos (`import numpy as np`).
6. Funções básicas com NumPy.
7. Criação e exploração de tabelas com Pandas.
8. Funções `len()`, `max()`, `min()`, `sum()` e `print()`.
9. Comparações práticas entre R e Python.
10. Introdução à organização das funções do Scanpy.

## Arquivos

```text
 aula-07-introducao-python/
 ├── aula_07_introducao_python.py
 ├── README.md
 └── .gitignore
```

## Requisitos

- Python 3.12 (versão usada como exemplo na aula; outras versões compatíveis com as bibliotecas também podem funcionar).
- Terminal do Windows, VS Code, Spyder ou Jupyter Notebook.
- NumPy e Pandas para executar os exemplos do script.

## Como executar no Windows

No **PowerShell** ou **Prompt de Comando**, na pasta do projeto:

```powershell
python --version
python -m venv scRNA_env
```

**Prompt de Comando (cmd):**

```bat
scRNA_env\Scripts\activate.bat
```

**PowerShell:**

```powershell
.\scRNA_env\Scripts\Activate.ps1
```

Em seguida, instale as dependências:

```powershell
python -m pip install --upgrade pip
python -m pip install numpy pandas matplotlib scipy scikit-learn jupyter
```

Para executar o script completo:

```powershell
python aula_07_introducao_python.py
```

> **Dica:** a aula foi pensada para executar o código **por partes**, acompanhando a saída de cada seção. No arquivo `.py`, os exemplos utilizam `print()` para que os resultados também apareçam ao executar o script completo no terminal.

Para sair do ambiente virtual:

```powershell
deactivate
```

### Pacotes para as próximas aulas

Quando avançarmos para a análise de scRNA-seq com Scanpy:

```powershell
python -m pip install scanpy anndata
```

Esses dois pacotes **não são necessários** para executar os exemplos básicos desta aula.

## Exemplos utilizados

```python
# Criar uma lista de genes
genes = ["CD3D", "MS4A1", "LYZ"]
print(genes[0])

# Média com NumPy
import numpy as np
print(np.mean([1, 2, 3, 4, 5]))

# Tabela com Pandas
import pandas as pd
dados = pd.DataFrame({
    "gene": ["CD3D", "MS4A1", "LYZ"],
    "expression": [5.2, 8.1, 6.4]
})
print(dados.head())
```

## R × Python: comparação rápida

| Operação | R | Python |
|---|---|---|
| Atribuição | `x <- 10` | `x = 10` |
| Vetor/lista simples | `c(1, 2, 3)` | `[1, 2, 3]` |
| Importar pacote | `library(dplyr)` | `import pandas as pd` |
| Acessar uma coluna | `dados$gene` | `dados["gene"]` |
| Booleanos | `TRUE`, `FALSE` | `True`, `False` |

## Próximos passos

Nas próximas aulas, utilizaremos **Scanpy** para trabalhar com dados de RNA de célula única, incluindo controle de qualidade, normalização, PCA, UMAP, clustering e identificação de genes marcadores.

## Videoaula

- **YouTube:** [Adicionar link da Aula 7]
- **Playlist — Transcriptoma de Célula Única:** [Adicionar link da playlist]

---

**Seif Bioinformatics — Da Biologia aos Dados**  
Conteúdo educacional gratuito sobre bioinformática, programação e análise de dados biológicos.
