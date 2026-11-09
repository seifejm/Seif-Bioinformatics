# Aula 01 — Exemplo didático: organização de um projeto em RStudio
# Este arquivo é um modelo, não o script original da aula.

# 1. Criar pastas a partir da raiz do projeto .Rproj
pastas <- c("data", "results", "results/figures", "functions", "scripts")
invisible(lapply(pastas, dir.create, recursive = TRUE, showWarnings = FALSE))

# 2. Usar caminhos relativos (evite caminhos absolutos do computador)
caminho_dados <- file.path("data", "meus_dados.csv")

# 3. Carregar funções reutilizáveis quando o arquivo existir
arquivo_funcoes <- file.path("functions", "funcoes_auxiliares.R")
if (file.exists(arquivo_funcoes)) source(arquivo_funcoes)

# 4. Salvar resultados para evitar reprocessamento
resumo_demo <- data.frame(etapa = "organização", concluida = TRUE)
saveRDS(resumo_demo, file.path("results", "resumo_demo.rds"))

# 5. Importante: arquivos grandes/dados sensíveis não devem ser enviados ao GitHub
