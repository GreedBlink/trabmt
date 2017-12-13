# Análise fatorial

require(psych)
require(corrplot)

dados <- read.csv2("C:\\Users\\leo-f\\Desktop\\Leonardo\\UFF\\6º Período\\Análise multivariada\\Trabalho\\Dados.csv",
                   encoding = "UTF-8")

corrplot(cor(dados))

# Teste de Bartlett

#H0: A matriz de correlações é uma matriz identidade
#H1: A matriz de correlações não é uma matriz identidade

# Assumindo alfa = 5%

cortest.bartlett(cor(dados), n = nrow(dados))
# p-valor = 0. Rejeitamos H0.
# Podemos então prosseguir com a análise fatorial.

# KMO

cor_parcial <- function (X, ...){
  R <- cor(X, ...)
  RI <- solve(R)
  D <- 1/sqrt(diag(RI))
  Rp <- -RI * (D %o% D)
  diag(Rp) <- 0
  rownames(Rp) <- colnames(Rp) <- colnames(X)
  Rp
}


matriz_cor_parcial <- cor_parcial(dados)

diagonal <- seq(1, by = ncol(dados) + 1, length.out = ncol(dados))
soma <- sum((cor(dados)[-diagonal])^2 )

(kmo <- soma / (soma + sum((matriz_cor_parcial[-diagonal])^2)))

# Valores aceitáveis: kmo \in (0,5; 1)
# Então podemos aplicar análise fatorial


#### Fatores via componentes principais ####

KMO(dados) # MAA = MSA

valores_vetores <- eigen(cor(dados))

componentes <- prcomp(dados, scale. = TRUE)
summary(componentes) # Mais de 90% da proporção da variância acumulada até PC6

plot(1:ncol(dados), (componentes$sdev)^2, type = "b", xlab = "Componentes", 
     ylab = "Variância", pch = 20, cex.axis = 1.3, cex.lab = 1.3)

plot(valores_vetores$values, type = "b", ylab = "Autovalores", main = "Screeplot",
     xlab = "Número de autovalores", pch = 20, axes = F) 
axis(2)
axis(1,at = 1:length(valores_vetores$values))

abline(h=1, col="royalblue", lwd=2)

# Por aqui parece ser melhor escolher até PC5

# Escolhendo o critério da proporção da variância acumulada,
# vamos definir 6 componentes principais.

n_comp_principais <- 6

cargas <- componentes$rotation[, 1:n_comp_principais] %*% 
  diag(componentes$sdev[1:n_comp_principais])
colnames(cargas) <- paste("Fator", 1:n_comp_principais)
cargas


comum <- rowSums(cargas^2)
vespec <- diag(matriz_cor_parcial) - comum
estimat <- cbind(comum, vespec, diag(matriz_cor_parcial))
rownames(estimat) <- colnames(dados)
colnames(estimat) <- c("Comunalidade", "Variância unica", "Variância")
estimat # estimacao das comunalidades e das variancias especificas

# Rotação
varimax(cargas)
