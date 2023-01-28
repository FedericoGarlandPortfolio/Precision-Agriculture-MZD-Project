# ----------------- VARIABILIDAD ESPACIAL DE SUELO Y DELINEACIÓN DE ZONAS DE MANEJO ----------------------
# Cargando paquetes
library(ggplot2)
library(tidyverse)
library(skimr)
library(e1071)
require(ppclust)
require(factoextra)
require(cluster)
require(fclust)
require(psych)
library(NbClust)
library(fpc)
library(corrplot)
library(Hmisc)
library(MASS)

# TERRENO A ----------------------------------------------------------------------
# Cargando data ---- 
data_lab <- read.csv("C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2/datos_lab_csv.csv", head = TRUE, sep = ";")
data_A <- data_lab[data_lab$Terreno == "A",]
data <- data_A[,-c(1:4)]
data


# Análisis exploratorio ---- 
psych::describe(data)

coef_var <- function(x){
  CV <- c()
  for(i in seq_along(x)){
    CV[i] <- sd(x[,i])/mean(x[,i]) * 100
    
  }
  names(CV) <- colnames(x)
  print(CV)
  
}
coef.var <- function(x){
  sd(x)/mean(x) * 100
}

CVs <- coef_var(data)

cormatrix <- cor(data)
rcorr(as.matrix(data), type = "pearson")
stargazer(cormatrix, type = "html", out = "C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2/corrA.html")

# Función para reducir decimales
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

# Gráficos
plot(data)
boxplot(scale(data), col = "orange", las = 2, main = "Terreno A", ylab= "Z score")

# Atípicos 
z <- list() 
for(i in seq_along(data)){
  z[[i]] <- (data[,i]-mean(data[,i]))/sd(data[,i])
}
names(z) <- colnames(data)
z

atipicos <- list()
for(i in 1:length(z)){
  atipicos[[i]] <- z[[i]][z[[i]] > 3 | z[[i]] < - 3]
  
}
names(atipicos) <- colnames(data)
atipicos

# Normalidad (Shapiro-Wilk)
apply(data, 2, shapiro.test)

#Transformación CE --
# Log
CE <- data[,"CE"]
log_CE <- log(CE)
shapiro.test(log_CE)

# Box Cox
bc<-boxcox(lm(CE ~ 1), lambda=seq(-2,2,l=100))
L<-with(bc,x[which.max(y)]);L	# Lambda = -1.27
CEtrans <- (CE^L-1)/ L
shapiro.test(CEtrans) 


# POST - INTERPOLACIÓN KRIGING EN ARCGIS -------------
# Cargar data 
datos_interpolados <- read.csv("C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2/Interpolados_A.csv", head = TRUE, sep = ";")
head(datos_interpolados)
class(datos_interpolados)
str(datos_interpolados)
length(datos_interpolados)

# Preparar data 
coordenadas <- datos_interpolados[,1:2] # creando objeto con las coordenadas 
coordenadas
data <- datos_interpolados[,3:length(datos_interpolados)] # creando objeto de trabajo -> Eliminar columnas de coordenadas (Este, Norte)
head(data)

# ANÁLISIS DE COMPONENTES PRINCIPALES - ACP ----
PCA <- prcomp(data, scale = TRUE) 
PCA
summary(PCA) # La Standard Deviation al cuadrado equivale a los Eigenvalues 
cumsum <- cumsum(PCA$sdev^2 / sum(PCA$sdev^2)) * 100
cumsum

# Elegir la cantidad de PC's que se usarán (puede ser los de Eigenvalue > 1 o los que logran explicar más del 70-80% de la variabilidad)
fviz_eig(PCA, addlabels = TRUE, ylim = c(0, 100))

eig.val <- get_eigenvalue(PCA)
eig.val #Eigenvalues > 1 significan que ese PC explica mas variabilidad que una de las variables originales 

write.table(eig.val, file = "C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2/PCA_A.csv", sep = ",")

# RELACIONES ENTRE VARIABLES Y PCs
var_pca <- get_pca_var(PCA)

# Correlación 
var_pca$cor 
corrplot(var_pca$cor[,1:3], is.corr = FALSE) # Gráfico de correlación 

write.table(var_pca$cor[,1:3], file = "C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2/PCAcor_A.csv", sep = ",")

# Contribución 
var_pca$contrib #Contribución porcentual de cada variable al PC
# Contributions of variables to PC1
fviz_contrib(PCA, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(PCA, choice = "var", axes = 2, top = 10)

# Calidad de Representación 
var_pca$cos2 # Calidad de representación de la variable por cada PC según Cos2
corrplot(var_pca$cos2, is.corr=FALSE) # Gráfico de Cos2

fviz_cos2(PCA, choice = "var", axes = 1:3) # Calidad de representación de cada variable entre los primeros 3 PC's

# Variable Correlation Plot 
fviz_pca_var(PCA, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE,
             ggtheme = theme_minimal()
)

# Añadir PC scores a la tabla de datos
str(PCA)
PCA$x
data_pca <- cbind(data, PC1 = PCA$x[, 1], PC2 = PCA$x[, 2], PC3 = PCA$x[,3]) #Ingresar solo las columnas de los PC's que se usarán
head(data_pca)
data_pca
length(data_pca)

#CORRELACIÓN ENTRE VARIABLES Y PCs escogidos
cor(data_pca)

#EXTRAER PC SCORES O VARIABLES DE INTERÉS
data_int <- data.frame(PC1 = data_pca$PC1, PC2 = data_pca$PC2, PC3 = data_pca$PC3)
head(data_int)


# ANÁLISIS CLUSTER ------------------------------------------------
#1. ESTANDARIZAR DATOS Y MUESTREO------------------------
# Estandarización
data_scale <- as.data.frame(scale(data_int))
head(data_scale)

# Muestreo
muestra <- sample_n(data_scale, 7000) #debe ser lo mas grande posible
muestra

#2. K-MEANS CLUSTERING 
#2.1. Método 1 - K-means
fitK <- kmeans(data_scale, 3, nstart = 25)
fitK
str(fitK)
plot(data_int, col=fitK$cluster)

#2.3. VALIDACIÓN DEL MODELO Y SELECCIÓN DE NÚMERO DE CLUSTERS "K" ---- 
#2.3.1. BetweenSS/TotSS y Elbow Method ------------------------------------
valoresk<-list()
for(i in 1:10){
  valoresk[[i]] <- kmeans(data_scale, i)
}
valoresk

betweenss_totss<-list()
for(i in 1:10){
  betweenss_totss[[i]]<-valoresk[[i]]$betweenss/valoresk[[i]]$totss
}
betweenss_totss

BetSS_TotSS<-c()
for(i in 1:10){
  BetSS_TotSS[i]<-valoresk[[i]]$betweenss/valoresk[[i]]$totss
}
BetSS_TotSS

within_ss <- list()
for(i in 1:10){
  within_ss[[i]] <- sum(valoresk[[i]]$withinss)
}
within_ss

plot(1:10, betweenss_totss, type="b", ylab="Variabilidad Explicada", xlab="Zonas", col = "blue", 
     main = "Variabilidad Explicada en distintas cantidades de Zonas", lwd = 5, cex =1.5, cex.axis = 1.3
     , cex.lab = 1.3, xaxp  = c(1, 10, 9))

plot(1:10, within_ss, type = "b", ylab = "Within_SS", xlab = "Clusters(k)")


#2.3.2. VARIANCE REDUCTION para "k" cantidad de clusters ---- 
# Creando función de Variance Reduction VR() para multiples variables ---- 
VR <- function(x, y){
  y <- as.factor(y)
  var_reduction <- list()
  ar_total <- length(y)
  clusters <- as.numeric(levels(y))
  ar_clusters <- c()
  ar_relativas <- c()
  for(i in as.numeric(levels(y))){
    ar_clusters[i] <- length(x[y == i, ncol(x)])
    ar_relativas[i] <- ar_clusters[i]/length(y)
  }
  resultados_prop <- data.frame(Proportion = ar_relativas, 
                                Count = ar_clusters, 
                                Total = rep(ar_total, length(ar_clusters)), 
                                Cluster = clusters)
  
  vars_clusters <- list()
  for(i in seq_along(x[,-length(x)])){
    vars_clusters[[i]] <- tapply(X = x[,i], INDEX = y, FUN = var)
  }
  var_total <- apply(x, 2, var)
  
  var_reduction <- list()
  for(i in seq_along(x[,-length(x)])){
    var_reduction[[i]] <- (1-sum(resultados_prop$Proportion*vars_clusters[[i]])/var_total[i])*100 
  }
  var_reduction
} #Argumentos: x = data frame con variables y clusters, y = columna de clusters
# Creando función de Variance Reduction VR_k() para multiples variables en "k" clusters---- 
VR_k <- function(x, k){ # Argumentos: x = datos de interés crudos, k = rango de clusters (ej: 2:5)
  cl <- list()
  vr <- list()
  datos_clusters <- list()
  for(i in k){
    cl[[i]] <- kmeans(scale(x), i)
    datos_clusters[[i]] <- cbind(x, Cluster = cl[[i]]$cluster)
    vr[[i]] <- VR(datos_clusters[[i]], datos_clusters[[i]]$Cluster)
    
  }
  print(vr)
}
# Creando función de Variance Reduction VR_vector() para 1 variable ---- 
VR_vector <- function(x, y){ # Argumentos: x = data frame de variables de interés y clusters, y = columna de clusters
  y <- as.factor(y)
  var_reduction <- list()
  ar_total <- length(y)
  clusters <- as.numeric(levels(y))
  ar_clusters <- c()
  ar_relativas <- c()
  for(i in as.numeric(levels(y))){
    ar_clusters[i] <- length(x[y == i, ncol(x)])
    ar_relativas[i] <- ar_clusters[i]/length(y)
  }
  resultados_prop <- data.frame(Proportion = ar_relativas, 
                                Count = ar_clusters, 
                                Total = rep(ar_total, length(ar_clusters)), 
                                Cluster = clusters)
  
  vars_clusters <- tapply(X = x[,1], INDEX = y, FUN = var)
  var_total <- apply(x[1], 2, var)
  
  
  var_reduction <- (1-sum(resultados_prop$Proportion*vars_clusters)/var_total)*100 
  var_reduction
}
# Creando función de Variance Reduction VR_vector_k() para 1 variable en "k" clusters  ---- 
VR_vector_k <- function(x, k){ # Argumentos: x = data de interés, k = rango de clusters (ej: 2:5)
  cl <- list()
  vr <- c()
  datos_clusters <- list()
  for(i in k){
    cl[[i]] <- kmeans(scale(x), i)
    datos_clusters[[i]] <- as.data.frame(cbind(x, Cluster = as.factor(cl[[i]]$cluster)))
    vr[i] <- VR_vector(datos_clusters[[i]], datos_clusters[[i]]$Cluster)
    
  }
  print(vr)
} 
#---- 
var_reduc <- VR_k(data_int, 2:6)


#2.3.3 SILHOUETTE INDEX  ---- 
# Silhouette Index para distintos valores de k
km <- list()
sil_info <- list()
indices_sil <- c()
for(i in 2:6){
  km[[i]] <- eclust(muestra, "kmeans", k = i, nstart = 25, graph = FALSE) 
  sil_info[[i]] <- km[[i]]$silinfo
  indices_sil[i] <- sil_info[[i]]$avg.width
}
indices_sil

# Gráficos Silhouette plot 
for(i in 2:6){
  print(fviz_silhouette(km[[i]]))
}

# Gráfico Elbow Method para Silhouette Index 
fviz_nbclust(muestra, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

#2.3.4. DUNN INDEX ---- 
km_stats <- list()
km <- list()
indices_dunn <- c()
for(i in 2:6){
  km[[i]] <- eclust(muestra, "kmeans", k = i, nstart = 25, graph = FALSE) 
  km_stats[[i]] <- cluster.stats(dist(muestra), km[[i]]$cluster)
  indices_dunn[i] <- km_stats[[i]]$dunn
}
indices_dunn 

km_stats #todas las estadísticas de cada modelo de k-means clustering

plot(1:6, indices_dunn, type="b", ylab = "Dunn Index", xlab="Clusters(k)")

#2.3.5. GAP STATISTIC 
fviz_nbclust(muestra, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

#2.3.4. MÉTODO DE LA MAYORÍA (más de 30 indices para determinar numero óptimo de clusters) ---- 
#Ojo: Necesita bastante tiempo y GB de memoria ram para ejecutarse, pero entrega directamente el "k" óptimo
nb <- NbClust(muestra, distance = "euclidean", min.nc = 2,
              max.nc = 6, method = "kmeans", index ="all")

fviz_nbclust(nb) + theme_minimal()

# 2.4. RESUMEN DE INDICADORES PARA DISTINTOS VALORES DE "K" ---- 
indicadores_k <- data.frame(Clusters = c(1:6),
                            Bss_Tss = BetSS_TotSS[1:6], 
                            Dunn = indices_dunn, 
                            Silhouette = indices_sil)
indicadores_k

#3. FUZZY C MEANS CLUSTERING ---- 
# 3.3.1. FUZZINES PERFORMANCE INDEX (FPI) y NORMALIZED CLASSIFICATION ENTROPY (NCE) ---- INDICADORES CLAVE
FPI <- function(cmem){
  c <- ncol(cmem)
  n <- nrow(cmem)
  
  1 - (c / (c - 1)) * (1 - sum(cmem^2) / n)
}
NCE <- function(cmem){
  c <- ncol(cmem)
  n <- nrow(cmem)
  
  (n / (n - c)) * (- sum(cmem * log(cmem)) / n)
}

# prepare variables
cl <- list()
fpi <- nce <- NULL

# cycle through the desired number of clusters
for(i in 2:6){
  cl[[i]] <- cmeans(data_scale, i, 20, method = "cmeans")
  fpi <- c(fpi, FPI(cl[[i]]$membership))
  nce <- c(nce, NCE(cl[[i]]$membership))
}

# add space for the second axis label
par(mar = c(5,4,1,4) + .1)

# plot FPI
plot(2:6, fpi, lty = 2, pch = 18, type = "b", xlab = "Número de zonas", ylab = "FPI", col = "red")

# plot NCE, manually adding the second axis
par(new = TRUE)
plot(2:6, nce, lty = 1, pch = 15, type = "b", xlab = "", ylab = "", axes = FALSE, col = "blue")
axis(4, at = pretty(range(nce)))
mtext("NCE", side = 4, line = 3)

# add legend
legend("bottom", legend = c("FPI", "NCE"), pch = c(18,15), lty = c(2,1), horiz = TRUE, col = c("red", "blue"))


# 3.3.6. BetweenSS/TotSS y Elbow Method 
valoresc<-list()
for(i in 1:6){
  valoresc[[i]] <- fcm(muestra, i)
}
valoresc

betweenss_totss_c <- list()
for(i in 1:6){
  betweenss_totss_c[[i]]<-valoresc[[i]]$sumsqrs$between.ss/valoresc[[i]]$sumsqrs$tot.ss
}
betweenss_totss_c

BetSS_TotSS_c<-c()
for(i in 1:6){
  BetSS_TotSS_c[i]<-valoresc[[i]]$sumsqrs$between.ss/valoresc[[i]]$sumsqrs$tot.ss
}
BetSS_TotSS_c

within_ss_c <- list()
for(i in 1:6){
  within_ss_c[[i]] <- sum(valoresc[[i]]$sumsqrs$within.ss)
}
within_ss_c

plot(1:6, betweenss_totss_c, type="b", ylab="Between_SS/Tot_SS", xlab="Clusters(c)")

plot(1:6, within_ss_c, type = "b", ylab = "Within_SS", xlab = "Clusters(c)")

#3.4. RESUMEN DE INDICADORES CLAVE PARA DISTINTOS VALORES DE "C"
indicadores_c <- data.frame(Clusters = c(2:6),
                            Bss_Tss = BetSS_TotSS_c[2:6], 
                            FPI = fpi, 
                            NCE = nce)
indicadores_c
indicadores_k

# 4. SELECCIÓN DEL MODELO FINAL, DELIMITACIÓN DE ZONAS Y RESULTADOS ---- 
# Modelo final (escoger uno de los 2 y establecer valor de "i")
modelo <- kmeans(data_scale, i, nstart = 25)
modelo <- fcm(data_scale, 2, nstart = 5)

# Añadir clusters a la tabla de datos 
zonas <- as.factor(modelo$cluster)
datos_zonas <- cbind(datos_interpolados, Zona = zonas)
datos_zonas <- datos_zonas %>% arrange(Zona)
datos_zonas %>% head

# Mapa de Zonas de Manejo 
ggplot(datos_zonas, aes(x = Este, y = Norte, col = Zona)) + 
  geom_point()

# Biplot de Componentes Principales coloreados por Cluster
data_plot <- cbind(data_pca, Zona = datos_zonas$Zona)
colores <- c("#FF0000", "#FFFF00", "#00FF00")
ggplot(data_plot, aes(x = PC1, y = PC2, col = Zona)) + 
  geom_point() + 
  labs(x = paste("PC1 (", specify_decimal(cumsum[1], 2),"%)"), y = paste("PC2 (", specify_decimal(cumsum[2]-cumsum[1], 2),"%)")) + 
  ggtitle(label = "Gráfico de Clusters") + 
  scale_color_manual(values = colores) +
  theme(axis.line.y = element_line(),
        axis.line.x = element_line(), 
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 13, angle = 90, hjust = 0.5),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        legend.title= element_text(size= 15),
        legend.text= element_text(size= 13),
        legend.position = c(0.15, 0.2),
        legend.background = element_rect(fill="lightblue", 
                                         size=0.5, linetype="solid",
                                         colour ="black"))

# Exportar datos de zonas 
write.csv(datos_zonas, "C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2/datos_zonas_A.csv", row.names = FALSE)
datos_zonas <- read.csv("C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2/datos_zonas_A.csv", head = TRUE, sep = ",")
datos_zonas


# POST - SELECCIÓN DE PUNTOS DE MUESTREO EN CADA ZONA EN ARCGIS ------------------------------------------------------------------------ 
# Puntos de muestreo dentro de cada zona 
puntos_Z1 <- read.csv("C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2/puntos_A_zona1.csv", head = TRUE, sep = ",")
puntos_Z1$Zona <- c(rep(1, nrow(puntos_Z1))) %>% as.factor
puntos_Z1

puntos_Z2 <- read.csv("C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2/puntos_A_zona2.txt", head = TRUE, sep = ",")
puntos_Z2$Zona <- c(rep(2, nrow(puntos_Z2))) %>% as.factor
puntos_Z2

datos_zonas <- rbind(puntos_Z1, puntos_Z2)
datos_zonas <- datos_zonas[,6:length(datos_zonas)]
datos_zonas

# Análisis Exploratorio de Clusters/Zonas ----
#Estadísticos descriptivos 
resumen_zonas <- describeBy(datos_zonas[,3:(length(datos_zonas)-1)], datos_zonas$Zona)
resumen_zonas

psych::describe(datos_zonas)

# Coeficientes de variación 
data_A <- datos_zonas[,-c(1, 2, length(datos_zonas))]
CV_campo <- coef_var(data_A)
CV_campo_vector <- c()
for(i in 1:length(CV_campo)){
  CV_campo_vector[i] <- CV_campo[[i]]
}
CV_campo_vector

CV_zonas <- list()
for(i in 3:(length(datos_zonas)-1)){
  CV_zonas[[i]] <- tapply(datos_zonas[,i], datos_zonas[,length(datos_zonas)], coef.var)
}
names(CV_zonas) <- colnames(datos_zonas[,1:(length(datos_zonas)-1)])
CV_zonas

CV_Zona1 <- c()
for(i in 3:length(CV_zonas)){
  CV_Zona1[i] <- CV_zonas[[i]][1]
}
CV_Zona1


CV_Zona2 <- c()
for(i in 3:length(CV_zonas)){
  CV_Zona2[i] <- CV_zonas[[i]][2]
}
CV_Zona2

# Boxplots 
# Creando función  que entrega las medias y los boxplots de cada variable en las distintas zonas ----
mean_treatment_boxplot <- function(x, t, k = 1, j = ncol(x)){
  colores <- c("#FF0000", "#FFFF00", "#00FF00")
  medias <- list()
  boxplot <- list()
  for(i in k:j){
    medias[[i]] <- tapply(x[,i], t, mean)
    boxplot[[i]] <- ggplot(x, aes_string(as.factor(t), colnames(x[i]), fill = t)) + 
      geom_boxplot(show.legend = FALSE) + 
      stat_summary(fun = mean, geom = "point", shape = 15, color = "black", size = 2, show.legend = FALSE) + 
      labs(x = "Zona", y = colnames(x[i]), fill = "Zonas") + 
      ggtitle(label = paste("Distribución de ",colnames(x[i]), "en cada Zona")) +
      scale_fill_manual(values = colores) + 
      theme(plot.title = element_text(hjust = 0.5, size = 18),
            panel.border = element_rect(colour = "black", fill = NA),
            plot.background = element_rect(colour = "black"), 
            axis.text.y = element_text(size = 14, angle = 90, hjust = 0.5),
            axis.text.x = element_text(size = 15),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16))
  }
  names(medias) <- colnames(x)
  print(medias)
  print(boxplot)
}
#----
mean_treatment_boxplot(datos_zonas[,3:length(datos_zonas)-1], datos_zonas$Zona)

# Comparación entre Clusters/Zonas 
# ANOVA
anovas <- list()
resumenes <- list()
p_values <- c()
for(i in 3:(length(datos_zonas)-1)){
  anovas[[i]] <- aov(datos_zonas[,i] ~ Zona, datos_zonas) 
  resumenes[[i]] <- anovas[[i]] %>% summary
  p_values[i] <- resumenes[[i]][[1]]$'Pr(>F)'[1]
}
names(resumenes) <- colnames(datos_zonas[,1:(length(datos_zonas)-1)])
resumenes
names(p_values) <- colnames(datos_zonas[,1:(length(datos_zonas)-1)])
p_values

# Tukey 
tukey <- list()
for(i in 3:(length(datos_zonas)-1)){
  tukey[[i]] <- TukeyHSD(anovas[[i]])
}
names(tukey) <- colnames(datos_zonas[,1:(length(datos_zonas)-1)])
tukey

# Kruskal-Wallis 
kruskal <- list()
for(i in 3:(length(datos_zonas)-1)){
  kruskal[[i]] <- kruskal.test(datos_zonas[,i], datos_zonas$Zona)
}
names(kruskal) <- colnames(datos_zonas[,-length(datos_zonas)])
kruskal

# t-test 
t <- list()
t.p_values <- c()
for(i in 3:(length(datos_zonas)-1)){
  t[[i]] <- t.test(datos_zonas[datos_zonas$Zona == 1, i], datos_zonas[datos_zonas$Zona == 2, i])
  t.p_values[i] <- t[[i]]$p.value
}
names(t) <- colnames(datos_zonas[,-length(datos_zonas)])
names(t.p_values) <- colnames(datos_zonas[,-length(datos_zonas)])
t
t.p_values

# Mann-Whitney 
mann_whitney <- list()
for(i in 3:9){
  mann_whitney[[i]] <- wilcox.test(datos_zonas[datos_zonas$Zona == 1, i], datos_zonas[datos_zonas$Zona == 2, i], paired = FALSE)
}
names(mann_whitney) <- colnames(datos_zonas[,1:9])
mann_whitney


# Resumen de Zonas
resumen_final <- data.frame(Variable = colnames(data_A),
                            'Zona 1' = resumen_zonas$'1'$mean, 
                            'Zona 2' = resumen_zonas$'2'$mean,
                            CV_Zona1 = CV_Zona1[3:length(CV_Zona1)], 
                            CV_Zona2 = CV_Zona2[3:length(CV_Zona2)],
                            CV_campo = CV_campo_vector, 
                            t.test_p_value = t.p_values[3:length(p_values)])
resumen_final
write.table(resumen_final, "C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2//resumen_final_A.csv")



# TERRENO B ----------------------------------------------------------------------
# Preparando data ---- 
data_lab
data_lab$CE <- replace_na(data_lab$CE, 155.85)
data_B <- data_lab[data_lab$Terreno == "B",]
data <- data_B[,-c(1:4)]
data


# Análisis exploratorio ---- 
psych::describe(data)

CVs <- coef_var(data)

cormatrix <- cor(data)
rcorr(as.matrix(data), type = "pearson")
stargazer(cormatrix, type = "html", out = "C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2/corrB.html")

# Función para reducir decimales
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

# Gráficos
plot(data)
boxplot(scale(data), col = "orange", las = 2, main = "Terreno B", ylab= "Z score")

# Atípicos 
z <- list() 
for(i in seq_along(data)){
  z[[i]] <- (data[,i]-mean(data[,i]))/sd(data[,i])
}
names(z) <- colnames(data)
z

atipicos <- list()
for(i in 1:length(z)){
  atipicos[[i]] <- z[[i]][z[[i]] > 3 | z[[i]] < - 3]
  
}
names(atipicos) <- colnames(data)
atipicos

# Normalidad (Shapiro-Wilk)
apply(data, 2, shapiro.test)

#Transformación CE --
# Log
CE <- data[,"CE"]
log_CE <- log(CE)
shapiro.test(log_CE)

# Box Cox
bc<-boxcox(lm(CE ~ 1), lambda=seq(-2,2,l=100))
L<-with(bc,x[which.max(y)]);L	# Lambda = -1.76
CEtrans <- (CE^L-1)/ L
shapiro.test(CEtrans) 


# POST - INTERPOLACIÓN KRIGING EN ARCGIS -------------
# Cargar data 
datos_interpolados <- read.csv("C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2/Interpolados_B.csv", head = TRUE, sep = ";")
head(datos_interpolados)
class(datos_interpolados)
str(datos_interpolados)
length(datos_interpolados)

# Preparar data 
coordenadas <- datos_interpolados[,1:2] # creando objeto con las coordenadas 
coordenadas
data <- datos_interpolados[,3:length(datos_interpolados)] # creando objeto de trabajo -> Eliminar columnas de coordenadas (Este, Norte)
head(data)

# ANÁLISIS DE COMPONENTES PRINCIPALES - ACP ----
PCA <- prcomp(data, scale = TRUE) 
PCA
summary(PCA) # La Standard Deviation al cuadrado equivale a los Eigenvalues 
cumsum <- cumsum(PCA$sdev^2 / sum(PCA$sdev^2)) * 100
cumsum

# Elegir la cantidad de PC's que se usarán (puede ser los de Eigenvalue > 1 o los que logran explicar más del 70-80% de la variabilidad)
fviz_eig(PCA, addlabels = TRUE, ylim = c(0, 100))

eig.val <- get_eigenvalue(PCA)
eig.val #Eigenvalues > 1 significan que ese PC explica mas variabilidad que una de las variables originales 

write.table(eig.val, file = "C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2/PCA_B.csv", sep = ",")

# RELACIONES ENTRE VARIABLES Y PCs
var_pca <- get_pca_var(PCA)

# Correlación 
var_pca$cor 
corrplot(var_pca$cor[,1:3], is.corr = FALSE) # Gráfico de correlación 

write.table(var_pca$cor[,1:3], file = "C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2/PCAcor_B.csv", sep = ",")

# Contribución 
var_pca$contrib #Contribución porcentual de cada variable al PC
# Contributions of variables to PC1
fviz_contrib(PCA, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(PCA, choice = "var", axes = 2, top = 10)

# Calidad de Representación 
var_pca$cos2 # Calidad de representación de la variable por cada PC según Cos2
corrplot(var_pca$cos2, is.corr=FALSE) # Gráfico de Cos2

fviz_cos2(PCA, choice = "var", axes = 1:3) # Calidad de representación de cada variable entre los primeros 3 PC's

# Variable Correlation Plot 
fviz_pca_var(PCA, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE,
             ggtheme = theme_minimal()
)

# Añadir PC scores a la tabla de datos
str(PCA)
PCA$x
data_pca <- cbind(data, PC1 = PCA$x[, 1], PC2 = PCA$x[, 2], PC3 = PCA$x[,3]) #Ingresar solo las columnas de los PC's que se usarán
head(data_pca)
data_pca
length(data_pca)

#CORRELACIÓN ENTRE VARIABLES Y PCs escogidos
cor(data_pca)

#EXTRAER PC SCORES O VARIABLES DE INTERÉS
data_int <- data.frame(PC1 = data_pca$PC1, PC2 = data_pca$PC2, PC3 = data_pca$PC3)
head(data_int)


# ANÁLISIS CLUSTER ------------------------------------------------
#1. ESTANDARIZAR DATOS Y MUESTREO------------------------
# Estandarización
data_scale <- as.data.frame(scale(data_int))
head(data_scale)

# Muestreo
muestra <- sample_n(data_scale, 7000) #debe ser lo mas grande posible
muestra

#2. K-MEANS CLUSTERING 
#2.1. Método 1 - K-means
fitK <- kmeans(data_scale, 3, nstart = 25)
fitK
str(fitK)
plot(data_int, col=fitK$cluster)

#2.3. VALIDACIÓN DEL MODELO Y SELECCIÓN DE NÚMERO DE CLUSTERS "K" ---- 
#2.3.1. BetweenSS/TotSS y Elbow Method ------------------------------------
valoresk<-list()
for(i in 1:10){
  valoresk[[i]] <- kmeans(data_scale, i)
}
valoresk

betweenss_totss<-list()
for(i in 1:10){
  betweenss_totss[[i]]<-valoresk[[i]]$betweenss/valoresk[[i]]$totss
}
betweenss_totss

BetSS_TotSS<-c()
for(i in 1:10){
  BetSS_TotSS[i]<-valoresk[[i]]$betweenss/valoresk[[i]]$totss
}
BetSS_TotSS

within_ss <- list()
for(i in 1:10){
  within_ss[[i]] <- sum(valoresk[[i]]$withinss)
}
within_ss

plot(1:10, betweenss_totss, type="b", ylab="Variabilidad Explicada", xlab="Zonas", col = "blue", 
     main = "Variabilidad Explicada en distintas cantidades de Zonas", lwd = 5, cex =1.5, cex.axis = 1.3
     , cex.lab = 1.3, xaxp  = c(1, 10, 9))

plot(1:10, within_ss, type = "b", ylab = "Within_SS", xlab = "Clusters(k)")


#---- Variance Reduction
var_reduc <- VR_k(data_int, 2:6)


#2.3.3 SILHOUETTE INDEX  ---- 
# Silhouette Index para distintos valores de k
km <- list()
sil_info <- list()
indices_sil <- c()
for(i in 2:6){
  km[[i]] <- eclust(muestra, "kmeans", k = i, nstart = 25, graph = FALSE) 
  sil_info[[i]] <- km[[i]]$silinfo
  indices_sil[i] <- sil_info[[i]]$avg.width
}
indices_sil

# Gráficos Silhouette plot 
for(i in 2:6){
  print(fviz_silhouette(km[[i]]))
}

# Gráfico Elbow Method para Silhouette Index 
fviz_nbclust(muestra, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

#2.3.4. DUNN INDEX ---- 
km_stats <- list()
km <- list()
indices_dunn <- c()
for(i in 2:6){
  km[[i]] <- eclust(muestra, "kmeans", k = i, nstart = 25, graph = FALSE) 
  km_stats[[i]] <- cluster.stats(dist(muestra), km[[i]]$cluster)
  indices_dunn[i] <- km_stats[[i]]$dunn
}
indices_dunn 

km_stats #todas las estadísticas de cada modelo de k-means clustering

plot(1:6, indices_dunn, type="b", ylab = "Dunn Index", xlab="Clusters(k)")

#2.3.5. GAP STATISTIC 
fviz_nbclust(muestra, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

#2.3.4. MÉTODO DE LA MAYORÍA (más de 30 indices para determinar numero óptimo de clusters) ---- 
#Ojo: Necesita bastante tiempo y GB de memoria ram para ejecutarse, pero entrega directamente el "k" óptimo
nb <- NbClust(muestra, distance = "euclidean", min.nc = 2,
              max.nc = 6, method = "kmeans", index ="all")

fviz_nbclust(nb) + theme_minimal()

# 2.4. RESUMEN DE INDICADORES PARA DISTINTOS VALORES DE "K" ---- 
indicadores_k <- data.frame(Clusters = c(1:6),
                            Bss_Tss = BetSS_TotSS[1:6], 
                            Dunn = indices_dunn, 
                            Silhouette = indices_sil)
indicadores_k

#3. FUZZY C MEANS CLUSTERING ---- 
# 3.3.1. FUZZINES PERFORMANCE INDEX (FPI) y NORMALIZED CLASSIFICATION ENTROPY (NCE) ---- INDICADORES CLAVE
FPI <- function(cmem){
  c <- ncol(cmem)
  n <- nrow(cmem)
  
  1 - (c / (c - 1)) * (1 - sum(cmem^2) / n)
}
NCE <- function(cmem){
  c <- ncol(cmem)
  n <- nrow(cmem)
  
  (n / (n - c)) * (- sum(cmem * log(cmem)) / n)
}

# prepare variables
cl <- list()
fpi <- nce <- NULL

# cycle through the desired number of clusters
for(i in 2:6){
  cl[[i]] <- cmeans(data_scale, i, 20, method = "cmeans")
  fpi <- c(fpi, FPI(cl[[i]]$membership))
  nce <- c(nce, NCE(cl[[i]]$membership))
}

# add space for the second axis label
par(mar = c(5,4,1,4) + .1)

# plot FPI
plot(2:6, fpi, lty = 2, pch = 18, type = "b", xlab = "Número de zonas", ylab = "FPI", col = "red")

# plot NCE, manually adding the second axis
par(new = TRUE)
plot(2:6, nce, lty = 1, pch = 15, type = "b", xlab = "", ylab = "", axes = FALSE, col = "blue")
axis(4, at = pretty(range(nce)))
mtext("NCE", side = 4, line = 3)

# add legend
legend("bottom", legend = c("FPI", "NCE"), pch = c(18,15), lty = c(2,1), horiz = TRUE, col = c("red", "blue"))


# 3.3.6. BetweenSS/TotSS y Elbow Method 
valoresc<-list()
for(i in 1:6){
  valoresc[[i]] <- fcm(muestra, i)
}
valoresc

betweenss_totss_c <- list()
for(i in 1:6){
  betweenss_totss_c[[i]]<-valoresc[[i]]$sumsqrs$between.ss/valoresc[[i]]$sumsqrs$tot.ss
}
betweenss_totss_c

BetSS_TotSS_c<-c()
for(i in 1:6){
  BetSS_TotSS_c[i]<-valoresc[[i]]$sumsqrs$between.ss/valoresc[[i]]$sumsqrs$tot.ss
}
BetSS_TotSS_c

within_ss_c <- list()
for(i in 1:6){
  within_ss_c[[i]] <- sum(valoresc[[i]]$sumsqrs$within.ss)
}
within_ss_c

plot(1:6, betweenss_totss_c, type="b", ylab="Between_SS/Tot_SS", xlab="Clusters(c)")

plot(1:6, within_ss_c, type = "b", ylab = "Within_SS", xlab = "Clusters(c)")

#3.4. RESUMEN DE INDICADORES CLAVE PARA DISTINTOS VALORES DE "C"
indicadores_c <- data.frame(Clusters = c(2:6),
                            Bss_Tss = BetSS_TotSS_c[2:6], 
                            FPI = fpi, 
                            NCE = nce)
indicadores_c
indicadores_k

# 4. SELECCIÓN DEL MODELO FINAL, DELIMITACIÓN DE ZONAS Y RESULTADOS ---- 
# Modelo final (escoger uno de los 2 y establecer valor de "i")
modelo <- kmeans(data_scale, i, nstart = 25)
modelo <- fcm(data_scale, 2, nstart = 5)

# Añadir clusters a la tabla de datos 
zonas <- as.factor(modelo$cluster)
datos_zonas <- cbind(datos_interpolados, Zona = zonas)
datos_zonas <- datos_zonas %>% arrange(Zona)
datos_zonas %>% head

# Mapa de Zonas de Manejo 
ggplot(datos_zonas, aes(x = Este, y = Norte, col = Zona)) + 
  geom_point()

# Biplot de Componentes Principales coloreados por Cluster
data_plot <- cbind(data_pca, Zona = datos_zonas$Zona)
colores <- c("#FF0000", "#FFFF00", "#00FF00")
ggplot(data_plot, aes(x = PC1, y = PC2, col = Zona)) + 
  geom_point() + 
  labs(x = paste("PC1 (", specify_decimal(cumsum[1], 2),"%)"), y = paste("PC2 (", specify_decimal(cumsum[2]-cumsum[1], 2),"%)")) + 
  ggtitle(label = "Gráfico de Clusters") + 
  scale_color_manual(values = colores) +
  theme(axis.line.y = element_line(),
        axis.line.x = element_line(), 
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 13, angle = 90, hjust = 0.5),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), 
        legend.title= element_text(size= 15),
        legend.text= element_text(size= 13),
        legend.position = c(0.15, 0.2),
        legend.background = element_rect(fill="lightblue", 
                                         size=0.5, linetype="solid",
                                         colour ="black"))

# POST - SELECCIÓN DE PUNTOS DE MUESTREO EN CADA ZONA EN ARCGIS ------------------------------------------------------------------------ 
# Puntos de muestreo dentro de cada zona 
puntos_Z1 <- read.csv("C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2/puntos_B_zona1.csv", head = TRUE, sep = ";")
puntos_Z1$Zona <- c(rep(1, nrow(puntos_Z1))) %>% as.factor
puntos_Z1

puntos_Z2 <- read.csv("C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2/puntos_B_zona2.csv", head = TRUE, sep = ";")
puntos_Z2$Zona <- c(rep(2, nrow(puntos_Z2))) %>% as.factor
puntos_Z2

datos_zonas <- rbind(puntos_Z1, puntos_Z2)
datos_zonas <- datos_zonas[,6:length(datos_zonas)]
datos_zonas

# Análisis Exploratorio de Clusters/Zonas ----
#Estadísticos descriptivos 
resumen_zonas <- describeBy(datos_zonas[,3:(length(datos_zonas)-1)], datos_zonas$Zona)
resumen_zonas

# Coeficientes de variación 
data_B <- datos_zonas[,-c(1, 2, length(datos_zonas))]
data_B
CV_campo <- coef_var(data_B)
CV_campo_vector <- c()
for(i in 1:length(CV_campo)){
  CV_campo_vector[i] <- CV_campo[[i]]
}
CV_campo_vector

CV_zonas <- list()
for(i in 3:(length(datos_zonas)-1)){
  CV_zonas[[i]] <- tapply(datos_zonas[,i], datos_zonas[,length(datos_zonas)], coef.var)
}
CV_zonas

CV_Zona1 <- c()
for(i in 3:length(CV_zonas)){
  CV_Zona1[i] <- CV_zonas[[i]][1]
}
CV_Zona1


CV_Zona2 <- c()
for(i in 3:length(CV_zonas)){
  CV_Zona2[i] <- CV_zonas[[i]][2]
}
CV_Zona2

# Boxplots 
# Creando función  que entrega las medias y los boxplots de cada variable en las distintas zonas ----
mean_treatment_boxplot <- function(x, t, k = 1, j = ncol(x)){
  colores <- c("#FF0000", "#FFFF00", "#00FF00")
  medias <- list()
  boxplot <- list()
  for(i in k:j){
    medias[[i]] <- tapply(x[,i], t, mean)
    boxplot[[i]] <- ggplot(x, aes_string(as.factor(t), colnames(x[i]), fill = t)) + 
      geom_boxplot(show.legend = FALSE) + 
      stat_summary(fun = mean, geom = "point", shape = 15, color = "black", size = 2, show.legend = FALSE) + 
      labs(x = "Zona", y = colnames(x[i]), fill = "Zonas") + 
      ggtitle(label = paste("Distribución de ",colnames(x[i]), "en cada Zona")) +
      scale_fill_manual(values = colores) + 
      theme(plot.title = element_text(hjust = 0.5, size = 18),
            panel.border = element_rect(colour = "black", fill = NA),
            plot.background = element_rect(colour = "black"), 
            axis.text.y = element_text(size = 14, angle = 90, hjust = 0.5),
            axis.text.x = element_text(size = 15),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16))
  }
  names(medias) <- colnames(x)
  names(boxplot) <- colnames(x)
  print(medias)
  print(boxplot)
}
#----
mean_treatment_boxplot(datos_zonas[,3:length(datos_zonas)-1], datos_zonas$Zona)

# ANOVA entre Clusters/Zonas 
anovas <- list()
resumenes <- list()
p_values <- c()
for(i in 3:(length(datos_zonas)-1)){
  anovas[[i]] <- aov(datos_zonas[,i] ~ Zona, datos_zonas) 
  resumenes[[i]] <- anovas[[i]] %>% summary
  p_values[i] <- resumenes[[i]][[1]]$'Pr(>F)'[1]
}
names(resumenes) <- colnames(datos_zonas[,1:(length(datos_zonas)-1)])
resumenes
names(p_values) <- colnames(datos_zonas[,1:(length(datos_zonas)-1)])
p_values

# Tukey 
tukey <- list()
for(i in 3:(length(datos_zonas)-1)){
  tukey[[i]] <- TukeyHSD(anovas[[i]])
}
names(tukey) <- colnames(datos_zonas[,1:(length(datos_zonas)-1)])
tukey

# Kruskal-Wallis 
kruskal <- list()
for(i in 3:(length(datos_zonas)-1)){
  kruskal[[i]] <- kruskal.test(datos_zonas[,i], datos_zonas$Zona)
}
names(kruskal) <- colnames(datos_zonas[,-length(datos_zonas)])
kruskal

# t-test 
t <- list()
t.p_values <- c()
for(i in 3:(length(datos_zonas)-1)){
  t[[i]] <- t.test(datos_zonas[datos_zonas$Zona == 1, i], datos_zonas[datos_zonas$Zona == 2, i])
  t.p_values[i] <- t[[i]]$p.value
}
names(t) <- colnames(datos_zonas[,-length(datos_zonas)])
names(t.p_values) <- colnames(datos_zonas[,-length(datos_zonas)])
t
t.p_values

# Mann-Whitney 
mann_whitney <- list()
for(i in 3:9){
  mann_whitney[[i]] <- wilcox.test(datos_zonas[datos_zonas$Zona == 1, i], datos_zonas[datos_zonas$Zona == 2, i], paired = FALSE)
}
names(mann_whitney) <- colnames(datos_zonas[,1:9])
mann_whitney

# Resumen de Zonas
resumen_final <- data.frame(Variable = colnames(data_B),
                            'Zona 1' = resumen_zonas$'1'$mean, 
                            'Zona 2' = resumen_zonas$'2'$mean,
                            CV_Zona1 = CV_Zona1[3:length(CV_Zona1)], 
                            CV_Zona2 = CV_Zona2[3:length(CV_Zona2)],
                            CV_campo = CV_campo_vector, 
                            t.test_p_value = t.p_values[3:length(p_values)])

resumen_final
write.table(resumen_final, "C:/Users/User/Documents/INGENIERIA AGROFORESTAL/Campo 2/resumen_final_B.csv")
