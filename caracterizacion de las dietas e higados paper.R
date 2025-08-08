# Ejecutar la recolección de basura para liberar memoria
gc()
# Eliminar todas las variables del entorno de trabajo
rm(list = ls())
# Establecer el directorio de trabajo en el escritorio (opcional)
setwd("~/Desktop/R/Doctorad")
# Verificar que se haya establecido correctamente
getwd()

# Cargar librerías
library(dplyr)
library(readxl)
library(ggplot2)
library(tidyr)


# Leer los datos
planilla <- read_excel("planilla.xlsx", sheet = "mg_g")

# Crear un nuevo dataframe df con los datos
df <- planilla

# 1) Primero extraemos los nombres de las columnas de n-3 y n-6
n3_cols <- grep("n3", names(df), value = TRUE)
n6_cols <- grep("n6", names(df), value = TRUE)


# 2) Luego creamos un nuevo data frame con las dos variables total_n3 y total_n6
df2 <- df %>%
  mutate(
    total_n3 = rowSums(across(all_of(n3_cols)), na.rm = TRUE),
    total_n6 = rowSums(across(all_of(n6_cols)), na.rm = TRUE)
  )

df2$n6_n3 <- df2$total_n6 / df2$total_n3

# Crear el dataframe original 'dietas' filtrando los elementos que contienen "DIETA" en la columna 'MUESTRA'
dietas <- df2 %>%
  filter(grepl("DIETA", MUESTRA))

# Seleccionar columnas sin NA y sin ceros
dietas <- dietas %>%
  dplyr::select(where(~ !all(is.na(.)) & !all(. == 0, na.rm = TRUE))) %>%
  dplyr::select(where(~ mean(is.na(.) | . == 0, na.rm = TRUE) <= 0.60))

# Excluir columnas 1, 3, 4
dietas_pca <- dietas[, -c(1, 3, 4)]

# Eliminar columnas no numéricas
dietas_pca <- dietas_pca %>%
  select_if(is.numeric)

# PCA
pca_result <- prcomp(dietas_pca, scale. = TRUE)

# Visualizar resumen del PCA
summary(pca_result)

# Calcular el porcentaje de variabilidad explicada
var_explained <- summary(pca_result)$importance[2, ] * 100  # El segundo renglón contiene la variabilidad explicada

# Graficar PCA usando DIETA como etiquetas
pca_data <- as.data.frame(pca_result$x)
pca_data$DIETA <- dietas$DIETA  # Añadir la columna DIETA

# Especificar el orden deseado de las dietas
orden_dietas1 <- c("FO", "SH", "SHCA", "CA", "COCA", "CO")

ggplot(pca_data, aes(x = PC1, y = PC2, color = DIETA)) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c("FO" = "blue", "SH" = "orange", "SHCA" = "#FFA07A", 
                                "CA" = "darkgreen", "COCA" = "green", "CO" = "lightgreen"), 
                     breaks = orden_dietas1) +
  labs(title = "PCA Dietas", x = paste0("PC1 (", round(var_explained[1], 2), "%)"), 
       y = paste0("PC2 (", round(var_explained[2], 2), "%)"), color = "") +
  theme_light() +
  theme(legend.title = element_blank())

# Filtrar los datos de hígado
higado <- df2 %>%
  filter(grepl("HIGADO", MUESTRA))
higado <- higado %>%
  filter(DIETA != "INI")

# Seleccionar columnas sin NA y sin ceros
higado <- higado %>%
  dplyr::select(where(~ !all(is.na(.)) & !all(. == 0, na.rm = TRUE))) %>%
  dplyr::select(where(~ mean(is.na(.) | . == 0, na.rm = TRUE) <= 0.60))

# Excluir columnas 1, 3, 4
higado_pca <- higado[, -c(1, 3, 4)]

# Eliminar columnas no numéricas
higado_pca <- higado_pca %>%
  select_if(is.numeric)

# PCA para los datos de hígado
pca_higado <- prcomp(higado_pca, scale. = TRUE)

# Obtener los resultados del PCA
pca_higado_data <- as.data.frame(pca_higado$x)

pca_higado_data$MUESTRA <- higado$MUESTRA  # Añadir la columna MUESTRA

pca_higado_data$DIETA <- higado$DIETA  # Añadir la columna DIETA

# Añadir columnas vacías a pca_data para que coincidan con pca_higado_data
pca_data[, paste0("PC", 3:18)] <- NA

# Añadir la columna MUESTRA a pca_data
pca_data$MUESTRA <- dietas$MUESTRA

# Combinar los dataframes pca_data y pca_higado_data
pca_combined_data <- rbind(
  cbind(pca_data, Tipo = "Dieta"),
  cbind(pca_higado_data, Tipo = "Higado")
)

# Graficar PCA para dietas y hígado en el mismo plot, con leyenda en inglés
ggplot(pca_combined_data, aes(x = PC1, y = PC2, color = DIETA, shape = Tipo)) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c("FO" = "blue", "SH" = "orange", "SHCA" = "#FFA07A", 
                                "CA" = "darkgreen", "COCA" = "green", "CO" = "lightgreen"), 
                     breaks = orden_dietas1) +
  scale_shape_manual(values = c("Dieta" = 15, "Higado" = 17), 
                     labels = c("Dieta" = "Diet", "Higado" = "Liver")) +  # Traducir las etiquetas en la leyenda
  labs(title = "PCA of Fatty Acid Profiles of Diets and Liver", x = paste0("PC1 (", round(var_explained[1], 2), "%)"), 
       y = paste0("PC2 (", round(var_explained[2], 2), "%)"), color = "") +
  theme_light() +
  theme(legend.title = element_blank()) +
  guides(shape = guide_legend(override.aes = list(shape = c(15, 17)))) + # Cambiar la forma de la leyenda según el tipo
  guides(color = guide_legend(override.aes = list(shape = 17)))  # Mantener la forma de círculo para la leyenda de color


# Mapeo de los nombres de los ácidos grasos en español a inglés
label_mapping <- c(
  "C14:0 Acido Mirístico" = "Myristic Acid",
  "C15:0 Acido Pentadecanoico" = "Pentadecanoic Acid",
  "C16:0 Acido Palmítico" = "Palmitic Acid",
  "C17:0 Acido Heptadecanoico" = "Heptadecanoic Acid",
  "C18:0 Acido Esteárico" = "Stearic Acid",
  "C21:0 Acido Heneicosanoico" = "Heneicosanoic Acid",
  "C22:0 Acido Behénico" = "Behenic Acid",
  "C23:0 Acido Tricosanoico" = "Tricosanoic Acid",
  "Total de SAFAs" = "Total SFAs",
  "C16:1 Acido Palmitoléico" = "Palmitoleic Acid",
  "C17:1 Acido Cis-10-Heptadecenoico" = "Heptadecenoic Acid",
  "C18:1n9t Acido Elaidico" = "Elaidic Acid",
  "C18:1n9c Acido Oléico" = "Oleic Acid",
  "C20:1n9 Acido Cis-11-Eicosenoico" = "Eicosenoic Acid",
  "Total de MUFAs" = "Total MUFAs",
  "C18:2n6c Acido Linoléico" = "Linoleic Acid",
  "C18:3n3 Acido Linolénico" = "Linolenic Acid (n3)",
  "C18:3n6 Acido g-linolénico" = "Gamma Linolenic Acid",
  "C20:2 Acido Cis-11,14-Eicosadienoico" = "Eicosadienoic Acid",
  "C20:4n6 Acido Araquidónico" = "Arachidonic Acid",
  "C20:5n3 Acido Cis-5,8,11,14,17-Eicosapentaenoico" = "Eicosapentaenoic Acid",
  "C22:2 Acido Cis-13,16-Docosadienoico" = "Docosadienoic Acid",
  "C22:6n3 Acido Cis-4,7,10,13,16,19-Docosahexaenoico" = "Docosahexaenoic Acid",
  "Total de PUFAs" = "Total PUFAs",
  "total_n3" = "Total n3",
  "total_n6" = "Total n6",
  "n6_n3" = "n6:n3 Ratio",
  "% Lipidos totales" = "Total Lipids"
)

# Extraer las cargas de los componentes principales desde pca_result
pca_loadings <- pca_result$rotation  # Las cargas de las variables (componentes)

# Obtener las variables (ácidos grasos) y los coeficientes de PC1 y PC2
variables <- rownames(pca_loadings)
coeff_PC1 <- pca_loadings[, 1]
coeff_PC2 <- pca_loadings[, 2]

# Traducir los nombres de las variables a inglés usando el mapeo
variables_in_english <- label_mapping[variables]

# Crear las ecuaciones con los nombres de las variables en inglés y sus coeficientes
equation_PC1 <- paste("PC1 = ", paste(variables_in_english, "(", round(coeff_PC1, 2), ")", collapse = " + "), sep = "")
equation_PC2 <- paste("PC2 = ", paste(variables_in_english, "(", round(coeff_PC2, 2), ")", collapse = " + "), sep = "")

# Imprimir las ecuaciones
print(equation_PC1)
print(equation_PC2)


# Graficar PCA con las ecuaciones de los componentes fuera del gráfico
ggplot(pca_combined_data, aes(x = PC1, y = PC2, color = DIETA, shape = Tipo)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c("FO" = "blue", "SH" = "orange", "SHCA" = "#FFA07A", 
                                "CA" = "darkgreen", "COCA" = "green", "CO" = "lightgreen"), 
                     breaks = orden_dietas1) +
  scale_shape_manual(values = c("Dieta" = 15, "Higado" = 17), 
                     labels = c("Dieta" = "Diet", "Higado" = "Liver")) +  # Traducir las etiquetas en la leyenda
  labs(title = "PCA of Fatty Acid Profiles of Diets and Liver", 
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"), 
       y = paste0("PC2 (", round(var_explained[2], 2), "%)"), 
       color = "") +
  theme_light() +
  theme(legend.title = element_blank()) +
  guides(shape = guide_legend(override.aes = list(shape = c(15, 17)))) + 
  guides(color = guide_legend(override.aes = list(shape = 17)))

higado1 <- df2 %>%
  filter(grepl("HIGADO", MUESTRA))
# Excluir columnas 1, 3, 4
higado1_pca <- higado1[, -c(1, 3, 4)]

# Seleccionar solo las columnas numéricas para el PCA
higado1_pca <- higado1_pca %>%
  dplyr::select(where(is.numeric))  # Selecciona solo columnas numéricas

# Realizar el PCA
pca_result <- prcomp(higado1_pca, scale. = TRUE)

# Resumen del PCA
summary(pca_result)

# Obtener el porcentaje de variabilidad explicada por cada componente principal
var_explained <- summary(pca_result)$importance[2, ] * 100

# Añadir las componentes principales al dataframe para visualización
pca_data <- as.data.frame(pca_result$x)

# Añadir la columna DIETA al dataframe pca_data para poder diferenciar las dietas en el gráfico
pca_data$DIETA <- higado1$DIETA  # Añadir la columna DIETA a los datos PCA

# Calcular el porcentaje de variabilidad explicada por cada componente
var_explained <- summary(pca_result)$importance[2, ] * 100  # El segundo renglón contiene la variabilidad explicada

# Visualización del PCA utilizando ggplot2
ggplot(pca_data, aes(x = PC1, y = PC2, color = DIETA)) +
  geom_point(size = 3) +  # Graficar los puntos
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +  # Línea horizontal en cero
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +  # Línea vertical en cero
  scale_color_manual(values = c("FO" = "blue", "SH" = "#FFa11f", "SHCA" = "#FFd17A", 
                                "CA" = "darkgreen", "COCA" = "green", "CO" = "lightgreen")) +  # Colores para las dietas
  labs(title = "PCA of Liver Data Including INI Condition", 
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"), 
       y = paste0("PC2 (", round(var_explained[2], 2), "%)"), 
       color = "Diet Condition") +  # Etiquetas y título
  theme_minimal() +  # Estilo de tema minimalista
  theme(legend.title = element_blank())  # Eliminar título de la leyenda


# 1. Realizar clustering k-means sobre los primeros dos componentes principales (PC1, PC2)
set.seed(123)  # Fijar semilla para reproducibilidad
k_max <- 10  # Definir un máximo de k (número de clústeres)
wss <- numeric(k_max)

# Calcular el WSS para diferentes valores de k
for (k in 1:k_max) {
  kmeans_result <- kmeans(pca_data[, c("PC1", "PC2")], centers = k, nstart = 25)
  wss[k] <- sum(kmeans_result$withinss)
}

# Graficar el método del codo
ggplot(data.frame(k = 1:k_max, WSS = wss), aes(x = k, y = WSS)) +
  geom_line() +
  geom_point() +
  labs(title = "Elbow Method for Determining Optimal Number of Clusters", 
       x = "Number of Clusters", y = "Within-Cluster Sum of Squares (WSS)") +
  theme_minimal()

# A partir del gráfico, puedes decidir el número de clústeres donde el WSS disminuye rápidamente (codo)
# Supongamos que elegimos k = 3 como el número de clústeres óptimo
# Especificar el orden de las dietas en la leyenda
orden_dietas1 <- c("FO", "SH", "SHCA", "CA", "COCA", "CO", "INI")

# Convertir la columna DIETA a factor y especificar el orden
pca_data$DIETA <- factor(pca_data$DIETA, levels = orden_dietas1)

# 2. Realizar el clustering k-means con el número de clústeres óptimo (k=4)
kmeans_result <- kmeans(pca_data[, c("PC1", "PC2")], centers = 3, nstart = 25)

# 3. Añadir la variable de clúster al dataframe PCA para visualizarlo
pca_data$cluster <- factor(kmeans_result$cluster)
# Eliminar los NA en las columnas PC1 y PC2 si existen
pca_data <- pca_data %>% drop_na(PC1, PC2)

# Asegúrate de que 'cluster' sea un factor
pca_data$cluster <- as.factor(pca_data$cluster)

# Graficar PCA con elipses de agrupamiento, ajustando nivel de confianza a 0.95
ggplot(pca_data, aes(x = PC1, y = PC2, color = DIETA)) +
  geom_point(size = 4) +  # Aumentar el tamaño de los puntos
  stat_ellipse(aes(group = cluster), level = 0.95, linetype = "dotted") +  # Elipses punteadas
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c("FO" = "blue", "SH" = "orange", "SHCA" = "#FFA07A", 
                                "CA" = "darkgreen", "COCA" = "green", "CO" = "lightgreen", 
                                "INI" = "grey")) +  # Colores para las dietas y gris para INI
  labs(title = "Clusters based on k-means (PCA)", 
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"), 
       y = paste0("PC2 (", round(var_explained[2], 2), "%)"), 
       color = "Diet") +
  theme_light() +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)))  # Ajustar tamaño de los puntos en la leyenda
