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
library(reshape2)
library(car)  # For ANOVA
library(ggpubr)  # For multiple comparison tests
library(factoextra)  # For PCA
library(ggrepel)
library(agricolae)
library(rlang)
library(stringr)
library(multcomp)
library(purrr)
library(broom)
library(FSA)
library(ggsignif)
library(openxlsx)
library(ggfortify)
library(vegan)
library(ggpubr)
library(gridExtra)

# Definir los colores para cada tratamiento
treatment_colors <- c("FO" = "blue", 
                      "CO" = "#98FB98", 
                      "SH" = "#FFA07A", 
                      "CA" = "darkgreen", 
                      "SHCA" = "#FFCC99",
                      "COCA" = "#32CD32"
)

# Asegúrate de que los tratamientos están en el orden correcto
treatment_order <- c("FO", "SH", "SHCA", "COCA", "CO", "CA")

# Leer los datos iniciales
planilla <- read_excel("planilla.xlsx", sheet = "mg_g")
tabla_genes_R <- read_excel("~/Desktop/UBAINT Doctoral/Ensayo/tabla genes R.xlsx")
resultados_microbiota <- read_excel("resultados microbiota.xlsx", 
                                    sheet = "Track")
microbiota_abundancia <- read_excel("resultados microbiota.xlsx", 
                                    sheet = "Tabla S9")
microbiota_gen <- read_excel("resultados microbiota.xlsx", 
                             sheet = "Tabla S13")


# Tratando los datos de microbiota

str(resultados_microbiota)
summary(resultados_microbiota)
str(microbiota_abundancia)
summary(microbiota_abundancia)

# Filtrar datos experimentales (eliminando controles internos)
resultados_microbiota_filtrado <- resultados_microbiota %>%
  filter(!grepl("Control", SampleID, ignore.case = TRUE))

# Verificar que el filtrado fue correcto
table(resultados_microbiota_filtrado$SampleID)

# Agregar proporciones al dataframe filtrado
resultados_microbiota_filtrado <- resultados_microbiota_filtrado %>%
  mutate(
    Prop_Filtered = filtered / input,
    Prop_NonChim = nonchim / input
  )

# Inspeccionar los nuevos cálculos
head(resultados_microbiota_filtrado)

# Resumir datos por tratamiento
resumen_microbiota <- resultados_microbiota_filtrado %>%
  group_by(SampleID) %>%
  summarize(
    Mean_Input = mean(input, na.rm = TRUE),
    Mean_Filtered = mean(filtered, na.rm = TRUE),
    Mean_NonChim = mean(nonchim, na.rm = TRUE),
    Prop_Filtered = mean(Prop_Filtered, na.rm = TRUE),
    Prop_NonChim = mean(Prop_NonChim, na.rm = TRUE),
    .groups = "drop"
  )

# Ver el resumen
print(resumen_microbiota)

# Ensure SampleID respects the desired order
resultados_microbiota_filtrado$SampleID <- factor(resultados_microbiota_filtrado$SampleID, 
                                                  levels = treatment_order)

# Generate the boxplot with ggplot2
ggplot(resultados_microbiota_filtrado, aes(x = SampleID, y = Prop_NonChim, fill = SampleID)) +
  geom_boxplot() +
  labs(
    title = "Proportion of Nonchim Sequences by Treatment",
    x = "",
    y = "Proportion of Nonchim Sequences"
  ) +
  scale_fill_manual(values = treatment_colors) +  # Apply custom colors
  theme_minimal() +
  theme(
    legend.position = "right",  # Adjust legend position
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1)  # Rotate x-axis labels
  )

# Ensure SampleID respects the desired order
resultados_microbiota_filtrado$SampleID <- factor(resultados_microbiota_filtrado$SampleID, 
                                                  levels = treatment_order)

# Generate the boxplot with ggplot2
ggplot(resultados_microbiota_filtrado, aes(x = SampleID, y = input, fill = SampleID)) +
  geom_boxplot() +
  labs(
    title = "Number of Input Sequences by Treatment",
    x = "",
    y = "Input Sequences"
  ) +
  scale_fill_manual(values = treatment_colors) +  # Apply custom colors
  theme_minimal() +
  theme(
    legend.position = "none",  # Hide the legend
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1)  # Rotate x-axis labels
  )

# Para la variable Prop_NonChim
anova_result <- aov(Prop_NonChim ~ SampleID, data = resultados_microbiota_filtrado)
residuos_prop_nonchim <- residuals(anova_result)

# Gráfico Q-Q
qqnorm(residuos_prop_nonchim)
qqline(residuos_prop_nonchim, col = "red")

# Prueba de normalidad de Shapiro-Wilk
shapiro.test(residuos_prop_nonchim)

# Para la variable input
anova_input <- aov(input ~ SampleID, data = resultados_microbiota_filtrado)
residuos_input <- residuals(anova_input)

# Gráfico Q-Q
qqnorm(residuos_input)
qqline(residuos_input, col = "red")

# Prueba de normalidad de Shapiro-Wilk
shapiro.test(residuos_input)


anova_result <- aov(Prop_NonChim ~ SampleID, data = resultados_microbiota_filtrado)
summary(anova_result)

anova_input <- aov(input ~ SampleID, data = resultados_microbiota_filtrado)
summary(anova_input)

# Post-hoc con Tukey
tukey_input <- TukeyHSD(anova_input)
print(tukey_input)

microbiota_abundancia <- microbiota_abundancia %>%
  filter(!Trat %in% c("T01", "T02", "T03"))

# Extraer solo las columnas de abundancia (sin la columna Phylum)
abundancia <- microbiota_abundancia[, -1]

# Calcular el índice de Shannon
shannon_index <- diversity(abundancia, index = "shannon")

# Calcular el índice de Simpson
simpson_index <- diversity(abundancia, index = "simpson")

# Calcular la riqueza de especies (número de taxa presentes)
chao1 <- specnumber(abundancia)

# Agregar los índices calculados al dataframe original
microbiota_abundancia$Shannon <- shannon_index
microbiota_abundancia$Simpson <- simpson_index
microbiota_abundancia$Chao1 <- chao1 # Richness = Chao1

# Ver los resultados
head(microbiota_abundancia)
colnames(microbiota_abundancia)


### pruebas de normalidad previas - supuestos

# Para Shannon
anova_shannon <- aov(Shannon ~ Trat, data = microbiota_abundancia)
residuos_shannon <- residuals(anova_shannon)

# Para Chao1
anova_chao1 <- aov(Chao1 ~ Trat, data = microbiota_abundancia)
residuos_chao1 <- residuals(anova_chao1)

# Para Simpson
anova_simpson <- aov(Simpson ~ Trat, data = microbiota_abundancia)
residuos_simpson <- residuals(anova_simpson)

# 1. Histograma de residuos para cada índice
hist(residuos_shannon, main = "Histograma de Residuos - Shannon", xlab = "Residuos")
hist(residuos_chao1, main = "Histograma de Residuos - Chao1", xlab = "Residuos")
hist(residuos_simpson, main = "Histograma de Residuos - Simpson", xlab = "Residuos")

# 2. Q-Q plot de residuos para cada índice
qqnorm(residuos_shannon)
qqline(residuos_shannon, col = "red")
qqnorm(residuos_chao1)
qqline(residuos_chao1, col = "red")
qqnorm(residuos_simpson)
qqline(residuos_simpson, col = "red")

# 3. Test de normalidad Shapiro-Wilk para cada índice
shapiro.test(residuos_shannon)
shapiro.test(residuos_chao1)
shapiro.test(residuos_simpson)


dist_bray <- vegdist(abundancia, method = "bray")

pca <- prcomp(dist_bray)
plot(pca)

# Extraer las puntuaciones de las muestras en los dos primeros componentes principales
pca_scores <- as.data.frame(pca$x)

# Calcular la varianza explicada de cada componente
pca_variance <- summary(pca)$importance[2, 1:2]  # Para PC1 y PC2
pca_variance_percent <- pca_variance * 100

#Aseguramos que las filas en 'pca_scores' coinciden con el orden de las muestras
# Si 'Trat' está en el dataframe original, lo agregamos a 'pca_scores'
pca_scores$Trat <- microbiota_abundancia$Trat  # Suponiendo que 'Trat' está en 'microbiota_abundancia'

# También puedes asegurarte de que la columna 'Trat' esté en el mismo orden que los datos en 'pca_scores'
# Asegúrate de que los tratamientos estén en el mismo orden para que las correspondencias sean correctas
pca_scores$Trat <- factor(pca_scores$Trat, levels = treatment_order)

# Ahora que tenemos la columna 'Trat', podemos graficar
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Trat)) +
  geom_point(size = 3) +
  scale_color_manual(values = treatment_colors) +
  labs(x = paste("PC1 (", round(pca_variance_percent[1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(pca_variance_percent[2], 2), "%)", sep = ""),
       title = "PCA - First Two Components") +
  theme_light() +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(order = 1))  # Para que la leyenda se ordene

# CON PCoA
# Realizar PCoA utilizando cmdscale
pcoa <- cmdscale(dist_bray, k = 2, eig = TRUE)  # k = 2 para obtener los dos primeros componentes principales

# Extraer las puntuaciones de las muestras en los dos primeros componentes
pcoa_scores <- as.data.frame(pcoa$points)

# Calcular la varianza explicada de cada componente
pcoa_variance <- pcoa$eig[1:2]  # Para los dos primeros componentes
pcoa_variance_percent <- (pcoa_variance / sum(pcoa$eig)) * 100

# Asegurarse de que las filas en 'pcoa_scores' coincidan con el orden de las muestras
pcoa_scores$Trat <- microbiota_abundancia$Trat  # Suponiendo que 'Trat' está en 'microbiota_abundancia'

# Asegurarse de que la columna 'Trat' esté en el mismo orden que los datos de PCoA
pcoa_scores$Trat <- factor(pcoa_scores$Trat, levels = treatment_order)

# Graficar los resultados de PCoA
ggplot(pcoa_scores, aes(x = V1, y = V2, color = Trat)) +
  geom_point(size = 3) +
  scale_color_manual(values = treatment_colors) +
  labs(x = paste("PCoA1 (", round(pcoa_variance_percent[1], 2), "%)", sep = ""),
       y = paste("PCoA2 (", round(pcoa_variance_percent[2], 2), "%)", sep = ""),
       title = "PCoA - Bray-Curtis distance") +
  theme_light() +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(order = 1))  # Para que la leyenda se ordene



pielou_index <- diversity(abundancia) / log(specnumber(abundancia))

# Supongamos que 'Trat' es la variable con los tratamientos
pielou_data <- data.frame(Pielou = pielou_index, Trat = microbiota_abundancia$Trat)

# ANOVA para comparar el índice de Pielou entre tratamientos
anova_pielou <- aov(Pielou ~ Trat, data = pielou_data)
summary(anova_pielou)
# Boxplot del índice de Pielou por tratamientos
boxplot(Pielou ~ Trat, data = pielou_data, 
        main = "Índice de Pielou por Tratamiento", 
        ylab = "Índice de Pielou", 
        col = "lightblue")

# Calcular la desviación estándar del índice de Pielou dentro de cada grupo
pielou_summary <- pielou_data %>%
  group_by(Trat) %>%
  summarise(mean_pielou = mean(Pielou), sd_pielou = sd(Pielou))

# Mostrar resumen
print(pielou_summary)

# Gráfico de dispersión entre Pielou y Shannon
plot(pielou_index, diversity(abundancia), 
     xlab = "Índice de Pielou", ylab = "Diversidad de Shannon", 
     main = "Relación entre Pielou y Shannon")

# Asegurarte de que Chao1 está en tu tabla y se llama correctamente
chao1_values <- microbiota_abundancia$Chao1

# Gráfico de dispersión entre Pielou y Chao1
plot(pielou_index, chao1_values, 
     xlab = "Índice de Pielou", ylab = "Índice de Chao1", 
     main = "Relación entre Pielou y Chao1")

# Gráfico de dispersión entre Pielou y Simpson
plot(pielou_index, diversity(abundancia, index = "simpson"), 
     xlab = "Índice de Pielou", ylab = "Índice de Simpson", 
     main = "Relación entre Pielou y Simpson")


# Supongamos que tienes una columna de factores llamada 'grupo' en tu dataframe
boxplot(Shannon ~ Trat, data = microbiota_abundancia, main = "Índice de Shannon por Grupo")

hist(microbiota_abundancia$Shannon, main = "Distribución del Índice de Shannon", xlab = "Índice de Shannon")

ggplot(microbiota_abundancia, aes(x = Trat, y = Shannon)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Diversidad Alfa: Índice de Shannon")


# Para Shannon
anova_shannon <- aov(Shannon ~ Trat, data = microbiota_abundancia)
summary(anova_shannon)

# Para Chao1
anova_chao1 <- aov(Chao1 ~ Trat, data = microbiota_abundancia)
summary(anova_chao1)

# Para Simpson
anova_simpson <- aov(Simpson ~ Trat, data = microbiota_abundancia)
summary(anova_simpson)

# Ver resultados del análisis post-hoc de Tukey para Shannon
tukey_shannon <- TukeyHSD(anova_shannon)
# Resultados detallados de las comparaciones de Tukey
print(tukey_shannon)


library(rlang)

# Definir las variables de respuesta
response_vars <- c("Shannon", "Simpson", "Chao1")

plot_boxplots <- function(df, response_vars, treatment_order, treatment_colors) {
  lapply(response_vars, function(var) {
    # Realizar la prueba ANOVA para obtener el valor p
    anova_result <- aov(eval_tidy(ensym(var)) ~ Trat, data = df)  # Evalúa correctamente la variable `var`
    p_value <- summary(anova_result)[[1]]$`Pr(>F)`[1]  # Extraer el p-value de la ANOVA
    
    # Crear el gráfico de caja
    p <- ggplot(df, aes(x = Trat, y = eval_tidy(ensym(var)), fill = Trat)) +  # Evalúa correctamente la variable `var`
      geom_boxplot() +
      scale_fill_manual(values = treatment_colors) +
      scale_x_discrete(limits = treatment_order) +
      labs(x = NULL, y = var) +
      theme_light() +
      theme(axis.text.x = element_blank())  # Quitar los rótulos del eje X
    
    # Agregar el valor p en la parte superior del gráfico
    p + annotate("text", x = 3, y = max(df[[var]], na.rm = TRUE) + 0.5, 
                 label = paste("ANOVA, p =", round(p_value, 4)), 
                 size = 4, color = "black")
  })
}

# Llamar la función
plots <- plot_boxplots(microbiota_abundancia, response_vars, treatment_order, treatment_colors)

# Instalar y cargar el paquete cowplot si no lo tienes instalado
# install.packages("cowplot")
library(cowplot)

# Crear la lista de gráficos sin los rótulos del eje X
plots_no_x_labels <- lapply(plots, function(p) 
  p + 
    theme(legend.position = "none",  # Sin leyenda en cada gráfico
          axis.text.x = element_blank())  # Quitar los rótulos del eje X
)

# Asegúrate de que 'Trat' sea un factor con el orden correcto
microbiota_abundancia$Trat <- factor(microbiota_abundancia$Trat, levels = treatment_order)

# Crear la leyenda sin título
legend <- get_legend(
  ggplot(microbiota_abundancia, aes(x = Trat, y = Shannon, fill = Trat)) +  # Un gráfico de ejemplo
    geom_boxplot() + 
    scale_fill_manual(values = treatment_colors) +  # Definir los colores manualmente
    theme_minimal() +
    theme(legend.position = "right") +  # Posición de la leyenda
    guides(fill = guide_legend(title = NULL, override.aes = list(shape = 15, size = 4)))  # Sin título en la leyenda
)

# Crear los gráficos sin leyenda y sin los rótulos del eje X
plots_no_x_labels <- lapply(plots, function(p) 
  p + 
    theme(legend.position = "none",  # Sin leyenda en cada gráfico
          axis.text.x = element_blank())  # Quitar los rótulos del eje X
)

# Organizar los gráficos con la leyenda compartida
combined_plot <- plot_grid(
  plot_grid(plotlist = plots_no_x_labels, ncol = 3),  # Gráficos sin los rótulos del eje X
  legend,  # Leyenda al lado derecho
  ncol = 2,  # Dos columnas, gráficos a la izquierda y leyenda a la derecha
  rel_widths = c(0.85, 0.15)  # Ajustar el tamaño de los gráficos y la leyenda
)

# Mostrar el gráfico combinado con la leyenda sin título
print(combined_plot)

# 1. Eliminar las columnas Chao1, Shannon y Simpson
microbiota_abundancia_phyla <- microbiota_abundancia[, !colnames(microbiota_abundancia) %in% c("Chao1", "Shannon", "Simpson")]

# 2. Transformar los datos a formato largo
microbiota_long <- microbiota_abundancia_phyla %>%
  pivot_longer(cols = -Trat, names_to = "Phylum", values_to = "Abundance")

# 3. Calcular la abundancia promedio por Trat y Phylum
top_phyla <- microbiota_long %>%
  group_by(Trat, Phylum) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance)) %>%
  top_n(5, mean_abundance)  # Seleccionar los 5 phyla con mayor abundancia promedio


# Ordenar 'Trat' según 'treatment_order'
top_phyla$Trat <- factor(top_phyla$Trat, levels = treatment_order)

# Crear el gráfico de barras con el orden en el eje X especificado
ggplot(top_phyla, aes(x = Trat, y = mean_abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "dodge") +  # Barras para cada phylum
  theme_light() +
  labs(x = "", y = "Relative abundance (%)", fill = "Phylum") +  # Etiquetas
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotar los rótulos del eje X para mejor visualización
  annotate("text", x = 1, y = max(top_phyla$mean_abundance) * 1.1,  # Colocar la "A" fuera de los límites de las barras
           label = "A", size = 4, fontface = "bold", color = "black")  # Configuración de la letra "A"


# Normalizar las abundancias al 100% por tratamiento
top_phyla_normalized <- top_phyla %>%
  group_by(Trat) %>%
  mutate(normalized_abundance = mean_abundance / sum(mean_abundance) * 100)  # Convertir a porcentaje

# Crear el gráfico con barras más angostas y normalizadas al 100%
ggplot(top_phyla_normalized, aes(x = Trat, y = normalized_abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +  # Barras más angostas con width ajustado
  theme_light() +
  labs(x = "", y = "Relative abundance (%)", fill = "Phylum") +  # Etiquetas
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotar los rótulos del eje X para mejor visualización
  annotate("text", x = 0.7, y = 105,  # Colocar la "A" fuera de los límites de las barras
           label = "A", size = 4, fontface = "bold", color = "black")  # Configuración de la letra "A"

# Paso 1: Filtrar los tratamientos T01, T02 y T03
microbiota_gen <- microbiota_gen %>%
  filter(!Trat %in% c("T01", "T02", "T03"))

# Paso 2: Convertir los datos de formato ancho a largo
microbiota_gen_long <- microbiota_gen %>%
  gather(key = "Genus", value = "Abundance", -Trat)

# Paso 3: Filtrar géneros con más del 1% de abundancia relativa por réplica
microbiota_gen_long <- microbiota_gen_long %>%
  filter(Abundance >= 0.15)

replica_counts <- microbiota_gen_long %>%
  group_by(Trat, Genus) %>%
  summarise(Replica_Count = n(), .groups = "drop")

# Paso 5: Filtrar géneros presentes en al menos 2 réplicas por tratamiento
valid_genera <- replica_counts %>%
  filter(Replica_Count >= 2) %>%
  dplyr::select(Trat, Genus)

# Paso 6: Filtrar los datos originales con los géneros válidos
filtered_microbiota <- microbiota_gen_long %>%
  semi_join(valid_genera, by = c("Trat", "Genus"))


# Paso 6: Filtrar el conjunto de datos principal para incluir solo géneros válidos
filtered_microbiota <- microbiota_gen_long %>%
  semi_join(valid_genera, by = c("Trat", "Genus"))

# Paso 7: Calcular la suma total de abundancia por tratamiento antes de filtrar
total_abundance_before_filtering <- microbiota_gen_long %>%
  group_by(Trat) %>%
  summarise(Total_Abundance_Before_Filtering = sum(Abundance, na.rm = TRUE))

# Paso 8: Unir la suma total calculada con el dataset filtrado
filtered_microbiota <- filtered_microbiota %>%
  left_join(total_abundance_before_filtering, by = "Trat")

# Paso 9: Normalizar las abundancias relativas usando el total previo al filtrado
filtered_microbiota <- filtered_microbiota %>%
  mutate(Abundance_Percentage = (Abundance / Total_Abundance_Before_Filtering) * 100)

# Paso 10: Asegurarte de que los valores sumen 100% basados en el total previo al filtrado
filtered_microbiota %>%
  group_by(Trat) %>%
  summarise(Total_Abundance = sum(Abundance_Percentage))

# Asegurar el orden de los tratamientos
filtered_microbiota$Trat <- factor(filtered_microbiota$Trat, levels = treatment_order)

# Crear el gráfico de barras apiladas
ggplot(filtered_microbiota, aes(x = Trat, y = Abundance_Percentage, fill = Genus)) +
  geom_bar(stat = "identity", width = 0.6) +  # Ajusta el ancho de las barras
  theme_light() +  # Estilo minimalista
  labs(x = "", y = "Relative abundance (%)", fill = "Genus") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotar etiquetas del eje X
  annotate("text", x = 0.7, y = max(filtered_microbiota$Abundance_Percentage) * 2.95, 
           label = "B", size = 4, fontface = "bold", color = "black")  # Letra "B" en negrita

library(RColorBrewer)
n <- length(unique(filtered_microbiota$Genus))  # Número de géneros únicos
pal <- unlist(brewer.pal.info[1:5, "maxcolors"])  # Extrae más de una paleta y las une
colors <- colorRampPalette(brewer.pal(9, "Set1"))(n)  # Crea una paleta personalizada

ggplot(filtered_microbiota, aes(x = Trat, y = Abundance_Percentage, fill = Genus)) +
  geom_bar(stat = "identity", width = 0.6) + 
  scale_fill_manual(values = colors) +  # Usa la paleta personalizada
  theme_light() + 
  labs(x = "", y = "Relative abundance (%)", fill = "Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Usar viridis para colores diferenciados
  annotate("text", x = 0.7, y = max(filtered_microbiota$Abundance_Percentage) * 2.95, 
           label = "B", size = 4, fontface = "bold", color = "black")  # Letra "B" en negrita

