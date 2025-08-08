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
library(patchwork)

# Definir los colores para cada tratamiento
treatment_colors <- c("FO" = "blue", 
                      "CO" = "#98FB98", 
                      "SH" = "#FFA07A", 
                      "CA" = "darkgreen", 
                      "SHCA" = "#FFCC99",
                      "COCA" = "#32CD32")

# Asegúrate de que los tratamientos están en el orden correcto
treatment_order <- c("FO", "SH", "SHCA", "COCA", "CO", "CA")
treatment_order1 <- c("FO", "SH", "SHCA", "CO")

# Leer los datos iniciales
planilla <- read_excel("planilla.xlsx", sheet = "mg_g")
tabla_genes_R <- read_excel("~/Desktop/UBAINT Doctoral/Ensayo/tabla genes R.xlsx")


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

####################### EXPRESION GENICA #######################

summary(tabla_genes_R)

tabla_promedio_referencia <- tabla_genes_R %>%
  filter(Gen %in% c("bactin", "elf1")) %>%
  group_by(Trat, Rep, Pez) %>%
  summarize(Ct_Ref = mean(Ct_Value, na.rm = TRUE), .groups = "drop")

tabla_delta_ct <- tabla_genes_R %>%
  filter(!Gen %in% c("bactin", "elf1")) %>%
  left_join(tabla_promedio_referencia, by = c("Trat", "Rep", "Pez")) %>%
  mutate(Delta_Ct = Ct_Value - Ct_Ref)

control_promedio <- tabla_delta_ct %>%
  filter(Trat == "FO") %>%
  group_by(Gen) %>%
  summarize(Delta_Ct_Control = mean(Delta_Ct, na.rm = TRUE), .groups = "drop")

tabla_delta_delta_ct <- tabla_delta_ct %>%
  left_join(control_promedio, by = "Gen") %>%
  mutate(Delta_Delta_Ct = Delta_Ct - Delta_Ct_Control)

tabla_fold_change <- tabla_delta_delta_ct %>%
  mutate(Fold_Change = 2^(-Delta_Delta_Ct))

resumen_fold_change <- tabla_fold_change %>%
  group_by(Trat, Gen) %>%
  summarize(
    Mean_Fold_Change = mean(Fold_Change, na.rm = TRUE),
    SD_Fold_Change = sd(Fold_Change, na.rm = TRUE),
    .groups = "drop"
  )

# Filtrar solo los genes de interés (elovl2, elovl5)
tabla_fold_change_filtrada <- tabla_fold_change %>%
  filter(Gen %in% c("elovl2", "elovl5", "cpt1"))


# Realizar la transformación de raíz cuadrada
tabla_fold_change_filtrada$Fold_Change_sqrt <- sqrt(tabla_fold_change_filtrada$Fold_Change)

# Graficar histograma de la variable transformada
ggplot(tabla_fold_change_filtrada, aes(x = Fold_Change_sqrt)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  labs(title = "Histograma de Fold Change con Transformación Raíz Cuadrada",
       x = "Fold Change (Raíz Cuadrada)", y = "Frecuencia") +
  theme_minimal()

# Realizar la prueba de Shapiro-Wilk para normalidad
shapiro_test_sqrt <- shapiro.test(tabla_fold_change_filtrada$Fold_Change_sqrt)
# Mostrar los resultados del test
shapiro_test_sqrt


#### Los datos no siguen distribucion normal ####

# Prueba de normalidad para cada gen
tabla_fold_change_filtrada %>%
  group_by(Gen) %>%
  do({
    shapiro_test <- shapiro.test(.$Fold_Change_sqrt)
    data.frame(Gen = unique(.$Gen), p_value = shapiro_test$p.value)
  })

# 1. Prueba de Levene y normalidad para cada gen
pruebas_diagnostico <- tabla_fold_change_filtrada %>%
  mutate(Trat = as.factor(Trat)) %>%  # Convertir Trat a factor
  group_by(Gen) %>%
  summarise(
    # Prueba de normalidad
    shapiro_p_value = shapiro.test(Fold_Change_sqrt)$p.value,
    
    # Prueba de homogeneidad de varianzas (Levene)
    levene_p_value = leveneTest(Fold_Change_sqrt ~ Trat, data = .)$`Pr(>F)`[1]
  )

# 2. Análisis de residuos: ajustar un modelo y verificar normalidad de residuos
# Supongamos que tienes un modelo lineal:
modelo <- lm(Fold_Change_sqrt ~ Trat * Gen, data = tabla_fold_change_filtrada)

# Extraer residuos
residuos <- residuals(modelo)

# Prueba de normalidad de residuos
shapiro_residuos <- shapiro.test(residuos)

# Imprimir resultados
pruebas_diagnostico
shapiro_residuos

plot(residuos)
qqnorm(residuos)
qqline(residuos)

#### tras transformacion sqrt cpt1 y elovl5 siguen dist normal



# Visualización con boxplot, ordenando los tratamientos
ggplot(tabla_fold_change_filtrada, aes(x = factor(Trat, levels = c("FO", "SH", "SHCA", "COCA", "CO", "CA")), 
                                       y = Fold_Change, fill = Trat)) +
  geom_boxplot() +
  scale_fill_manual(values = treatment_colors) +
  scale_fill_manual(values = treatment_colors, 
                    breaks = c("FO", "SH", "SHCA", "COCA", "CO", "CA")) +  # Ordenar tratamientos en la leyenda
  scale_x_discrete(limits = c("FO", "SH", "SHCA", "COCA", "CO", "CA")) +  # Ordenar tratamientos en las barras
  facet_wrap(~Gen) +  # Para separar los gráficos por gen (elovl2 y elovl5)
  labs(title = "Distribución de Fold Change por Tratamiento",
       x = "", y = "Fold Change") +
  theme_minimal()

#### eliminando outliers

# Eliminar las filas con NA en la columna Fold_Change_sqrt
tabla_fold_change_filtrada_sin_NA <- tabla_fold_change_filtrada %>%
  filter(!is.na(Fold_Change_sqrt))

# Ahora realizar el cálculo de cuartiles y el filtrado de outliers en la tabla sin NA
Q1 <- quantile(tabla_fold_change_filtrada_sin_NA$Fold_Change_sqrt, 0.25)
Q3 <- quantile(tabla_fold_change_filtrada_sin_NA$Fold_Change_sqrt, 0.75)
IQR <- Q3 - Q1

limite_inferior <- Q1 - 1.5 * IQR
limite_superior <- Q3 + 1.5 * IQR

tabla_sin_outliers <- tabla_fold_change_filtrada_sin_NA %>%
  filter(Fold_Change_sqrt >= limite_inferior & Fold_Change_sqrt <= limite_superior)

# Ver los primeros registros de la tabla sin outliers
head(tabla_sin_outliers)

# Visualización con boxplot, ordenando los tratamientos
ggplot(tabla_sin_outliers, aes(x = factor(Trat, levels = c("FO", "SH", "SHCA", "COCA", "CO", "CA")), 
                               y = Fold_Change, fill = Trat)) +
  geom_boxplot() +
  scale_fill_manual(values = treatment_colors) +
  scale_fill_manual(values = treatment_colors, 
                    breaks = c("FO", "SH", "SHCA", "COCA", "CO", "CA")) +  # Ordenar tratamientos en la leyenda
  scale_x_discrete(limits = c("FO", "SH", "SHCA", "COCA", "CO", "CA")) +  # Ordenar tratamientos en las barras
  facet_wrap(~Gen) +  # Para separar los gráficos por gen (elovl2 y elovl5)
  labs(title = "Distribución de Fold Change por Tratamiento",
       x = "", y = "Fold Change") +
  theme_minimal()



# ANOVA para elovl2 y elovl5 sobre la transformación cuadrática
anova_results_sqrt <- tabla_sin_outliers %>%
  group_by(Gen) %>%
  do({
    aov_result <- aov(Fold_Change_sqrt ~ Trat, data = .)
    tidy(aov_result)
  })

# Ver resultados del ANOVA con transformación cuadrática
anova_results_sqrt

# Filtrar los datos para excluir los tratamientos CA y COCA, y quedarse solo con cpt1 en la columna Gen
tabla_cpt1 <- tabla_sin_outliers %>%
  filter(!Trat %in% c("CA", "COCA"), Gen == "cpt1")

# Ver los primeros registros de la tabla filtrada
head(tabla_cpt1)

# prueba de Kruskal-Wallis sobre los datos filtrados
kruskal_test <- kruskal.test(Fold_Change_sqrt ~ Trat, data = tabla_cpt1)

# Ver los resultados
print(kruskal_test)


# Test de Tukey solo para el gen 'elovl5' sobre la transformación cuadrática
tukey_results_elovl5_sqrt <- tabla_sin_outliers %>%
  filter(Gen == "elovl5") %>%  # Filtrar solo para el gen elovl5
  do({
    aov_result <- aov(Fold_Change_sqrt ~ Trat, data = .)
    tukey_result <- TukeyHSD(aov_result)
    tukey_result <- HSD.test(aov_result, "Trat", group = TRUE)
    tukey_result$groups  # Extraer las letras de los grupos
  })

# Mostrar los resultados del test de Tukey para elovl5
tukey_results_elovl5_sqrt

# Test de Tukey solo para el gen 'elovl5' sobre la transformación cuadrática
tukey_results_elovl2_sqrt <- tabla_sin_outliers %>%
  filter(Gen == "elovl2") %>%  # Filtrar solo para el gen elovl5
  do({
    aov_result <- aov(Fold_Change_sqrt ~ Trat, data = .)
    tukey_result <- TukeyHSD(aov_result)
    tukey_result <- HSD.test(aov_result, "Trat", group = TRUE)
    tukey_result$groups  # Extraer las letras de los grupos
  })

# Mostrar los resultados del test de Tukey para elovl5
tukey_results_elovl2_sqrt



#### hasta aca sale bien todo, se comprueba normalidad, se transforma sqrt, se comprueba normalidad nuevamente
#### se procede con elovl5 y elovl2 (tukey) y cpt1 (kruskal) no dio diferencias

#### GRAFICO DE BARRAS PARA CPT1 QUE NO DIO DIFERENCIAS

# Calcular medias y errores estándar por tratamiento
tabla_resumen <- tabla_cpt1 %>%
  group_by(Trat) %>%
  summarise(
    Media = mean(Fold_Change_sqrt, na.rm = TRUE),
    SE = sd(Fold_Change_sqrt, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup()

# Crear el gráfico de barras
ggplot(tabla_resumen, aes(x = factor(Trat, levels = treatment_order), 
                          y = Media, 
                          fill = factor(Trat, levels = treatment_order))) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +  # Barras
  geom_errorbar(aes(ymin = Media - SE, ymax = Media + SE), 
                width = 0.2, color = "black") +  # Barras de error
  scale_fill_manual(values = treatment_colors) +  # Colores personalizados
  labs(title = "cpt1", 
       x = "", 
       y = "Fold Change (sqrt)") +
  theme_minimal() +
  theme(
    legend.position = "right", 
    legend.title = element_blank(), 
    legend.text = element_text(size = 10), 
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5))  # Centrar el título


#### GRAFICOS DE BARRAS CON LETRAS TUKEY PARA GENES ELOVL 2 Y 5

###elovl2###
# Obtener los resultados de Tukey para el gen elovl5 con transformación cuadrática

tukey_letters_elovl2_sqrt <- tabla_sin_outliers %>%
  filter(Gen == "elovl2")

# Calcular los promedios de 'Pez' agrupados por 'Trat' y 'Rep'
promedios_pez <- tabla_sin_outliers %>%
  filter(Gen == "elovl2") %>%
  group_by(Trat, Rep) %>%
  summarise(
    Promedio_Pez = mean(Pez, na.rm = TRUE),
    .groups = "drop"  # Desagrupar después de calcular los promedios
  )

# Calcular medias y errores estándar de Fold_Change_sqrt por Trat (ya lo tienes hecho)
resumen_elovl2_sqrt <- tabla_sin_outliers %>%
  filter(Gen == "elovl2") %>%
  group_by(Trat) %>%
  summarise(
    Media = mean(Fold_Change_sqrt, na.rm = TRUE),
    SE = sd(Fold_Change_sqrt, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup()

# Realizar la prueba de Tukey y añadir las letras de Tukey
tukey_letters_elovl2_sqrt <- tabla_sin_outliers %>%
  filter(Gen == "elovl2") %>%
  do({
    aov_result <- aov(Fold_Change_sqrt ~ Trat, data = .)
    tukey_result <- HSD.test(aov_result, "Trat", group = TRUE)
    tukey_result$groups %>%
      as.data.frame() %>%
      mutate(Trat = rownames(.))  # Mantener los tratamientos como columna
  })

# Unir los resultados de Tukey con las medias y errores estándar
tabla_fold_change_filtrada_with_letters_sqrt2 <- resumen_elovl2_sqrt %>%
  left_join(tukey_letters_elovl2_sqrt, by = "Trat")

# Verifica el dataframe resultante para asegurarte de que Fold_Change_sqrt, Media y SE están presentes
head(tabla_fold_change_filtrada_with_letters_sqrt2)

# Crear el gráfico con ggplot y título centrado, sin título en la leyenda
ggplot(tabla_fold_change_filtrada_with_letters_sqrt2, aes(x = Trat, y = Media, fill = Trat)) +
  geom_bar(stat = "identity", position = "dodge") +  # Barras de media
  geom_errorbar(aes(ymin = Media - SE, ymax = Media + SE), width = 0.2, color = "black") +  # Líneas de error
  scale_fill_manual(values = treatment_colors, 
                    breaks = c("FO", "SH", "SHCA", "COCA", "CO", "CA")) +  # Ordenar tratamientos en la leyenda
  scale_x_discrete(limits = c("FO", "SH", "SHCA", "COCA", "CO", "CA")) +  # Ordenar tratamientos en las barras
  labs(title = "elovl2",
       x = "", y = "Average Fold Change") +
  theme_minimal() +
  geom_text(aes(label = groups), 
            position = position_dodge(width = 0.8), 
            vjust = -0.4, # Mover las letras por encima de las barras
            hjust = -0.4,
            size = 4) +  # Añadir letras de Tukey sobre las barras
  theme(legend.position = "right",  # Mostrar la leyenda a la derecha
        legend.title = element_blank(),  # Eliminar el título de la leyenda
        plot.title = element_text(hjust = 0.5, face = "italic"))  # Título en cursiva



###elovl5###
# Obtener los resultados de Tukey para el gen elovl5 con transformación cuadrática
tukey_letters_elovl5_sqrt <- tabla_sin_outliers %>%
  filter(Gen == "elovl5")

# Calcular los promedios de 'Pez' agrupados por 'Trat' y 'Rep'
promedios_pez1 <- tabla_sin_outliers %>%
  filter(Gen == "elovl5") %>%
  group_by(Trat, Rep) %>%
  summarise(
    Promedio_Pez = mean(Pez, na.rm = TRUE),
    .groups = "drop"  # Desagrupar después de calcular los promedios
  )

# Calcular medias y errores estándar de Fold_Change_sqrt por Trat (ya lo tienes hecho)
resumen_elovl5_sqrt <- tabla_sin_outliers %>%
  filter(Gen == "elovl5") %>%
  group_by(Trat) %>%
  summarise(
    Media = mean(Fold_Change_sqrt, na.rm = TRUE),
    SE = sd(Fold_Change_sqrt, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup()

# Realizar la prueba de Tukey y añadir las letras de Tukey
tukey_letters_elovl5_sqrt <- tabla_sin_outliers %>%
  filter(Gen == "elovl5") %>%
  do({
    aov_result <- aov(Fold_Change_sqrt ~ Trat, data = .)
    tukey_result <- HSD.test(aov_result, "Trat", group = TRUE)
    tukey_result$groups %>%
      as.data.frame() %>%
      mutate(Trat = rownames(.))  # Mantener los tratamientos como columna
  })

# Unir los resultados de Tukey con las medias y errores estándar
tabla_fold_change_filtrada_with_letters_sqrt5 <- resumen_elovl5_sqrt %>%
  left_join(tukey_letters_elovl5_sqrt, by = "Trat")

# Verifica el dataframe resultante para asegurarte de que Fold_Change_sqrt, Media y SE están presentes
head(tabla_fold_change_filtrada_with_letters_sqrt5)

# Crear el gráfico con ggplot y título centrado, sin título en la leyenda
ggplot(tabla_fold_change_filtrada_with_letters_sqrt5, aes(x = Trat, y = Media, fill = Trat)) +
  geom_bar(stat = "identity", position = "dodge") +  # Barras de media
  geom_errorbar(aes(ymin = Media - SE, ymax = Media + SE), width = 0.2, color = "black") +  # Líneas de error
  scale_fill_manual(values = treatment_colors, 
                    breaks = c("FO", "SH", "SHCA", "COCA", "CO", "CA")) +  # Ordenar tratamientos en la leyenda
  scale_x_discrete(limits = c("FO", "SH", "SHCA", "COCA", "CO", "CA")) +  # Ordenar tratamientos en las barras
  labs(title = "elovl5",
       x = "", y = "Average Fold Change") +
  theme_minimal() +
  geom_text(aes(label = groups), 
            position = position_dodge(width = 0.8), 
            vjust = -0.4, # Mover las letras por encima de las barras
            hjust = -0.4,
            size = 4) +  # Añadir letras de Tukey sobre las barras
  theme(legend.position = "right",  # Mostrar la leyenda a la derecha
        legend.title = element_blank(),  # Eliminar el título de la leyenda
        plot.title = element_text(hjust = 0.5, face = "italic"))  # Título en cursiva




#### SALIDA DE GRAFICO COMPARTIENDO LEYENDA


# Crear el gráfico de barras para elovl5 con líneas de error
plot1 <- ggplot(tabla_fold_change_filtrada_with_letters_sqrt2, aes(x = Trat, y = Media, fill = Trat)) +
  geom_bar(stat = "identity", position = "dodge") +  # Barras de media
  geom_errorbar(aes(ymin = Media - SE, ymax = Media + SE), width = 0.2, color = "black") +  # Líneas de error
  scale_fill_manual(values = treatment_colors, 
                    breaks = c("FO", "SH", "SHCA", "COCA", "CO", "CA")) +  # Ordenar tratamientos en la leyenda
  scale_x_discrete(limits = c("FO", "SH", "SHCA", "COCA", "CO", "CA")) +  # Ordenar tratamientos en las barras
  labs(title = "elovl2",
       x = "", y = "Average Fold Change") +
  theme_minimal() +
  geom_text(aes(label = groups), 
            position = position_dodge(width = 0.8), 
            vjust = -0.4, # Mover las letras por encima de las barras
            hjust = -0.4,
            size = 4) +  # Añadir letras de Tukey sobre las barras
  theme(legend.position = "none",  # Mostrar la leyenda a la derecha
        legend.title = element_blank(),  # Eliminar el título de la leyenda
        plot.title = element_text(hjust = 0.5, face = "italic"))  # Título en cursiva

# Crear el gráfico de barras para elovl2 con líneas de error
plot2 <- ggplot(tabla_fold_change_filtrada_with_letters_sqrt5, aes(x = Trat, y = Media, fill = Trat)) +
  geom_bar(stat = "identity", position = "dodge") +  # Barras de media
  geom_errorbar(aes(ymin = Media - SE, ymax = Media + SE), width = 0.2, color = "black") +  # Líneas de error
  scale_fill_manual(values = treatment_colors, 
                    breaks = c("FO", "SH", "SHCA", "COCA", "CO", "CA")) +  # Ordenar tratamientos en la leyenda
  scale_x_discrete(limits = c("FO", "SH", "SHCA", "COCA", "CO", "CA")) +  # Ordenar tratamientos en las barras
  labs(title = "elovl5",
       x = "", y = "") +
  theme_minimal() +
  geom_text(aes(label = groups), 
            position = position_dodge(width = 0.8), 
            vjust = -0.4, # Mover las letras por encima de las barras
            hjust = -0.4,
            size = 4) +  # Añadir letras de Tukey sobre las barras
  theme(legend.position = "right",  # Mostrar la leyenda a la derecha
        legend.title = element_blank(),  # Eliminar el título de la leyenda
        plot.title = element_text(hjust = 0.5, face = "italic"))  # Título en cursiva


# Combinar los gráficos en una sola salida con las líneas de error estándar
plot1 + plot2 + plot_layout(guides = 'collect')



####################### ACIDOS GRASOS #######################

# Primero, filtrar los datos y organizar el dataframe
df_higado <- df2 %>% filter(grepl("HIGADO", MUESTRA))
df_higado_gral <- df_higado %>%
  filter(DIETA != "INI") %>%
  dplyr::select(-MUESTRA, -CONDICION, -REPLICA)


# Crear una lista para guardar resultados de ANOVA
anova_results_gral <- list()

# Realizar ANOVA para cada variable en df_higado_percent y almacenar resultados
for (variable in fatty_acid_variables) {
  # Crear una fórmula usando el nombre original de la variable, rodeando con backticks
  formula <- as.formula(paste("`", variable, "` ~ DIETA", sep = ""))
  
  # Realizar ANOVA y almacenar resultados
  anova_results_gral[[variable]] <- summary(aov(formula, data = df_higado_gral))
}

# Ver resultados de ANOVA
anova_results_gral

# Inicializar una lista para almacenar los resultados de las comparaciones
tukey_results_gral<- list()

for (variable in fatty_acid_variables) {
  # Crear una fórmula usando substitute para evitar problemas con caracteres especiales
  formula <- as.formula(substitute(y ~ DIETA, list(y = as.name(variable))))
  
  # Realizar ANOVA
  aov_result_gral <- aov(formula, data = df_higado_gral)
  
  # Realizar el test de Tukey
  tukey_results_gral[[variable]] <- HSD.test(aov_result_gral, "DIETA", group = TRUE)
}

# Mostrar resultados
tukey_results_gral






# Assuming df_higado is already filtered for HIGADO samples and does not include SHCA and COCA
# Convert DIETA and other relevant columns to factors
df_higado$DIETA <- as.factor(df_higado$DIETA)

# Lista de variables de interés
fatty_acid_variables <- c(
  "C16:0 Acido Palmítico",         # Saturado
  "C18:0 Acido Esteárico",         # Saturado
  "C14:0 Acido Mirístico",         # Saturado
  "C17:0 Acido Heptadecanoico",
  "C21:0 Acido Heneicosanoico",# Saturado
  "C22:0 Acido Behénico",          # Saturado
  "C23:0 Acido Tricosanoico",       # Saturado
  "C16:1 Acido Palmitoléico",      # Monoinsaturado
  "C17:1 Acido Cis-10-Heptadecenoico",
  "C18:1n9t Acido Elaidico",# Monoinsaturado
  "C18:1n9c Acido Oléico",         # Monoinsaturado
  "C20:1n9 Acido Cis-11-Eicosenoico",  # Monoinsaturado
  "C18:2n6c Acido Linoléico",      # Poliinsaturado (PUFA)
  "C18:3n3 Acido Linolénico",
  "C18:3n6 Acido g-linolénico",# Poliinsaturado (PUFA)
  "C20:2 Acido Cis-11,14-Eicosadienoico",  # Poliinsaturado (PUFA)
  "C22:2 Acido Cis-13,16-Docosadienoico",  # Poliinsaturado (PUFA)
  "C20:4n6 Acido Araquidónico",   # Poliinsaturado (PUFA)
  "C20:5n3 Acido Cis-5,8,11,14,17-Eicosapentaenoico",  # Poliinsaturado (PUFA)
  "C22:6n3 Acido Cis-4,7,10,13,16,19-Docosahexaenoico",  # Poliinsaturado (PUFA)
  "Total de SAFAs",               # Total saturado
  "Total de MUFAs",               # Total monoinsaturado
  "Total de PUFAs",
  "total_n3",
  "total_n6",
  "n6_n3",
  "% Lipidos totales"             # Total de lípidos
)



# Calcular la media de cada variable en el tratamiento control "FO"
ini_means <- df_higado %>%
  filter(DIETA == "FO") %>%  # Filtrar solo las filas donde DIETA es "FO"
  summarise(across(all_of(fatty_acid_variables), \(x) mean(x, na.rm = TRUE)))  # Calcular la media de cada variable


# Crear un nuevo dataframe con las proporciones respecto a la media del tratamiento control "FO"
df_higado_percent <- df_higado %>%
  filter(DIETA != "INI", DIETA != "FO") %>%  # Excluir las condiciones "INI" y "FO"
  dplyr::select(DIETA, all_of(fatty_acid_variables)) %>%  # Seleccionar DIETA y las variables de interés
  mutate(across(all_of(fatty_acid_variables), 
                ~ ((. - ini_means[[cur_column()]]) / ini_means[[cur_column()]]) * 100))  # Calcular el % de diferencia respecto a "FO"


# Calcular la media y el desvío estándar para cada tratamiento y variable
df_higado_stats <- df_higado_percent %>%
  group_by(DIETA) %>%
  summarise(across(all_of(fatty_acid_variables), 
                   list(mean = ~ mean(.x, na.rm = TRUE), 
                        sd = ~ sd(.x, na.rm = TRUE)), 
                   .names = "{.col}_{.fn}"))

# Primero, transformar el dataframe a formato largo
df_higado_stats_long <- df_higado_stats %>%
  pivot_longer(cols = -DIETA,  # Mantener la columna DIETA
               names_to = "Variable",  # Crear columna 'Variable' con los nombres de las variables
               values_to = "Value")  # Columna con los valores de las variables

# Separar medias y desviaciones estándar
df_higado_stats_mean <- df_higado_stats_long %>%
  filter(str_detect(Variable, "_mean$")) %>%  # Filtrar para obtener solo las medias
  rename(Mean = Value) %>%  # Renombrar la columna Value a Mean
  mutate(Variable = str_remove(Variable, "_mean$"))  # Eliminar el sufijo "_mean" de las variables

df_higado_stats_sd <- df_higado_stats_long %>%
  filter(str_detect(Variable, "_sd$")) %>%  # Filtrar para obtener solo las desviaciones estándar
  rename(SD = Value) %>%  # Renombrar la columna Value a SD
  mutate(Variable = str_remove(Variable, "_sd$"))  # Eliminar el sufijo "_sd" de las variables

# Unir las medias y desviaciones estándar en un solo dataframe
df_higado_stats_combined <- df_higado_stats_mean %>%
  left_join(df_higado_stats_sd, by = c("DIETA", "Variable"))

head(df_higado_stats_combined)

# Crear una lista para guardar resultados de ANOVA
anova_results_percent <- list()

# Realizar ANOVA para cada variable en df_higado_percent y almacenar resultados
for (variable in fatty_acid_variables) {
  # Crear una fórmula usando el nombre original de la variable, rodeando con backticks
  formula <- as.formula(paste("`", variable, "` ~ DIETA", sep = ""))
  
  # Realizar ANOVA y almacenar resultados
  anova_results_percent[[variable]] <- summary(aov(formula, data = df_higado_percent))
}

# Ver resultados de ANOVA
anova_results_percent


# Inicializar una lista para almacenar los resultados de las comparaciones
tukey_results_percent <- list()

for (variable in fatty_acid_variables) {
  # Crear una fórmula usando substitute para evitar problemas con caracteres especiales
  formula <- as.formula(substitute(y ~ DIETA, list(y = as.name(variable))))
  
  # Realizar ANOVA
  aov_result_percent <- aov(formula, data = df_higado_percent)
  
  # Realizar el test de Tukey
  tukey_results_percent[[variable]] <- HSD.test(aov_result_percent, "DIETA", group = TRUE)
}

# Mostrar resultados
tukey_results_percent

# Inicializar un dataframe vacío para almacenar los resultados
anova_combined <- data.frame()

# Iterar sobre cada variable en los resultados de ANOVA
for (variable in names(anova_results_percent)) {
  # Extraer el resumen del ANOVA (el primer y único elemento de la lista)
  anova_summary <- anova_results_percent[[variable]][[1]]
  
  # Obtener las filas relevantes (df, sum sq, f value, p value)
  anova_data <- data.frame(
    Variable = variable,
    Df = anova_summary$Df[1],         # Grados de libertad para el factor
    Sum_Sq = anova_summary$`Sum Sq`[1], # Suma de cuadrados para el factor
    F_value = anova_summary$`F value`[1], # Valor F
    Pr_F = anova_summary$`Pr(>F)`[1],  # Valor p para el factor
    stringsAsFactors = FALSE
  )
  
  # Agregar los resultados al dataframe combinado
  anova_combined <- rbind(anova_combined, anova_data)
}

# Ver los resultados
head(anova_combined)

# Crear una lista vacía para almacenar los resultados
comparisons <- list()

# Recorrer cada variable y extraer las comparaciones y grupos
for (var in fatty_acid_variables) {
  comparisons[[var]] <- tukey_results_percent[[var]]$groups
}

# Mostrar las comparaciones para cada variable
comparisons

# Lista de las variables que tenemos
variable_list <- names(comparisons)

# Crear un dataframe vacío para almacenar los resultados
result_data <- data.frame()

# Iterar sobre cada variable en comparisons
for (var in variable_list) {
  # Extraer los datos de la comparación para cada variable
  temp <- comparisons[[var]]
  
  # Crear un dataframe temporal para esta variable
  temp_df <- data.frame(
    Variable = rep(var, nrow(temp)),  # Replicar el nombre de la variable para todas las filas
    DIETA = rownames(temp),           # Usar los nombres de las filas como los tratamientos (DIETA)
    Mean = temp[, 1],                # Tomar las medias de la primera columna
    Groups = temp[, 2]               # Tomar las letras de Tukey de la segunda columna
  )
  
  # Añadir al dataframe final
  result_data <- rbind(result_data, temp_df)
}

# Verificar el dataframe resultante
head(result_data)

# Definir el orden de los tratamientos
result_data$DIETA <- factor(result_data$DIETA, levels = c("FO", "SH", "SHCA", "COCA", "CO", "CA"))

# Crear un vector con los nuevos nombres de las variables
variable_labels <- c(
  # Ácidos Grasos Saturados (SAFAs)
  "C12:0 Acido Láurico" = "Lauric",
  "C14:0 Acido Mirístico" = "Myristic",
  "C15:0 Acido Pentadecanoico" = "Pentadecanoic",
  "C16:0 Acido Palmítico" = "Palmitic",
  "C17:0 Acido Heptadecanoico" = "Heptadecanoic",
  "C18:0 Acido Esteárico" = "Stearic",
  "C21:0 Acido Heneicosanoico" = "Heneicosanoic",
  "C22:0 Acido Behénico" = "Behenic",
  "C23:0 Acido Tricosanoico" = "Tricosanoic",
  "C24:0 Acido Lignocérico" = "Lignoceric",
  
  # Ácidos Grasos Monoinsaturados (MUFAs)
  "C16:1 Acido Palmitoléico" = "Palmitoleic",
  "C17:1 Acido Cis-10-Heptadecenoico" = "Heptadecenoic",
  "C18:1n9t Acido Elaidico" = "Elaidic",
  "C18:1n9c Acido Oléico" = "Oleic",
  "C20:1n9 Acido Cis-11-Eicosenoico" = "Eicosenoic",
  
  # Ácidos Grasos Poliinsaturados (PUFAs)
  "C18:2n6c Acido Linoléico" = "18:2 n6 Linoleic",
  "C18:3n6 Acido g-linolénico" = "Gamma Linoleic",
  "C18:3n6 Acido g-linolénico" = "18:3 n6 Linolenic",
  "C18:3n3 Acido Linolénico" = "18:3 n3 Linolenic",
  "C20:2 Acido Cis-11,14-Eicosadienoico" = "Eicosadienoic",
  "C20:4n6 Acido Araquidónico" = "ARA",
  "C20:5n3 Acido Cis-5,8,11,14,17-Eicosapentaenoico" = "EPA",
  "C22:2 Acido Cis-13,16-Docosadienoico" = "Docosadienoic",
  "C22:6n3 Acido Cis-4,7,10,13,16,19-Docosahexaenoico" = "DHA",
  
  # Totales y relaciones
  "Total de SAFAs" = "Total SFAs",
  "Total de MUFAs" = "Total MUFAs",
  "Total de PUFAs" = "Total PUFAs",
  "total_n3" = "Total n3",
  "total_n6"= "Total n6",
  "n6_n3" = "n6:n3 ratio"
)

variable_labels_esp <- c(
  # Ácidos Grasos Saturados (SAFAs)
  "C12:0 Acido Láurico" = "Láurico",
  "C14:0 Acido Mirístico" = "Mirístico",
  "C15:0 Acido Pentadecanoico" = "Pentadecanoico",
  "C16:0 Acido Palmítico" = "Palmítico",
  "C17:0 Acido Heptadecanoico" = "Heptadecanoico",
  "C18:0 Acido Esteárico" = "Esteárico",
  "C22:0 Acido Behénico" = "Behénico",
  "C23:0 Acido Tricosanoico" = "Tricosanoico",
  "C24:0 Acido Lignocérico" = "Lignocérico",
  
  # Ácidos Grasos Monoinsaturados (MUFAs)
  "C16:1 Acido Palmitoléico" = "Palmitoléico",
  "C17:1 Acido Cis-10-Heptadecenoico" = "Heptadecenoico",
  "C18:1n9t Acido Elaidico" = "Elaídico",
  "C18:1n9c Acido Oléico" = "Oléico",
  "C20:1n9 Acido Cis-11-Eicosenoico" = "Eicosenoico",
  
  # Ácidos Grasos Poliinsaturados (PUFAs)
  "C18:2n6c Acido Linoléico" = "Linoleico",
  "C18:3n6 Acido g-linolénico" = "gamma-Linolénico",
  "C18:3n3 Acido Linolénico" = "Linolénico",
  "C20:2 Acido Cis-11,14-Eicosadienoico" = "Eicosadienoico",
  "C20:4n6 Acido Araquidónico" = "Araquidónico",
  "C20:5n3 Acido Cis-5,8,11,14,17-Eicosapentaenoico" = "Eicosapentaenoico",
  "C22:2 Acido Cis-13,16-Docosadienoico" = "Docosadienoico",
  "C22:6n3 Acido Cis-4,7,10,13,16,19-Docosahexaenoico" = "Docosahexaenoico",
  
  # Totales y relaciones
  "Total de SAFAs" = "Total SFAs",
  "Total de MUFAs" = "Total MUFAs",
  "Total de PUFAs" = "Total PUFAs",
  "n6" = "n6",
  "n3" = "n3",
  "n6_n3" = "Relación n6:n3"
)


# Filtrar el dataframe para excluir la variable "% Lipidos totales"
result_data_filtered <- result_data %>%
  filter(Variable != "% Lipidos totales")

result_data_filtered <- result_data_filtered %>%
  mutate(Variable = factor(Variable, levels = fatty_acid_variables))


ggplot(result_data_filtered, aes(x = DIETA, y = Mean, fill = DIETA)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # Barras con medias
  geom_text(aes(label = Groups), 
            position = position_dodge(width = 0.8),  # Asegura que las letras se alineen con las barras
            vjust = 1.5,  # Ajusta la posición vertical de las letras Tukey por encima de las barras (más abajo)
            size = 4,  # Ajusta el tamaño de las letras
            fontface = "bold") +  # Hace las letras Tukey en negrita
  facet_wrap(~ Variable, scales = "free_y", labeller = labeller(Variable = variable_labels)) +  # Crear un gráfico por cada variable con nombres personalizados
  scale_fill_manual(values = treatment_colors) +  # Colores personalizados para los tratamientos
  labs(y = "Mean of liver % variance respect of control", x = "Tratamiento") +  # Etiquetas de los ejes
  theme_minimal() +  # Tema minimalista
  theme(axis.text.x = element_blank(),  # Eliminar nombres de los tratamientos en el eje X
        strip.text.x = element_text(size = 12, margin = margin(b = 10)),  # Ajustar el tamaño de los títulos de facetas y margen inferior
        plot.margin = margin(0.5, 0.5, 1, 0.5, "cm"),  # Aumentar margen inferior para evitar corte de etiquetas
        axis.title.x = element_blank(),  # Eliminar el título del eje X
        panel.spacing = unit(1, "lines"),  # Ajustar el espacio entre facetas
        legend.title = element_blank()) +  # Eliminar el título de la leyenda
  coord_cartesian(clip = "off")  # Permite que las etiquetas se muestren fuera del área del gráfico si es necesario


# Crear los gráficos de barras sin la variable "% Lipidos totales"
ggplot(result_data_filtered, aes(x = DIETA, y = Mean, fill = DIETA)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # Barras con medias
  geom_text(aes(label = Groups), vjust = 0.3, size = 4) +  # Añadir las letras de Tukey sobre las barras
  facet_wrap(~ Variable, scales = "free_y", labeller = labeller(Variable = variable_labels_esp)) +  # Crear un gráfico por cada variable con nombres personalizados
  scale_fill_manual(values = treatment_colors) +  # Colores personalizados para los tratamientos
  labs(y = "Media de la diferencia respecto del control (%)", x = "Tratamiento") +  # Etiquetas de los ejes
  theme_minimal() +  # Tema minimalista
  theme(axis.text.x = element_blank(),  # Eliminar nombres de los tratamientos en el eje X
        strip.text.x = element_text(size = 10),  # Ajustar el tamaño de los títulos de facetas
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),  # Ajustar márgenes
        axis.title.x = element_blank(),  # Eliminar el título del eje X
        panel.spacing = unit(0.5, "lines"),
        legend.title = element_blank())



#### CORRELACIONES AC GRASO HIGADO - EXPRESION


df_exp <- tabla_sin_outliers %>%
  filter(Gen != "cpt1")

# 1. Convertir a formato ancho
df_exp_wide <- df_exp %>%
  dplyr::select(Trat, Rep, Pez, Gen, Fold_Change) %>%  # Si preferís Fold_Change_sqrt, cámbialo aquí
  pivot_wider(
    names_from = Gen,
    values_from = Fold_Change
  )

# 2. Promediar por tratamiento y réplica
df_exp_summary <- df_exp_wide %>%
  group_by(Trat, Rep) %>%
  summarize(across(c(elovl2, elovl5), mean, na.rm = TRUE), .groups = "drop") %>%
  rename(DIETA = Trat, REPLICA = Rep) %>%
  rename_with(~ paste0(., ""), .cols = c(elovl2, elovl5))

df_dieta <- df2 %>% filter(grepl("DIETA", MUESTRA))

# Alternativa con subset
df_dieta_exp <- df_dieta %>%
  subset(select = -c(MUESTRA, CONDICION, REPLICA)) %>%
  filter(DIETA != "INI")

# Alternativa con subset
df_higado_exp <- df_higado %>%
  subset(select = -c(MUESTRA, CONDICION)) %>%
  filter(DIETA != "INI")


# 2. Renombrar columnas de df_higado
df_higado_renamed <- df_higado %>%
  rename_with(~ paste0(., "_higado"), .cols = -c(DIETA, REPLICA))

# 3. Expandir df_dieta (sin réplica) para empatar con 3 réplicas
df_dieta_expanded <- df_dieta %>%
  rename_with(~ paste0(., "_dieta"), .cols = -DIETA) %>%
  slice(rep(1:n(), each = 3)) %>%
  mutate(REPLICA = rep(1:3, times = nrow(df_dieta)))  # Agrega réplica ficticia

# 4. Unir todo en combined_data
combined_data <- df_higado_renamed %>%
  inner_join(df_exp_summary, by = c("DIETA", "REPLICA")) %>%
  inner_join(df_dieta_expanded, by = c("DIETA", "REPLICA"))

# Verificar estructura
glimpse(combined_data)

safe_cor <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (
    all(is.na(x)) || all(is.na(y)) ||
    is.na(sd(x, na.rm = TRUE)) || is.na(sd(y, na.rm = TRUE)) ||
    sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0
  ) {
    return(NA)
  } else {
    return(cor(x, y, use = "complete.obs"))
  }
}

# Asegurar que solo columnas numéricas estén incluidas
acidos_higado <- grep("_higado$", names(combined_data), value = TRUE)
acidos_higado <- acidos_higado[sapply(combined_data[acidos_higado], is.numeric)]


# Extraer expresión génica
expresion_genes <- c("elovl2", "elovl5")

cor_matrix <- map_dfr(
  expresion_genes,
  function(gene) {
    map_dbl(acidos_higado, ~ safe_cor(combined_data[[gene]], combined_data[[.x]])) %>%
      setNames(acidos_higado) %>%
      as.data.frame() %>%
      mutate(Gene = gene)
  }
)

# Reformatear
cor_long <- cor_matrix %>%
  pivot_longer(-Gene, names_to = "Variable", values_to = "Correlation") %>%
  mutate(Variable = gsub("_higado$", "", Variable))

# Crear función para obtener correlación y p-valor
get_cor_p <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (
    all(is.na(x)) || all(is.na(y)) ||
    is.na(sd(x, na.rm = TRUE)) || is.na(sd(y, na.rm = TRUE)) ||
    sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0
  ) {
    return(tibble(cor = NA, p = NA))
  }
  test <- cor.test(x, y, method = "pearson")
  tibble(cor = test$estimate, p = test$p.value)
}
# Aplicar a cada combinación
cor_table <- expand.grid(Gene = expresion_genes, FA = acidos_higado, stringsAsFactors = FALSE) %>%
  mutate(
    result = map2(Gene, FA, ~ get_cor_p(combined_data[[.x]], combined_data[[.y]]))
  ) %>%
  unnest(cols = result) %>%
  mutate(FA_clean = gsub("_higado$", "", FA))

# Ajustar p-valores (a nivel global para todas comparaciones)
cor_table <- cor_table %>%
  mutate(
    p_adj = p.adjust(p, method = "fdr"),  # Cambia a "bonferroni" si preferís
    signif = case_when(
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

cor_sig <- cor_table %>%
  filter(!is.na(p) & p < 0.05)


# Ajustar p-valores (a nivel global para todas comparaciones)
cor_table <- cor_table %>%
  mutate(
    p_adj = p.adjust(p, method = "bonferroni"),  # Cambia a "bonferroni" si preferís
    signif = case_when(
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

ggplot(cor_table, aes(x = Gene, y = FA_clean, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  geom_text(aes(label = ifelse(is.na(cor), "", paste0(sprintf("%.2f", cor), signif))), size = 3) +
  labs(title = "Correlation between Liver Fatty Acids and elovl2/elovl5 Expression",
       x = "Gene", y = "Fatty Acid") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5))

label_mapping <- c(
  "C12:0 Acido Láurico_higado" = "Lauric Acid",
  "C14:0 Acido Mirístico_higado" = "Myristic Acid",
  "C15:0 Acido Pentadecanoico_higado" = "Pentadecanoic Acid",
  "C16:0 Acido Palmítico_higado" = "Palmitic Acid",
  "C17:0 Acido Heptadecanoico_higado" = "Heptadecanoic Acid",
  "C18:0 Acido Esteárico_higado" = "Stearic Acid",
  "C21:0 Acido Heneicosanoico_higado" = "Heneicosanoic Acid",
  "C22:0 Acido Behénico_higado" = "Behenic Acid",
  "C23:0 Acido Tricosanoico_higado" = "Tricosanoic Acid",
  "C24:0 Acido Lignocérico_higado" = "Lignoceric Acid",
  
  # Ácidos Grasos Monoinsaturados (MUFAs)
  "C16:1 Acido Palmitoléico_higado" = "Palmitoleic Acid",
  "C17:1 Acido Cis-10-Heptadecenoico_higado" = "Heptadecenoic Acid",
  "C18:1n9t Acido Elaidico_higado" = "Elaidic Acid",
  "C18:1n9c Acido Oléico_higado" = "Oleic Acid",
  "C20:1n9 Acido Cis-11-Eicosenoico_higado" = "Eicosenoic Acid",
  
  # Ácidos Grasos Poliinsaturados (PUFAs)
  "C18:2n6c Acido Linoléico_higado" = "18:2 n6 Linoleic Acid",
  "C18:3n6 Acido g-linolénico_higado" = "Gamma Linoleic Acid",
  "C18:3n6 Acido g-linolénico_higado" = "18:3 n6 Linolenic Acid",
  "C18:3n3 Acido Linolénico_higado" = "18:3 n3 Linolenic Acid",
  "C20:2 Acido Cis-11,14-Eicosadienoico_higado" = "Eicosadienoic Acid",
  "C20:4n6 Acido Araquidónico_higado" = "Arachidonic Acid",
  "C20:5n3 Acido Cis-5,8,11,14,17-Eicosapentaenoico_higado" = "Eicosapentaenoic Acid",
  "C22:2 Acido Cis-13,16-Docosadienoico_higado" = "Docosadienoic Acid",
  "C22:6n3 Acido Cis-4,7,10,13,16,19-Docosahexaenoico_higado" = "Docosahexaenoic Acid",
  
  # Totales y relaciones
  "Total de SAFAs_higado" = "Total SFAs",
  "Total de MUFAs_higado" = "Total MUFAs",
  "Total de PUFAs_higado" = "Total PUFAs",
  "total_n3_higado" = "Total n3",
  "total_n6_higado" = "Total n6",
  "n6_n3_higado" = "n6:n3 Ratio",
  "% Lipidos totales_higado" = "Total lipids"
)

# Definir el orden de los ácidos grasos de más saturado a menos saturado, seguido por totales y ratio
custom_order <- c(
  # Saturados (SAFAs)
  "Lauric Acid", "Myristic Acid", "Pentadecanoic Acid", "Palmitic Acid", "Heptadecanoic Acid", 
  "Stearic Acid", "Heneicosanoic Acid", "Behenic Acid", "Tricosanoic Acid", "Lignoceric Acid",
  
  # Monoinsaturados (MUFAs)
  "Palmitoleic Acid", "Heptadecenoic Acid", "Elaidic Acid", "Oleic Acid", "Eicosenoic Acid",
  
  # Poliinsaturados (PUFAs)
  "18:2 n6 Linoleic Acid", "Gamma Linoleic Acid", "18:3 n6 Linolenic Acid", "18:3 n3 Linolenic Acid",
  "Eicosadienoic Acid", "Arachidonic Acid", "Eicosapentaenoic Acid", "Docosadienoic Acid", "Docosahexaenoic Acid",
  
  # Totales y ratio
  "Total SFAs", "Total MUFAs", "Total PUFAs", "Total n3", "Total n6", "n6:n3 Ratio", "Total lipids"
)


cor_table$FAname <- label_mapping[cor_table$FA]
# Convertir FAname a factor y asignar el orden definido
cor_table$FAname <- factor(cor_table$FAname, levels = rev(custom_order))

# Crear el gráfico con ggplot
ggplot(cor_table, aes(x = Gene, y = FAname, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  geom_text(aes(label = ifelse(is.na(cor), "", paste0(sprintf("%.2f", cor), signif))), size = 3) +
  labs(title = "",
       x = "Gene", y = "Fatty Acid") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, face = "italic"),
        plot.title = element_text(hjust = 0.5))

#LO MISMO PARA AC GRASO DIETA Y EXP GENICA

# 1. Renombrar columnas de df_dieta
df_dieta_renamed <- df_dieta_exp %>%
  rename_with(~ paste0(., "_dieta"), .cols = -DIETA)

# 2. Expandir df_dieta (sin réplica) para empatar con 3 réplicas
df_dieta_expanded <- df_dieta_renamed %>%
  slice(rep(1:n(), each = 3)) %>%
  mutate(REPLICA = rep(1:3, times = nrow(df_dieta_renamed)))  # Agrega réplica ficticia

# 3. Unir df_dieta_expanded con df_exp_summary para tener expresión génica y ácidos grasos de la dieta
combined_data_dieta <- df_dieta_expanded %>%
  inner_join(df_exp_summary, by = c("DIETA", "REPLICA"))

# Verificar la estructura del dataframe combinado
glimpse(combined_data_dieta)

# 4. Preparar la matriz de correlación para los ácidos grasos de la dieta
acidos_dieta <- grep("_dieta$", names(combined_data_dieta), value = TRUE)
acidos_dieta <- acidos_dieta[sapply(combined_data_dieta[acidos_dieta], is.numeric)]

# Extraer expresión génica
expresion_genes <- c("elovl2", "elovl5")

# Crear la matriz de correlación
cor_matrix_dieta <- map_dfr(
  expresion_genes,
  function(gene) {
    map_dbl(acidos_dieta, ~ safe_cor(combined_data_dieta[[.x]], combined_data_dieta[[gene]])) %>%
      setNames(acidos_dieta) %>%
      as.data.frame() %>%
      mutate(Gene = gene)
  }
)

# Reformatear la matriz de correlación
cor_long_dieta <- cor_matrix_dieta %>%
  pivot_longer(-Gene, names_to = "Variable", values_to = "Correlation") %>%
  mutate(Variable = gsub("_dieta$", "", Variable))

# 5. Crear función para obtener correlación y p-valor
get_cor_p <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (
    all(is.na(x)) || all(is.na(y)) ||
    is.na(sd(x, na.rm = TRUE)) || is.na(sd(y, na.rm = TRUE)) ||
    sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0
  ) {
    return(tibble(cor = NA, p = NA))
  }
  test <- cor.test(x, y, method = "pearson")
  tibble(cor = test$estimate, p = test$p.value)
}

# Aplicar a cada combinación de ácidos grasos de la dieta y expresión génica
cor_table_dieta <- expand.grid(Gene = expresion_genes, FA = acidos_dieta, stringsAsFactors = FALSE) %>%
  mutate(
    result = map2(Gene, FA, ~ get_cor_p(combined_data_dieta[[.y]], combined_data_dieta[[.x]]))
  ) %>%
  unnest(cols = result) %>%
  mutate(FA_clean = gsub("_dieta$", "", FA))

# Ajustar p-valores (a nivel global para todas comparaciones)
cor_table_dieta <- cor_table_dieta %>%
  mutate(
    p_adj = p.adjust(p, method = "bonferroni"),  # Cambia a "bonferroni" si preferís
    signif = case_when(
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# 6. Crear el gráfico de correlación entre ácidos grasos de la dieta y expresión génica
ggplot(cor_table_dieta, aes(x = Gene, y = FA_clean, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  geom_text(aes(label = ifelse(is.na(cor), "", paste0(sprintf("%.2f", cor), signif))), size = 3) +
  labs(title = "Correlation between Diet Fatty Acids and elovl2/elovl5 Expression",
       x = "Gene", y = "Fatty Acid") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, face = "italic"),
        plot.title = element_text(hjust = 0.5))

# 1. Definir el mapeo de los nombres de los ácidos grasos de la dieta
label_mapping_dieta <- c(
  "C14:0 Acido Mirístico_dieta" = "Myristic Acid",
  "C15:0 Acido Pentadecanoico_dieta" = "Pentadecanoic Acid",
  "C16:0 Acido Palmítico_dieta" = "Palmitic Acid",
  "C17:0 Acido Heptadecanoico_dieta" = "Heptadecanoic Acid",
  "C18:0 Acido Esteárico_dieta" = "Stearic Acid",
  "C21:0 Acido Heneicosanoico_dieta" = "Heneicosanoic Acid",
  "C22:0 Acido Behénico_dieta" = "Behenic Acid",
  "C23:0 Acido Tricosanoico_dieta" = "Tricosanoic Acid",
  "C24:0 Acido Lignocérico_dieta" = "Lignoceric Acid",
  
  # Ácidos Grasos Monoinsaturados (MUFAs)
  "C16:1 Acido Palmitoléico_dieta" = "Palmitoleic Acid",
  "C17:1 Acido Cis-10-Heptadecenoico_dieta" = "Heptadecenoic Acid",
  "C18:1n9t Acido Elaidico_dieta" = "Elaidic Acid",
  "C18:1n9c Acido Oléico_dieta" = "Oleic Acid",
  "C20:1n9 Acido Cis-11-Eicosenoico_dieta" = "Eicosenoic Acid",
  
  # Ácidos Grasos Poliinsaturados (PUFAs)
  "C18:2n6c Acido Linoléico_dieta" = "18:2 n6 Linoleic Acid",
  "C18:3n6 Acido g-linolénico_dieta" = "Gamma Linoleic Acid",
  "C18:3n6 Acido g-linolénico_dieta" = "18:3 n6 Linolenic Acid",
  "C18:3n3 Acido Linolénico_dieta" = "18:3 n3 Linolenic Acid",
  "C20:2 Acido Cis-11,14-Eicosadienoico_dieta" = "Eicosadienoic Acid",
  "C20:4n6 Acido Araquidónico_dieta" = "Arachidonic Acid",
  "C20:5n3 Acido Cis-5,8,11,14,17-Eicosapentaenoico_dieta" = "Eicosapentaenoic Acid",
  "C22:2 Acido Cis-13,16-Docosadienoico_dieta" = "Docosadienoic Acid",
  "C22:6n3 Acido Cis-4,7,10,13,16,19-Docosahexaenoico_dieta" = "Docosahexaenoic Acid",
  
  # Totales y relaciones
  "Total de SAFAs_dieta" = "Total SFAs",
  "Total de MUFAs_dieta" = "Total MUFAs",
  "Total de PUFAs_dieta" = "Total PUFAs",
  "total_n3_dieta" = "Total n3",
  "total_n6_dieta" = "Total n6",
  "n6_n3_dieta" = "n6:n3 Ratio",
  "% Lipidos totales_dieta" = "Total lipids"
)

# 2. Definir el orden de los ácidos grasos de la dieta, de más saturado a menos saturado
custom_order_dieta <- c(
  # Saturados (SAFAs)
  "Myristic Acid", "Pentadecanoic Acid", "Palmitic Acid", "Heptadecanoic Acid", "Stearic Acid", 
  "Heneicosanoic Acid", "Behenic Acid", "Tricosanoic Acid", "Lignoceric Acid",
  
  # Monoinsaturados (MUFAs)
  "Palmitoleic Acid", "Heptadecenoic Acid", "Elaidic Acid", "Oleic Acid", "Eicosenoic Acid",
  
  # Poliinsaturados (PUFAs)
  "18:2 n6 Linoleic Acid", "Gamma Linoleic Acid", "18:3 n6 Linolenic Acid", "18:3 n3 Linolenic Acid",
  "Eicosadienoic Acid", "Arachidonic Acid", "Eicosapentaenoic Acid", "Docosadienoic Acid", "Docosahexaenoic Acid",
  
  # Totales y ratio
  "Total SFAs", "Total MUFAs", "Total PUFAs", "Total n3", "Total n6", "n6:n3 Ratio", "Total lipids"
)

# 3. Asignar los nombres a la tabla de correlaciones con el mapeo
cor_table_dieta$FAname <- label_mapping_dieta[cor_table_dieta$FA]

# Convertir FAname a factor y asignar el orden definido
cor_table_dieta$FAname <- factor(cor_table_dieta$FAname, levels = rev(custom_order_dieta))

# 4. Crear el gráfico con ggplot
ggplot(cor_table_dieta, aes(x = Gene, y = FAname, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  geom_text(aes(label = ifelse(is.na(cor), "", paste0(sprintf("%.2f", cor), signif))), size = 3) +
  labs(title = "Correlation between Diet Fatty Acids and elovl2/elovl5 Expression",
       x = "Gene", y = "Fatty Acid") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, face = "italic"),
        plot.title = element_text(hjust = 0.5))

# Crear el gráfico para dieta
plot_dieta <- ggplot(cor_table_dieta, aes(x = Gene, y = FAname, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  geom_text(aes(label = ifelse(is.na(cor), "", paste0(sprintf("%.2f", cor), signif))), size = 3) +
  labs(title = "", x = "Gene", y = "Fatty Acid") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, face = "italic"),
        plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none")


# Crear el gráfico para hígado
plot_higado <- ggplot(cor_table, aes(x = Gene, y = FAname, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  geom_text(aes(label = ifelse(is.na(cor), "", paste0(sprintf("%.2f", cor), signif))), size = 3) +
  labs(title = "", x = "Gene", y = "Fatty Acid") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, face = "italic"),
        plot.title = element_text(hjust = 0.5))

# Combinar ambos gráficos en una sola imagen con dos columnas
final_plot <- plot_dieta | plot_higado

# Mostrar el gráfico final
final_plot