#!/bin/bash
# En este script implemento el pipeline completo de análisis de microbiota fecal
# para los niños con muestras en t0 y t12 usando QIIME2, desde QC hasta
# abundancia diferencial (ANCOM y ANCOM-BC) siguiendo las prácticas del máster
# y la documentación oficial de QIIME2.

# ---------------------------------------------------------------------------
# PREFILTRADO Y CONTROL DE CALIDAD (FASTQC / MULTIQC)
# ---------------------------------------------------------------------------

# Primero realizo el FastQC de todos los archivos FASTQ comprimidos.
# FastQC calcula métricas básicas de calidad por posición de base, contenido
# de GC, presencia de adaptadores, etc., lo que me permite inspeccionar si la
# calidad cae a partir de cierta base y decidir el truncado para DADA2. 
fastqc *.gz

# MultiQC agrega en un solo informe los resultados de FastQC de todas las
# muestras, facilitando ver patrones globales y detectar muestras problemáticas
# sin tener que abrir informe por informe. 
multiqc ./

# ---------------------------------------------------------------------------
# METADATA Y ACTIVACIÓN DEL ENTORNO QIIME2
# ---------------------------------------------------------------------------

# Aquí edito manualmente el manifest.tsv para asegurar que cada sample-id
# apunta correctamente al archivo FASTQ correspondiente, requisito para el
# formato SingleEndFastqManifestPhred33V2. 
# (Este paso lo hago con un editor de texto, no hay comando automático.)

# Activo el entorno de QIIME2 en conda para tener disponibles todos los plugins
# necesarios (q2-diversity, q2-dada2, q2-phylogeny, q2-composition, etc.).
conda activate qiime2

# Visualizo los metadatos en formato .qzv para comprobar que las columnas
# (sample-id, anti-exposed-all, delivery, month, sex, diet, studyid, etc.)
# están correctamente formateadas y tipadas (numéricas vs categóricas). 
qiime metadata tabulate \
  --m-input-file metadata.tsv \
  --o-visualization metadata.qzv

# ---------------------------------------------------------------------------
# IMPORTACIÓN DE SECUENCIAS Y RESUMEN DE DEMUX
# ---------------------------------------------------------------------------

# Importo las secuencias single-end a QIIME2 usando un manifest con rutas
# absolutas/relativas y calidades Phred33. El resultado es un artefacto
# SampleData[SequencesWithQuality] que usaré como entrada de DADA2. 
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path ./manifest.tsv \
  --output-path ./demux_seqs.qza

# Genero un resumen interactivo de las secuencias demultiplexadas para evaluar
# la distribución de longitudes y la calidad por posición; esto me guía para
# elegir el parámetro --p-trunc-len de DADA2. 
qiime demux summarize \
  --i-data ./demux_seqs.qza \
  --o-visualization ./demux_seqs.qzv

# ---------------------------------------------------------------------------
# DENOISING CON DADA2 (ASV) Y ESTADÍSTICAS
# ---------------------------------------------------------------------------

# Aplico DADA2 en modo single-end para denoising, corrección de errores y
# eliminación de quimeras. Este método infiere variantes de secuencia de
# amplicón (ASV) exactas, corrige errores de secuenciación y construye una
# tabla de recuentos ASV x muestra. 
# --p-trunc-len 145 corta las lecturas a 145 nt (según la caída de calidad),
# y --p-n-threads 0 utiliza todos los hilos disponibles.
qiime dada2 denoise-single \
  --i-demultiplexed-seqs ./demux_seqs.qza \
  --p-trunc-len 145 \
  --p-n-threads 0 \
  --o-table ./dada2_table.qza \
  --o-representative-sequences ./dada2_rep_set.qza \
  --o-denoising-stats ./dada2_stats.qza

# Las estadísticas de denoising me permiten ver, por muestra, cuántas lecturas
# pasan los filtros, cuántas se consideran quimeras, y el porcentaje de lecturas
# no quiméricas finales; así evalúo si el filtrado es razonable. 
qiime metadata tabulate \
  --m-input-file ./dada2_stats.qza \
  --o-visualization ./dada2_stats.qzv

# Resumo la tabla de frecuencias (ASV x muestra) para ver número de muestras,
# número de ASV, frecuencia mínima/máxima por muestra y distribución de recuentos,
# información clave para elegir la profundidad de rarefacción. 
qiime feature-table summarize \
  --i-table ./dada2_table.qza \
  --m-sample-metadata-file ./metadata.tsv \
  --o-visualization ./dada2_table.qzv

# Visualizo las secuencias representativas para comprobar la longitud de las ASV
# y/o lanzar BLAST contra NCBI si quiero validar alguna asignación. 
qiime feature-table tabulate-seqs \
  --i-data dada2_rep_set.qza \
  --o-visualization dada2_rep_set.qzv

# ---------------------------------------------------------------------------
# NORMALIZACIÓN POR RAREFACCIÓN Y ELECCIÓN DE PROFUNDIDAD
# ---------------------------------------------------------------------------

# Genero curvas de rarefacción alfa para explorar cómo crece la riqueza
# observada (observed_features) en función de la profundidad de secuenciación.
# Esto me permite ver dónde se estabilizan las curvas y estimar una
# profundidad de rarefacción que mantenga el máximo de muestras posible. 
qiime diversity alpha-rarefaction \
  --i-table ./dada2_table.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./alpha_rarefaction_curves.qzv \
  --p-min-depth 500 \
  --p-max-depth 15000

# A partir de estas curvas y del resumen de la tabla, decido una profundidad de
# rarefacción (en la práctica en torno a 8375). La rarefacción es un submuestreo
# sin reemplazo: primero elimina muestras por debajo de la profundidad elegida y
# luego recorta todas las demás a esa profundidad para igualar la profundidad de
# secuenciación entre muestras. 

# (La sección siguiente con rarefy + alpha/beta explícitos es ilustrativa; en la
# práctica, `core-metrics` ya hace esta rarefacción internamente.)

# Rarefacción explícita de la tabla para obtener una tabla normalizada
# (mismas lecturas por muestra) sobre la que calculo diversidad alfa/beta.
qiime feature-table rarefy \
  --i-table ./dada2_table.qza \
  --p-sampling-depth 8375 \
  --o-rarefied-table ./rar_table.qza

# Resumen de la tabla rarefactada para comprobar cuántas muestras se han perdido
# y verificar que todas las restantes tienen exactamente la misma profundidad.
qiime feature-table summarize \
  --i-table ./rar_table.qza \
  --m-sample-metadata-file ./metadata.tsv \
  --o-visualization ./rar_table.qzv

# ---------------------------------------------------------------------------
# DIVERSIDAD ALFA (OBSERVED, SHANNON, EVENNESS)
# ---------------------------------------------------------------------------

# La diversidad alfa mide la diversidad dentro de cada muestra. Aquí calculo:
# - observed_features: riqueza (número de ASV distintos).
# - shannon: índice de Shannon, que combina riqueza y equidad.
# - pielou_e: evenness de Pielou, equidad de la distribución. 

# Observed_features: cuenta cuántas ASV distintas hay en cada muestra.
qiime diversity alpha \
  --i-table ./rar_table.qza \
  --p-metric observed_features \
  --o-alpha-diversity ./obs_featur.qza

# Tabulo los valores de riqueza por muestra para su inspección y posterior
# asociación con covariables (Kruskal-Wallis en alpha-group-significance). 
qiime metadata tabulate \
  --m-input-file ./obs_featur.qza \
  --o-visualization ./obs_featur.qzv

# Índice de Shannon, sensible tanto al número de ASV como a su abundancia
# relativa; aumenta cuando la comunidad es más rica y más homogénea. 
qiime diversity alpha \
  --i-table ./rar_table.qza \
  --p-metric shannon \
  --o-alpha-diversity ./shannon.qza

qiime metadata tabulate \
  --m-input-file ./shannon.qza \
  --o-visualization ./shannon.qzv

# Evenness de Pielou, que normaliza Shannon por la riqueza máxima posible y
# mide cuán uniforme es la distribución de abundancias entre ASV. 
qiime diversity alpha \
  --i-table ./rar_table.qza \
  --p-metric pielou_e \
  --o-alpha-diversity ./pielou_e.qza

qiime metadata tabulate \
  --m-input-file ./pielou_e.qza \
  --o-visualization ./pielou_e.qzv

# ---------------------------------------------------------------------------
# DIVERSIDAD BETA (BRAY-CURTIS, JACCARD) Y ORDINACIÓN (PCoA + EMPEROR)
# ---------------------------------------------------------------------------

# La diversidad beta cuantifica diferencias de composición entre muestras.
# Aquí uso:
# - Bray-Curtis: distancia basada en abundancias relativas.
# - Jaccard: disimilitud basada en presencia/ausencia. 

# Bray-Curtis: genera una matriz de distancias entre todas las parejas de
# muestras según sus abundancias de ASV; valores cercanos a 1 indican comunidades
# muy diferentes. [cite:7]
qiime diversity beta \
  --i-table ./rar_table.qza \
  --p-metric braycurtis \
  --o-distance-matrix ./braycurtis.qza

# PCoA sobre la matriz de Bray-Curtis, que proyecta las distancias en un
# espacio de menor dimensión para visualizar patrones de agrupamiento global. 
qiime diversity pcoa \
  --i-distance-matrix ./braycurtis.qza \
  --o-pcoa ./braycurtis_pcoa.qza

# Emperor genera una representación interactiva 3D de la PCoA, donde coloreo las
# muestras por variables de metadatos (delivery, month, etc.) para explorar
# posibles separaciones de grupos. 
qiime emperor plot \
  --i-pcoa ./braycurtis_pcoa.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./braycurtis_emperor.qzv

# Jaccard: calcula una matriz de disimilitudes basadas solo en presencia/ausencia
# de ASV, ignorando las abundancias, útil para ver cambios en composición sin
# ponderar por recuentos. 
qiime diversity beta \
  --i-table ./rar_table.qza \
  --p-metric jaccard \
  --o-distance-matrix ./jaccard.qza

qiime diversity pcoa \
  --i-distance-matrix ./jaccard.qza \
  --o-pcoa ./jaccard_pcoa.qza

qiime emperor plot \
  --i-pcoa ./jaccard_pcoa.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./jaccard_emperor.qzv

# ---------------------------------------------------------------------------
# DIVERSIDAD: PIPELINE INTEGRADO CORE-METRICS (NO FILOGENÉTICO)
# ---------------------------------------------------------------------------

# `core-metrics` aplica en un solo paso la rarefacción y calcula un conjunto
# de índices alfa (observed_features, Shannon, evenness) y beta (Bray-Curtis y
# Jaccard), generando también las PCoA y los Emperor plots asociados. 

qiime diversity core-metrics \
  --i-table ./dada2_table.qza \
  --m-metadata-file ./metadata.tsv \
  --p-sampling-depth 8375 \
  --output-dir ./core-metrics-results

# Para cada índice alfa genero una prueba de significación por grupos con
# `alpha-group-significance`, que aplica pruebas no paramétricas tipo
# Kruskal-Wallis para comparar la distribución del índice entre categorías de
# metadatos (delivery, month, etc.). 
qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-results/observed_features_vector.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./core-metrics-results/observed_features_statistics.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-results/evenness_vector.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./core-metrics-results/pielou_e_statistics.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-results/shannon_vector.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./core-metrics-results/shannon_statistics.qzv

# ---------------------------------------------------------------------------
# DIVERSIDAD BETA: SIGNIFICACIÓN ESTADÍSTICA (PERMANOVA / ANOSIM / PERMDISP)
# ---------------------------------------------------------------------------

# `beta-group-significance` contrasta si las distancias intra-grupo son menores
# que las distancias entre grupos usando tests de permutación (PERMANOVA por
# defecto, también ANOSIM o PERMDISP). Aquí pruebo si la composición difiere
# entre grupos definidos por anti-exposed-all, delivery, diet, month. 

# Anti-exposed-all (sí/no antibióticos). PERMANOVA evalúa si el centroide de
# cada grupo en el espacio de distancias Bray-Curtis/Jaccard difiere de forma
# significativa al comparar la suma de cuadrados entre grupos vs dentro de
# grupos, utilizando permutaciones para obtener el valor p. 
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column anti-exposed-all \
  --o-visualization core-metrics-results/bray_curtis_anti-exposed-all_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column anti-exposed-all \
  --o-visualization core-metrics-results/jaccard_anti-exposed-all_significance.qzv \
  --p-pairwise

# Lo mismo para el tipo de parto (delivery): comparo la estructura de la
# microbiota entre nacidos por cesárea y por parto vaginal. 
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column delivery \
  --o-visualization core-metrics-results/bray_curtis_delivery_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column delivery \
  --o-visualization core-metrics-results/jaccard_delivery_significance.qzv \
  --p-pairwise

# Comparo también las comunidades según el tipo de dieta (diet) usando las
# mismas matrices de distancia, siempre con permutaciones y comparaciones
# por pares cuando hay más de dos categorías.
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column diet \
  --o-visualization core-metrics-results/bray_curtis_diet_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column diet \
  --o-visualization core-metrics-results/jaccard_diet_significance.qzv \
  --p-pairwise

# Para utilizar month como variable categórica (0 vs 12) en un test de grupos,
# he creado un metadata_fixed.tsv donde declaro explícitamente month como
# categórica en la fila de tipos de QIIME2; esto es importante para que el
# plugin interprete correctamente la variable. 
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata_fixed.tsv \
  --m-metadata-column month \
  --o-visualization core-metrics-results/bray_curtis_month_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file metadata_fixed.tsv \
  --m-metadata-column month \
  --o-visualization core-metrics-results/jaccard_month_significance.qzv \
  --p-pairwise

# ---------------------------------------------------------------------------
# FILOGENIA (ÁRBOL DE NOVO) Y DIVERSIDAD FILOGENÉTICA (FAITH PD, UNIFRAC)
# ---------------------------------------------------------------------------

# Aquí construyo un árbol filogenético de novo a partir de las secuencias
# representativas, usando el pipeline `align-to-tree-mafft-fasttree` que:
# 1) alinea con MAFFT, 2) enmascara regiones ambiguas, 3) infiere el árbol
# con FastTree, 4) lo enrraiza en el punto medio. 
# Este árbol es necesario para métricas como Faith PD y UniFrac. 
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ./dada2_rep_set.qza \
  --o-alignment ./aligned_dada2_rep_seq.qza \
  --o-masked-alignment ./masked_aligned_dada2_rep_seq.qza \
  --o-tree ./unrooted_tree.qza \
  --o-rooted-tree ./rooted_tree.qza

# Con el árbol enraizado y la tabla sin rarefacción, `core-metrics-phylogenetic`
# calcula tanto métricas no filogenéticas (observed, Shannon, evenness,
# Bray-Curtis, Jaccard) como métricas filogenéticas (Faith PD, UniFrac
# ponderado y no ponderado), integrando la historia evolutiva de las ASV en la
# estimación de la diversidad. 
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ./rooted_tree.qza \
  --i-table ./dada2_table.qza \
  --p-sampling-depth 8375 \
  --m-metadata-file ./metadata.tsv \
  --output-dir core-metrics-results

# Repito las pruebas de significación de diversidad alfa sobre los índices
# generados por `core-metrics-phylogenetic`, de forma análoga a antes.
qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-results/observed_features_vector.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./core-metrics-results/observed_features_statistics.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-results/evenness_vector.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./core-metrics-results/pielou_e_statistics.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-results/shannon_vector.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./core-metrics-results/shannon_statistics.qzv

# Y vuelvo a aplicar beta-group-significance, ahora sobre matrices filogenéticas
# (UniFrac) y no filogenéticas, para estudiar si la estructura filogenética de
# las comunidades difiere entre grupos clínicos. 
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column anti-exposed-all \
  --o-visualization ./core-metrics-results/bray_curtis_anti-exposed-all_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column anti-exposed-all \
  --o-visualization ./core-metrics-results/jaccard_anti-exposed-all_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column delivery \
  --o-visualization ./core-metrics-results/bray_curtis_delivery_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column delivery \
  --o-visualization ./core-metrics-results/jaccard_delivery_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column diet \
  --o-visualization ./core-metrics-results/bray_curtis_diet_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column diet \
  --o-visualization ./core-metrics-results/jaccard_diet_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column month \
  --o-visualization ./core-metrics-results/bray_curtis_month_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column month \
  --o-visualization ./core-metrics-results/jaccard_month_significance.qzv \
  --p-pairwise

# ---------------------------------------------------------------------------
# CLASIFICACIÓN TAXONÓMICA (NAIVE BAYES) Y GRÁFICOS DE BARRAS
# ---------------------------------------------------------------------------

# Clasifico las ASV con un clasificador Naive Bayes entrenado sobre Greengenes
# 13_8 recortado a la región V4 (515F-806R). `classify-sklearn` asigna a cada
# ASV una jerarquía taxonómica y una confianza de asignación. 
qiime feature-classifier classify-sklearn \
  --i-reads ./dada2_rep_set.qza \
  --i-classifier ./gg-13-8-99-515-806-nb-classifier.qza \
  --o-classification ./taxonomy.qza

# Visualizo la tabla taxonómica para inspeccionar las asignaciones y sus
# confidencias, y para identificar ASV de interés. 
qiime metadata tabulate \
  --m-input-file ./taxonomy.qza \
  --o-visualization ./taxonomy.qzv

# Genero un barplot taxonómico a partir de la tabla sin rarefacción para ver la
# composición relativa de filos, clases, etc., por muestra/grupo de metadatos.
qiime taxa barplot \
  --i-table ./dada2_table.qza \
  --i-taxonomy ./taxonomy.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./taxa_barplot.qzv

# Repito el barplot utilizando la tabla rarefactada por si quiero comparar
# resultados normalizados por profundidad de secuenciación.
qiime taxa barplot \
  --i-table ./rar_table.qza \
  --i-taxonomy ./taxonomy.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./taxa_barplot_rar.qzv

# Opcionalmente filtro muestras de baja profundidad (ej. >2020 lecturas) para
# generar barplots centrados en muestras de buena cobertura.
qiime feature-table filter-samples \
  --i-table ./dada2_table.qza \
  --p-min-frequency 2020 \
  --o-filtered-table ./table_2k.qza

qiime taxa barplot \
  --i-table ./table_2k.qza \
  --i-taxonomy ./taxonomy.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./taxa_barplot_2k.qzv

# ---------------------------------------------------------------------------
# ABUNDANCIA DIFERENCIAL: ANCOM (A NIVEL ASV Y COLAPSADO POR TAXONOMÍA)
# ---------------------------------------------------------------------------

# ANCOM es un test de abundancia diferencial diseñado para datos composicionales,
# que evita usar directamente tests clásicos tipo t-test/ANOVA sobre proporciones.
# Antes debo eliminar ceros añadiendo un pseudoconteo, porque ANCOM usa
# transformaciones logarítmicas que no toleran valores cero.
qiime composition add-pseudocount \
  --i-table ./dada2_table.qza \
  --o-composition-table ./table2k_comp.qza

# Aplico ANCOM para distintas variables de metadatos (anti-exposed-all,
# delivery, diet, month). Internamente, ANCOM compara log-cocientes de cada
# ASV frente a todos los demás y contabiliza cuántos cocientes muestran
# diferencias significativas (estadístico W); ASV con W alto se consideran
# diferencialmente abundantes entre grupos. 
qiime composition ancom \
  --i-table ./table2k_comp.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column anti-exposed-all \
  --o-visualization ./ancom_anti-exposed-all.qzv

qiime composition ancom \
  --i-table ./table2k_comp.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column delivery \
  --o-visualization ./ancom_delivery.qzv

qiime composition ancom \
  --i-table ./table2k_comp.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column diet \
  --o-visualization ./ancom_diet.qzv

qiime composition ancom \
  --i-table ./table2k_comp.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column month \
  --o-visualization ./ancom_month.qzv

# Para interpretar mejor los ASV detectados por ANCOM, colapso la tabla de
# abundancias al nivel taxonómico 6 (género) y repito el análisis. Esto produce
# resultados más legibles, pero introduce la incertidumbre derivada de la
# clasificación taxonómica.
qiime taxa collapse \
  --i-table ./dada2_table.qza \
  --i-taxonomy ./taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table ./tax-table.qza

qiime composition add-pseudocount \
  --i-table ./tax-table.qza \
  --o-composition-table ./comp_tax-table.qza

qiime composition ancom \
  --i-table ./comp_tax-table.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column month \
  --o-visualization ./ancom_ASVtaxo_month.qzv

# ---------------------------------------------------------------------------
# ABUNDANCIA DIFERENCIAL: ANCOM-BC (CORRECCIÓN DE SESGOS)
# ---------------------------------------------------------------------------

# ANCOM-BC extiende ANCOM para corregir sesgos de fracción de muestreo
# desigual entre muestras y eficiencias diferenciales de secuenciación entre
# taxones; implementa un modelo de regresión sobre abundancias absolutas
# estimadas y permite fórmulas complejas (covariables). 

# En estos bloques uso distintos archivos de metadata (metadata2.tsv,
# metadata3.tsv) para adaptar los nombres de las variables a las fórmulas
# (`antiexposedall`, `diet`, `month`) y defino un nivel de referencia para cada
# factor (por ejemplo, 'y' para exposición a antibióticos, 'Vaginal' para
# parto, 'bd' para dieta). 

qiime composition ancombc \
  --i-table ./dada2_table.qza \
  --m-metadata-file ./metadata2.tsv \
  --p-formula 'antiexposedall' \
  --p-p-adj-method 'holm' \
  --p-prv-cut 0.1 \
  --p-no-conserve \
  --p-alpha 0.05 \
  --p-reference-levels antiexposedall::'y' \
  --o-differentials ./ancombc_anti-exposed-all_ASV.qza

# Tabulo los coeficientes de ANCOM-BC (log-fold changes estimados, valores p
# ajustados) y genero barplots de abundancias diferenciales para los ASV
# significativos según un umbral de significación. 
qiime composition tabulate \
  --i-data ancombc_anti-exposed-all_ASV.qza \
  --o-visualization ancombc_anti-exposed-all_ASV.qzv

qiime composition da-barplot \
  --i-data ancombc_anti-exposed-all_ASV.qza \
  --p-significance-threshold 0.05 \
  --o-visualization plot_ancombc_ASV.qzv

qiime composition ancombc \
  --i-table ./dada2_table.qza \
  --m-metadata-file ./metadata.tsv \
  --p-formula 'delivery' \
  --p-p-adj-method 'holm' \
  --p-prv-cut 0.1 \
  --p-no-conserve \
  --p-alpha 0.05 \
  --p-reference-levels delivery::'Vaginal' \
  --o-differentials ./ancombc_delivery_ASV.qza

qiime composition tabulate \
  --i-data ancombc_delivery_ASV.qza \
  --o-visualization ancombc_delivery_ASV.qzv

qiime composition da-barplot \
  --i-data ancombc_delivery_ASV.qza \
  --p-significance-threshold 0.05 \
  --o-visualization plot_ancombc_delivery_ASV.qzv

qiime composition ancombc \
  --i-table ./dada2_table.qza \
  --m-metadata-file ./metadata2.tsv \
  --p-formula 'diet' \
  --p-p-adj-method 'holm' \
  --p-prv-cut 0.1 \
  --p-no-conserve \
  --p-alpha 0.05 \
  --p-reference-levels diet::'bd' \
  --o-differentials ./ancombc_diet_ASV.qza

qiime composition tabulate \
  --i-data ancombc_diet_ASV.qza \
  --o-visualization ancombc_diet_ASV.qzv

qiime composition da-barplot \
  --i-data ancombc_diet_ASV.qza \
  --p-significance-threshold 0.05 \
  --o-visualization plot_ancombc_diet_ASV.qzv

qiime composition ancombc \
  --i-table ./dada2_table.qza \
  --m-metadata-file ./metadata3.tsv \
  --p-formula 'month' \
  --p-p-adj-method 'holm' \
  --p-prv-cut 0.1 \
  --p-no-conserve \
  --p-alpha 0.05 \
  --p-reference-levels month::'0' \
  --o-differentials ./ancombc_month_ASV.qza

qiime composition tabulate \
  --i-data ancombc_month_ASV.qza \
  --o-visualization ancombc_month_ASV.qzv

qiime composition da-barplot \
  --i-data ancombc_month_ASV.qza \
  --p-significance-threshold 0.05 \
  --o-visualization plot_ancombc_month_ASV.qzv
