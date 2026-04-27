#!/bin/bash
# En este script parto de la tabla de ASV completa (t0+t12) y filtro para
# quedarme solo con las muestras de tiempo 0, repitiendo los análisis de
# diversidad, filogenia, taxonomía y abundancia diferencial pero restringidos
# a t0, lo que permite interpretar la microbiota neonatal sin el efecto del
# seguimiento longitudinal. [cite:2][cite:7]

# ---------------------------------------------------------------------------
# FILTRADO DE LA TABLA A MUESTRAS T0
# ---------------------------------------------------------------------------

# Filtro la tabla de ASV (dada2_table.qza) seleccionando únicamente las
# muestras con month = '0' usando una expresión SQL-like sobre el metadata.
# Esto genera una nueva tabla table_t0.qza con solo muestras t0. [cite:7]
qiime feature-table filter-samples \
  --i-table dada2_table.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[month]='0'" \
  --o-filtered-table table_t0.qza

# Resumo la nueva tabla filtrada para conocer el número de muestras t0, el
# rango de frecuencias y decidir posteriormente la profundidad de rarefacción. [cite:7]
qiime feature-table summarize \
  --i-table table_t0.qza \
  --m-sample-metadata-file metadata.tsv \
  --o-visualization table_t0.qzv

# ---------------------------------------------------------------------------
# NORMALIZACIÓN POR RAREFACCIÓN Y CORE-METRICS (T0)
# ---------------------------------------------------------------------------

# Genero curvas de rarefacción alfa solo para las muestras t0, con el mismo
# rango de profundidad que antes, para confirmar que la profundidad elegida
# sigue siendo razonable en este subconjunto. [cite:7]
qiime diversity alpha-rarefaction \
  --i-table ./table_t0.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./alpha_rarefaction_curves.qzv \
  --p-min-depth 500 \
  --p-max-depth 15000

# A partir de estas curvas y del resumen, concluyo que una profundidad de
# rarefacción en torno a 8375 sigue dando una buena cobertura para t0.

# Aplico `core-metrics` sobre la tabla de t0 para obtener índices alfa (observed,
# Shannon, evenness) y beta (Bray-Curtis, Jaccard) ya rarefactados.
qiime diversity core-metrics \
  --i-table ./table_t0.qza \
  --m-metadata-file ./metadata.tsv \
  --p-sampling-depth 8375 \
  --output-dir ./core-metrics-results

# ---------------------------------------------------------------------------
# DIVERSIDAD ALFA Y TESTS POR GRUPO (T0)
# ---------------------------------------------------------------------------

# Pruebas de significación de diversidad alfa para t0 usando Kruskal-Wallis
# (alpha-group-significance) sobre observed_features, evenness y Shannon,
# comparando, por ejemplo, por tipo de parto, dieta o sexo. [cite:7][cite:21]
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
# DIVERSIDAD BETA Y PERMANOVA (T0)
# ---------------------------------------------------------------------------

# Tests de significación de diversidad beta (PERMANOVA/ANOSIM/PERMDISP según
# el parámetro method, por defecto PERMANOVA) aplicados ahora solo a t0. [cite:19]
# Aquí no uso anti-exposed-all ni month porque para t0 todos los sujetos están
# en el mismo nivel de esas variables.

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

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column sex \
  --o-visualization core-metrics-results/bray_curtis_sex_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column sex \
  --o-visualization core-metrics-results/jaccard_sex_significance.qzv \
  --p-pairwise

# ---------------------------------------------------------------------------
# FILTRADO DE REP-SEQS A T0 Y FILOGENIA
# ---------------------------------------------------------------------------

# Para asegurar coherencia entre tabla y árbol, filtro las secuencias
# representativas para retener solo las ASV presentes en table_t0.qza.
qiime feature-table filter-seqs \
  --i-data dada2_rep_set.qza \
  --i-table table_t0.qza \
  --o-filtered-data rep_set_t0.qza

# Construyo el árbol filogenético de novo únicamente para las ASV de t0, usando
# el mismo pipeline MAFFT + enmascarado + FastTree + enrraizamiento. [cite:8]
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ./rep_set_t0.qza \
  --o-alignment ./aligned_t0_rep_seq.qza \
  --o-masked-alignment ./masked_aligned_t0_rep_seq.qza \
  --o-tree ./unrooted_tree.qza \
  --o-rooted-tree ./rooted_tree.qza

# Calculo de nuevo las métricas de diversidad (incluyendo Faith PD y UniFrac)
# usando `core-metrics-phylogenetic`, ahora solo en t0. [cite:8][cite:25]
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ./rooted_tree.qza \
  --i-table ./table_t0.qza \
  --p-sampling-depth 8375 \
  --m-metadata-file ./metadata.tsv \
  --output-dir core-metrics-results

# Repito las pruebas de diversidad alfa en este contexto filogenético.
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

# Y las pruebas de diversidad beta (PERMANOVA) solo en t0.
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
  --m-metadata-column sex \
  --o-visualization ./core-metrics-results/bray_curtis_sex_significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column sex \
  --o-visualization ./core-metrics-results/jaccard_sex_significance.qzv \
  --p-pairwise

# ---------------------------------------------------------------------------
# CLASIFICACIÓN TAXONÓMICA (T0) Y TAXA BARPLOT
# ---------------------------------------------------------------------------

# Clasificación taxonómica utilizando solo las ASV de t0, para obtener una
# visión específica de la microbiota neonatal. [cite:8]
qiime feature-classifier classify-sklearn \
  --i-reads ./rep_set_t0.qza \
  --i-classifier ./gg-13-8-99-515-806-nb-classifier.qza \
  --o-classification ./taxonomy.qza

qiime metadata tabulate \
  --m-input-file ./taxonomy.qza \
  --o-visualization ./taxonomy.qzv

qiime taxa barplot \
  --i-table ./table_t0.qza \
  --i-taxonomy ./taxonomy.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./taxa_barplot.qzv

# ---------------------------------------------------------------------------
# ABUNDANCIA DIFERENCIAL EN T0: ANCOM Y ANCOM-BC
# ---------------------------------------------------------------------------

# Añadir un pseudoconteo para evitar ceros antes de ANCOM y ANCOM-BC. [cite:8][cite:22]
qiime composition add-pseudocount \
  --i-table ./table_t0.qza \
  --o-composition-table ./table2k_comp.qza

# ANCOM para delivery, diet y sex en t0.
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
  --m-metadata-column sex \
  --o-visualization ./ancom_sex.qzv

# Colapso la tabla a nivel de género (nivel 6) y repito ANCOM para diet, lo que
# facilita interpretar qué géneros cambian entre categorías de dieta. [cite:8]
qiime taxa collapse \
  --i-table ./table_t0.qza \
  --i-taxonomy ./taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table ./tax-table.qza

qiime composition add-pseudocount \
  --i-table ./tax-table.qza \
  --o-composition-table ./comp_tax-table.qza

qiime composition ancom \
  --i-table ./comp_tax-table.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column diet \
  --o-visualization ./ancom_ASVtaxo_month.qzv

qiime composition ancom \
  --i-table ./comp_tax-table.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column diet \
  --o-visualization ./ancom_ASVtaxo_diet.qzv

# ANCOM-BC para delivery, diet y sex en t0, con ajustes de FDR y selección de
# niveles de referencia, igual que en el script t0+t12. [cite:22][cite:30]
qiime composition ancombc \
  --i-table ./table_t0.qza \
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
  --i-table ./table_t0.qza \
  --m-metadata-file ./metadata.tsv \
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
  --i-table ./table_t0.qza \
  --m-metadata-file ./metadata.tsv \
  --p-formula 'sex' \
  --p-p-adj-method 'holm' \
  --p-prv-cut 0.1 \
  --p-no-conserve \
  --p-alpha 0.05 \
  --p-reference-levels sex::'Female' \
  --o-differentials ./ancombc_sex_ASV.qza

qiime composition tabulate \
  --i-data ancombc_sex_ASV.qza \
  --o-visualization ancombc_sex_ASV.qzv

qiime composition da-barplot \
  --i-data ancombc_sex_ASV.qza \
  --p-significance-threshold 0.05 \
  --o-visualization plot_ancombc_sex_ASV.qzv
