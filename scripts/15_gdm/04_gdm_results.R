# Summarize and visualize GDM results

# Load libraries
library(tidyverse)
library(vegan)
library(gdm)
library(janitor)
library(patchwork)

# Load palettes
source("scripts/00_palettes.R")


#### Import and format data ####
# Import (generalized alpha = 0.5) UniFrac distance matrix
load("output/dissim/lp2017-2019watercolumn18s_protists_unifrac_gen0pt5alpha.rda")

# Import melted protist ASV table
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes_function.tsv")
length(unique(asv_melt$asv_code))

# Import metadata
source("scripts/05_metadata.R")


#### Compute dissimilarity matrix ####
# All protists ----
# Compute dissimilarity matrix
asv_bray <- asv_melt %>%
  dplyr::select(lake_id, asv_code, relseqs) %>%
  pivot_wider(names_from = asv_code, values_from = relseqs, values_fill = 0) %>%
  arrange(lake_id) %>%
  column_to_rownames("lake_id") %>%
  vegdist(method = "bray")

# Phototrophs ----
# Recalculate relative sequence abundances for phototrophs
phototrophs_melt <- asv_melt %>%
  filter(trophic_group == "phototrophs") %>%
  dplyr::select(-relseqs)

phototrophs_nseqs <- phototrophs_melt %>%
  group_by(lake_id) %>%
  summarize(nseqs_total = sum(nseqs))

phototrophs_melt <- phototrophs_melt %>%
  left_join(phototrophs_nseqs, by = "lake_id") %>%
  mutate(relseqs = nseqs/nseqs_total) %>%
  dplyr::select(-nseqs_total)
sum(phototrophs_melt$relseqs)  # Should evaluate to the number of lakes

# Compute dissimilarity matrix
phototrophs_bray <- phototrophs_melt %>%
  dplyr::select(lake_id, asv_code, relseqs) %>%
  pivot_wider(names_from = asv_code, values_from = relseqs, values_fill = 0) %>%
  arrange(lake_id) %>%
  column_to_rownames("lake_id") %>%
  vegdist(method = "bray")

# Consumers ----
# Recalculate relative sequence abundances for consumers
consumers_melt <- asv_melt %>%
  filter(trophic_group == "consumers") %>%
  dplyr::select(-relseqs)

consumers_nseqs <- consumers_melt %>%
  group_by(lake_id) %>%
  summarize(nseqs_total = sum(nseqs))

consumers_melt <- consumers_melt %>%
  left_join(consumers_nseqs, by = "lake_id") %>%
  mutate(relseqs = nseqs/nseqs_total) %>%
  dplyr::select(-nseqs_total)
sum(consumers_melt$relseqs)  # Should evaluate to the number of lakes

# Compute dissimilarity matrix
consumers_bray <- consumers_melt %>%
  dplyr::select(lake_id, asv_code, relseqs) %>%
  pivot_wider(names_from = asv_code, values_from = relseqs, values_fill = 0) %>%
  arrange(lake_id) %>%
  column_to_rownames("lake_id") %>%
  vegdist(method = "bray")

# Mixotrophs ----
# Recalculate relative sequence abundances for mixotrophs
mixotrophs_melt <- asv_melt %>%
  filter(trophic_group == "mixotrophs") %>%
  dplyr::select(-relseqs)

mixotrophs_nseqs <- mixotrophs_melt %>%
  group_by(lake_id) %>%
  summarize(nseqs_total = sum(nseqs))

mixotrophs_melt <- mixotrophs_melt %>%
  left_join(mixotrophs_nseqs, by = "lake_id") %>%
  mutate(relseqs = nseqs/nseqs_total) %>%
  dplyr::select(-nseqs_total)
sum(mixotrophs_melt$relseqs)  # Should evaluate to the number of lakes

# Compute dissimilarity matrix
mixotrophs_bray <- mixotrophs_melt %>%
  dplyr::select(lake_id, asv_code, relseqs) %>%
  pivot_wider(names_from = asv_code, values_from = relseqs, values_fill = 0) %>%
  arrange(lake_id) %>%
  column_to_rownames("lake_id") %>%
  vegdist(method = "bray")

# Define function to calculate n sequences per lake to weight site-pairs
# (to account for biases in sampling effort)
countSeqs <- function(asv_melt) {
  asv_melt %>%
    group_by(lake_id) %>%
    dplyr::summarize(weights = sum(nseqs)) %>%
    ungroup() %>%
    arrange(lake_id) %>%
    as.data.frame() %>%
    return()
}

lakes_nseqs_all <- countSeqs(asv_melt)
lakes_nseqs_phototrophs <- countSeqs(phototrophs_melt)
lakes_nseqs_consumers <- countSeqs(consumers_melt)
lakes_nseqs_mixotrophs <- countSeqs(mixotrophs_melt)


#### Import GDM permutation test results ####
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_protists_bray_geography_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_protists_bray_morphometry_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_protists_bray_landuse_soil_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_protists_bray_climate_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_protists_bray_waterchem_100perm.rda")

load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_phylogdm_protists_genunifrac_geography_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_phylogdm_protists_genunifrac_morphometry_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_phylogdm_protists_genunifrac_landuse_soil_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_phylogdm_protists_genunifrac_climate_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_phylogdm_protists_genunifrac_waterchem_100perm.rda")

load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_phototrophs_bray_geography_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_phototrophs_bray_morphometry_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_phototrophs_bray_landuse_soil_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_phototrophs_bray_climate_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_phototrophs_bray_waterchem_100perm.rda")

load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_consumers_bray_geography_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_consumers_bray_morphometry_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_consumers_bray_landuse_soil_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_consumers_bray_climate_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_consumers_bray_waterchem_100perm.rda")

load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_mixotrophs_bray_geography_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_mixotrophs_bray_morphometry_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_mixotrophs_bray_landuse_soil_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_mixotrophs_bray_climate_100perm.rda")
load("output/gdm/gdm_perm/lp2017-2019watercolumn18s_gdm_mixotrophs_bray_waterchem_100perm.rda")


#### Summarize GDM permutation test results ####
# Define function to compute GDMs based on selected variables
computeGDM <- function(gdm, dist, lakes_nseqs, envtable) {
  gdm_signif <- gdm[[3]] %>%
    as_tibble(rownames = "predictor") %>%
    clean_names()
  
  model_signif <- gdm_signif %>%
    pivot_longer(!predictor, names_to = "model", values_to = "p_value") %>%
    mutate(model_number = case_when(grepl("full_model_", model) ~ as.numeric(str_remove(model, "full_model_")),
                                    model == "full_model" ~ 0)) %>%
    group_by(model, model_number) %>%
    summarize(max_p_value = max(p_value, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(max_p_value < 0.05) %>%
    slice_min(model_number) %>%
    pull(model)
  
  if (length(model_signif) == 1) {
    predictors_signif <- gdm_signif %>%
      dplyr::select(predictor, p_value = model_signif) %>%
      filter(!is.na(p_value)) %>%
      pull(predictor)
    
    if ("Geographic" %in% predictors_signif) {
      geodist <- TRUE
    } else if (!"Geographic" %in% predictors_signif) {
      geodist <- FALSE
    }
    
    predictors_signif <- predictors_signif[!predictors_signif %in% "Geographic"]
    
    # Format dissimilarity matrix
    dissim <- dist %>%
      as.matrix() %>%
      as_tibble(rownames = "lake_id") %>%
      arrange(lake_id)
    
    # Filter significant predictors in explanatory data
    if ("latitude" %in% names(envtable) | "longitude" %in% names(envtable)) {
      envtable <- envtable %>%
        dplyr::select(lake_id, where(is.numeric)) %>%
        dplyr::select(lake_id, latitude, longitude, all_of(predictors_signif)) %>%
        arrange(lake_id)
    } else {
      envtable <- latlong %>%
        left_join(envtable, by = "lake_id") %>%
        dplyr::select(lake_id, where(is.numeric)) %>%
        dplyr::select(lake_id, latitude, longitude, all_of(predictors_signif)) %>%
        arrange(lake_id)
    }
    
    # Format site-pairs
    sitepairs <- formatsitepair(bioData = dissim, bioFormat = 3,
                                abundance = TRUE, siteColumn = "lake_id",
                                XColumn = "longitude", YColumn = "latitude",
                                predData = envtable,
                                weightType = "custom", custWeights = lakes_nseqs)
    sitepairs <- na.omit(sitepairs)
    
    # Compute GDM with significant predictors
    gdm_signif <- gdm(sitepairs, geo = geodist)
    
  } else if (length(model_signif) == 0) {
    print("No significant model")
    
    gdm_signif <- NULL
  }
  
  return(gdm_signif)
}

gdm_signif_protists_bray_geography <- computeGDM(gdm_geography, asv_bray, lakes_nseqs_all, geography)
gdm_signif_protists_bray_morphometry <- computeGDM(gdm_morphometry, asv_bray, lakes_nseqs_all, morphometry)
gdm_signif_protists_bray_landuse_soil <- computeGDM(gdm_landuse_soil, asv_bray, lakes_nseqs_all, landuse_soil)
gdm_signif_protists_bray_climate <- computeGDM(gdm_climate, asv_bray, lakes_nseqs_all, climate)
gdm_signif_protists_bray_waterchem <- computeGDM(gdm_waterchem, asv_bray, lakes_nseqs_all, waterchem)

gdm_signif_protists_genunifrac_geography <- computeGDM(phylogdm_geography, unifrac_gen0.5alpha_dist, lakes_nseqs_all, geography)
gdm_signif_protists_genunifrac_morphometry <- computeGDM(phylogdm_morphometry, unifrac_gen0.5alpha_dist, lakes_nseqs_all, morphometry)
gdm_signif_protists_genunifrac_landuse_soil <- computeGDM(phylogdm_landuse_soil, unifrac_gen0.5alpha_dist, lakes_nseqs_all, landuse_soil)
gdm_signif_protists_genunifrac_climate <- computeGDM(phylogdm_climate, unifrac_gen0.5alpha_dist, lakes_nseqs_all, climate)
gdm_signif_protists_genunifrac_waterchem <- computeGDM(phylogdm_waterchem, unifrac_gen0.5alpha_dist, lakes_nseqs_all, waterchem)

gdm_signif_phototrophs_bray_geography <- computeGDM(gdm_phototrophs_geography, phototrophs_bray, lakes_nseqs_phototrophs, geography)
gdm_signif_phototrophs_bray_morphometry <- computeGDM(gdm_phototrophs_morphometry, phototrophs_bray, lakes_nseqs_phototrophs, morphometry)
gdm_signif_phototrophs_bray_landuse_soil <- computeGDM(gdm_phototrophs_landuse_soil, phototrophs_bray, lakes_nseqs_phototrophs, landuse_soil)
gdm_signif_phototrophs_bray_climate <- computeGDM(gdm_phototrophs_climate, phototrophs_bray, lakes_nseqs_phototrophs, climate)
gdm_signif_phototrophs_bray_waterchem <- computeGDM(gdm_phototrophs_waterchem, phototrophs_bray, lakes_nseqs_phototrophs, waterchem)

gdm_signif_consumers_bray_geography <- computeGDM(gdm_consumers_geography, consumers_bray, lakes_nseqs_consumers, geography)
gdm_signif_consumers_bray_morphometry <- computeGDM(gdm_consumers_morphometry, consumers_bray, lakes_nseqs_consumers, morphometry)
gdm_signif_consumers_bray_landuse_soil <- computeGDM(gdm_consumers_landuse_soil, consumers_bray, lakes_nseqs_consumers, landuse_soil)
gdm_signif_consumers_bray_climate <- computeGDM(gdm_consumers_climate, consumers_bray, lakes_nseqs_consumers, climate)
gdm_signif_consumers_bray_waterchem <- computeGDM(gdm_consumers_waterchem, consumers_bray, lakes_nseqs_consumers, waterchem)

gdm_signif_mixotrophs_bray_geography <- computeGDM(gdm_mixotrophs_geography, mixotrophs_bray, lakes_nseqs_mixotrophs, geography)
gdm_signif_mixotrophs_bray_morphometry <- computeGDM(gdm_mixotrophs_morphometry, mixotrophs_bray, lakes_nseqs_mixotrophs, morphometry)
gdm_signif_mixotrophs_bray_landuse_soil <- computeGDM(gdm_mixotrophs_landuse_soil, mixotrophs_bray, lakes_nseqs_mixotrophs, landuse_soil)
gdm_signif_mixotrophs_bray_climate <- computeGDM(gdm_mixotrophs_climate, mixotrophs_bray, lakes_nseqs_mixotrophs, climate)
gdm_signif_mixotrophs_bray_waterchem <- computeGDM(gdm_mixotrophs_waterchem, mixotrophs_bray, lakes_nseqs_mixotrophs, waterchem)


#### Summarize GDM permutation test results ####
summarizeGDM <- function(gdm, dist, lakes_nseqs, envtable, response_vars, explanatory_vars) {
  gdm_signif <- gdm[[3]] %>%
    as_tibble(rownames = "predictor") %>%
    clean_names()
  
  model_signif <- gdm_signif %>%
    pivot_longer(!predictor, names_to = "model", values_to = "p_value") %>%
    mutate(model_number = case_when(grepl("full_model_", model) ~ as.numeric(str_remove(model, "full_model_")),
                                    model == "full_model" ~ 0)) %>%
    group_by(model, model_number) %>%
    summarize(max_p_value = max(p_value, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(max_p_value < 0.05) %>%
    slice_min(model_number) %>%
    pull(model)
  
  if (length(model_signif) == 1) {
    predictors_signif <- gdm_signif %>%
      dplyr::select(predictor, p_value = model_signif) %>%
      filter(!is.na(p_value)) %>%
      pull(predictor)
    
    if ("Geographic" %in% predictors_signif) {
      geodist <- TRUE
    } else if (!"Geographic" %in% predictors_signif) {
      geodist <- FALSE
    }
    
    predictors_signif <- predictors_signif[!predictors_signif %in% "Geographic"]
    
    # Format dissimilarity matrix
    dissim <- dist %>%
      as.matrix() %>%
      as_tibble(rownames = "lake_id") %>%
      arrange(lake_id)
    
    # Filter study lakes in explanatory data
    if ("latitude" %in% names(envtable) | "longitude" %in% names(envtable)) {
      envtable <- envtable %>%
        filter(lake_id %in% lakes) %>%
        dplyr::select(lake_id, where(is.numeric)) %>%
        dplyr::select(lake_id, latitude, longitude, all_of(predictors_signif)) %>%
        arrange(lake_id)
    } else {
      envtable <- latlong %>%
        filter(lake_id %in% lakes) %>%
        left_join(envtable, by = "lake_id") %>%
        dplyr::select(lake_id, where(is.numeric)) %>%
        dplyr::select(lake_id, latitude, longitude, all_of(predictors_signif)) %>%
        arrange(lake_id)
    }
    
    # Format site-pairs
    sitepairs <- formatsitepair(bioData = dissim, bioFormat = 3,
                                abundance = TRUE, siteColumn = "lake_id",
                                XColumn = "longitude", YColumn = "latitude",
                                predData = envtable,
                                weightType = "custom", custWeights = lakes_nseqs)
    sitepairs <- na.omit(sitepairs)
    
    # Compute GDM with significant predictors
    gdm_signif <- gdm(sitepairs, geo = geodist)
    
    
    # Extract GDM splines
    gdm_splines <- bind_cols(isplineExtract(gdm_signif)$x %>%
                               as_tibble() %>%
                               pivot_longer(everything(), names_to = "predictor", values_to = "x"),
                             isplineExtract(gdm_signif)$y %>%
                               as_tibble() %>%
                               pivot_longer(everything(), names_to = "predictor", values_to = "turnover") %>%
                               dplyr::select(turnover))
    
    # Summarize GDM coefficients
    gdm_coefficients <- tibble(predictor = gdm_signif$predictors,
                               coefficient1 = gdm_signif$coefficients[seq(1, length(gdm_signif$coefficients), 3)],
                               coefficient2 = gdm_signif$coefficients[seq(2, length(gdm_signif$coefficients), 3)],
                               coefficient3 = gdm_signif$coefficients[seq(3, length(gdm_signif$coefficients), 3)]) %>%
      group_by(predictor) %>%
      summarize(sum_coefficients = sum(coefficient1, coefficient2, coefficient3)) %>%
      ungroup() %>%
      arrange(-sum_coefficients) %>%
      mutate(predictor_coefficient = paste0(predictor, " (", round(sum_coefficients, 2), ")"))
    
    deviance_explained <- round(gdm_signif$explained, 1)
    
    gdm_summary <- tibble(response_vars = response_vars,
                          explanatory_vars = explanatory_vars,
                          model_number = model_signif,
                          deviance_explained_pct = deviance_explained,
                          predictors = toString(gdm_coefficients$predictor_coefficient))
    
  } else if (length(model_signif) == 0) {
    print("No significant model")
    
    gdm_summary <- tibble(response_vars = response_vars,
                          explanatory_vars = explanatory_vars,
                          model_number = "NS",
                          deviance_explained_pct = NA,
                          predictors = "NS")
  }
  
  return(gdm_summary)
}

gdm_summary <- summarizeGDM(gdm_geography, asv_bray, lakes_nseqs_all, geography, "protists|bray", "geography") %>%
  bind_rows(summarizeGDM(gdm_morphometry, asv_bray, lakes_nseqs_all, morphometry, "protists|bray", "morphometry")) %>%
  bind_rows(summarizeGDM(gdm_landuse_soil, asv_bray, lakes_nseqs_all, landuse_soil, "protists|bray", "landuse_soil")) %>%
  bind_rows(summarizeGDM(gdm_climate, asv_bray, lakes_nseqs_all, climate, "protists|bray", "climate")) %>%
  bind_rows(summarizeGDM(gdm_waterchem, asv_bray, lakes_nseqs_all, waterchem, "protists|bray", "waterchem")) %>%
  
  bind_rows(summarizeGDM(phylogdm_geography, unifrac_gen0.5alpha_dist, lakes_nseqs_all, geography, "protists|genunifrac", "geography")) %>%
  bind_rows(summarizeGDM(phylogdm_morphometry, unifrac_gen0.5alpha_dist, lakes_nseqs_all, morphometry, "protists|genunifrac", "morphometry")) %>%
  bind_rows(summarizeGDM(phylogdm_landuse_soil, unifrac_gen0.5alpha_dist, lakes_nseqs_all, landuse_soil, "protists|genunifrac", "landuse_soil")) %>%
  bind_rows(summarizeGDM(phylogdm_climate, unifrac_gen0.5alpha_dist, lakes_nseqs_all, climate, "protists|genunifrac", "climate")) %>%
  bind_rows(summarizeGDM(phylogdm_waterchem, unifrac_gen0.5alpha_dist, lakes_nseqs_all, waterchem, "protists|genunifrac", "waterchem")) %>%
  
  bind_rows(summarizeGDM(gdm_phototrophs_geography, phototrophs_bray, lakes_nseqs_phototrophs, geography, "phototrophs|bray", "geography")) %>%
  bind_rows(summarizeGDM(gdm_phototrophs_morphometry, phototrophs_bray, lakes_nseqs_phototrophs, morphometry, "phototrophs|bray", "morphometry")) %>%
  bind_rows(summarizeGDM(gdm_phototrophs_landuse_soil, phototrophs_bray, lakes_nseqs_phototrophs, landuse_soil, "phototrophs|bray", "landuse_soil")) %>%
  bind_rows(summarizeGDM(gdm_phototrophs_climate, phototrophs_bray, lakes_nseqs_phototrophs, climate, "phototrophs|bray", "climate")) %>%
  bind_rows(summarizeGDM(gdm_phototrophs_waterchem, phototrophs_bray, lakes_nseqs_phototrophs, waterchem, "phototrophs|bray", "waterchem")) %>%
  
  bind_rows(summarizeGDM(gdm_consumers_geography, consumers_bray, lakes_nseqs_consumers, geography, "consumers|bray", "geography")) %>%
  bind_rows(summarizeGDM(gdm_consumers_morphometry, consumers_bray, lakes_nseqs_consumers, morphometry, "consumers|bray", "morphometry")) %>%
  bind_rows(summarizeGDM(gdm_consumers_landuse_soil, consumers_bray, lakes_nseqs_consumers, landuse_soil, "consumers|bray", "landuse_soil")) %>%
  bind_rows(summarizeGDM(gdm_consumers_climate, consumers_bray, lakes_nseqs_consumers, climate, "consumers|bray", "climate")) %>%
  bind_rows(summarizeGDM(gdm_consumers_waterchem, consumers_bray, lakes_nseqs_consumers, waterchem, "consumers|bray", "waterchem")) %>%
  
  bind_rows(summarizeGDM(gdm_mixotrophs_geography, mixotrophs_bray, lakes_nseqs_mixotrophs, geography, "mixotrophs|bray", "geography")) %>%
  bind_rows(summarizeGDM(gdm_mixotrophs_morphometry, mixotrophs_bray, lakes_nseqs_mixotrophs, morphometry, "mixotrophs|bray", "morphometry")) %>%
  bind_rows(summarizeGDM(gdm_mixotrophs_landuse_soil, mixotrophs_bray, lakes_nseqs_mixotrophs, landuse_soil, "mixotrophs|bray", "landuse_soil")) %>%
  bind_rows(summarizeGDM(gdm_mixotrophs_climate, mixotrophs_bray, lakes_nseqs_mixotrophs, climate, "mixotrophs|bray", "climate")) %>%
  bind_rows(summarizeGDM(gdm_mixotrophs_waterchem, mixotrophs_bray, lakes_nseqs_mixotrophs, waterchem, "mixotrophs|bray", "waterchem"))

# Write GDM summary table to file
# gdm_summary %>%
#   write_csv("output/gdm/lp2017-2019watercolumn18s_gdm_summary.csv", col_names = TRUE)

# Summarize deviance explained by each model
gdm_deviance_summary <- gdm_summary %>%
  dplyr::select(response_vars, explanatory_vars, deviance_explained_pct) %>%
  pivot_wider(names_from = explanatory_vars, values_from = deviance_explained_pct) %>%
  as.data.frame()
gdm_deviance_summary[is.na(gdm_deviance_summary)] <- "NS"

# Write deviance explained summary to file
# gdm_deviance_summary %>%
#   write_csv("output/gdm/lp2017-2019watercolumn18s_gdm_deviance_summary.csv", col_names = TRUE)

# Summarize deviance explained by each model (transposed summary matrix)
gdm_deviance_summary_t <- gdm_summary %>%
  mutate(deviance_explained_pct = round(deviance_explained_pct, 0)) %>%
  dplyr::select(response_vars, explanatory_vars, deviance_explained_pct) %>%
  pivot_wider(names_from = response_vars, values_from = deviance_explained_pct) %>%
  arrange(factor(explanatory_vars, levels = c("waterchem", "landuse_soil", "morphometry", "climate", "geography"))) %>%
  as.data.frame()
gdm_deviance_summary_t[is.na(gdm_deviance_summary_t)] <- "NS"

# Write deviance explained (transposed) summary to file
# gdm_deviance_summary_t %>%
#   write_csv("output/gdm/lp2017-2019watercolumn18s_gdm_deviance_summary_t.csv", col_names = TRUE)


#### Plot predictor partial effects ####
trophic <- metadata %>%
  filter(lake_id %in% lakes) %>%
  dplyr::select(lake_id, trophic_status)

# Define function to plot (phylo)GDM splines representing both taxonomic and phylogenetic turnover
plotGDMs <- function(gdm, phylogdm, explanatory_vars, plot_nrows = 1) {
  gdm_splines <- bind_cols(isplineExtract(gdm)$x %>%
                             as_tibble() %>%
                             pivot_longer(everything(), names_to = "predictor", values_to = "x"),
                           isplineExtract(gdm)$y %>%
                             as_tibble() %>%
                             pivot_longer(everything(), names_to = "predictor", values_to = "taxonomic") %>%
                             dplyr::select(taxonomic))
  
  phylogdm_splines <- bind_cols(isplineExtract(phylogdm)$x %>%
                                  as_tibble() %>%
                                  pivot_longer(everything(), names_to = "predictor", values_to = "x"),
                                isplineExtract(phylogdm)$y %>%
                                  as_tibble() %>%
                                  pivot_longer(everything(), names_to = "predictor", values_to = "phylogenetic") %>%
                                  dplyr::select(phylogenetic))
  
  all_splines <- gdm_splines %>%
    full_join(phylogdm_splines, by = c("predictor", "x")) %>%
    mutate(predictor = case_when(predictor == "latitude_y" ~ "latitude",
                                 predictor == "longitude_x" ~ "longitude",
                                 TRUE ~ predictor))
  
  predictors_all <- unique(all_splines$predictor)
  predictors_all <- predictors_all[!predictors_all %in% "Geographic"]
  
  gdm_splines_max <- max(all_splines$taxonomic, na.rm = TRUE)
  phylogdm_splines_max <- max(all_splines$phylogenetic, na.rm = TRUE)
  
  #coeff <- round(phylogdm_splines_max/gdm_splines_max, 1) + 0.1
  coeff <- phylogdm_splines_max/gdm_splines_max
  
  all_splines %>%
    ggplot(aes(x = x)) +
    geom_line(aes(y = phylogenetic / coeff), size = 1, colour = "#A0A0A0") +
    geom_line(aes(y = taxonomic), size = 1, colour = "#7CACD3") +
    geom_point(aes(x = value, y = -0.05, colour = trophic_status),
               data = metadata %>%
                 dplyr::select(lake_id, all_of(predictors_all)) %>%
                 pivot_longer(!lake_id, names_to = "predictor", values_to = "value") %>%
                 left_join(trophic, by = "lake_id"),
               alpha = 0.7, shape = 39) +
    facet_wrap(~predictor, scales = "free_x", nrow = plot_nrows, labeller = labeller(predictor = variable_labels)) +
    scale_colour_manual(values = palette_trophic) +
    scale_y_continuous(name = "Taxonomic turnover",
                       sec.axis = sec_axis(~.*coeff, name = "Phylogenetic turnover")) +
    #ggtitle(str_to_sentence(explanatory_vars)) +
    theme_bw() %+replace%
    theme(legend.position = "none",
          strip.background = element_rect(fill = "black"),
          strip.text = element_text(colour = "white", face = "bold"),
          strip.text.x = element_text(face = "bold"),
          axis.title.y = element_text(colour = "#7CACD3", angle = 90),
          axis.title.y.right = element_text(colour = "#A0A0A0", angle = 270),
          axis.title.x = element_blank(),
          axis.text.y = element_text(colour = "#7CACD3"),
          axis.ticks.y = element_line(colour = "#7CACD3"),
          axis.text.y.right = element_text(colour = "#A0A0A0"),
          axis.ticks.y.right = element_line(colour = "#A0A0A0"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
}

(gdm_plots_protists_geography <- plotGDMs(gdm_signif_protists_bray_geography, gdm_signif_protists_genunifrac_geography, "Geography", plot_nrows = 1))
#ggsave("figures/lp2017-2019watercolumn18s_gdms_protists_geography.pdf", gdm_plots_protists_geography, "pdf", width = 6, height = 3, units = "in")

(gdm_plots_protists_morphometry <- plotGDMs(gdm_signif_protists_bray_morphometry, gdm_signif_protists_genunifrac_morphometry, "Morphometry", plot_nrows = 1))
#ggsave("figures/lp2017-2019watercolumn18s_gdms_protists_morphometry.pdf", gdm_plots_protists_morphometry, "pdf", width = 6, height = 3, units = "in")

(gdm_plots_protists_landuse_soil <- plotGDMs(gdm_signif_protists_bray_landuse_soil, gdm_signif_protists_genunifrac_landuse_soil, "Land use/soil", plot_nrows = 1))
#ggsave("figures/lp2017-2019watercolumn18s_gdms_protists_landuse_soil.pdf", gdm_plots_protists_landuse_soil, "pdf", width = 8, height = 3, units = "in")

(gdm_plots_protists_climate <- plotGDMs(gdm_signif_protists_bray_climate, gdm_signif_protists_genunifrac_climate, "Weather", plot_nrows = 1))
#ggsave("figures/lp2017-2019watercolumn18s_gdms_protists_climate.pdf", gdm_plots_protists_climate, "pdf", width = 4, height = 3, units = "in")

(gdm_plots_protists_waterchem <- plotGDMs(gdm_signif_protists_bray_waterchem, gdm_signif_protists_genunifrac_waterchem, "Physicochemistry", plot_nrows = 3))
#(gdm_plots_protists_waterchem <- plotGDMs(gdm_signif_protists_bray_waterchem, gdm_signif_protists_genunifrac_waterchem, "Physicochemistry", plot_nrows = 1))
#ggsave("figures/lp2017-2019watercolumn18s_gdms_protists_waterchem.pdf", gdm_plots_protists_waterchem, "pdf", width = 8, height = 9, units = "in")

# Combine GDM spline plots
gdm_protists_plots_all_layout <- "
AAAAAAAAAAA
BBBB#######
CCC########
"
(gdm_protists_plots_all <- gdm_plots_protists_waterchem /
    gdm_plots_protists_landuse_soil /
    gdm_plots_protists_morphometry + 
    plot_layout(design = gdm_protists_plots_all_layout, guides = "collect"))
#ggsave("figures/lp2017-2019watercolumn18s_gdms_protists_all.pdf", gdm_protists_plots_all, "pdf", width = 22, height = 7.5, units = "in")

# Define function to plot GDM splines representing phototroph, heterotroph, and mixotroph turnover
plotGDMsPhotoConsMixo <- function(gdm_phototrophs, gdm_consumers, gdm_mixotrophs, plot_nrows = NULL) {
  gdm_splines_phototrophs <- bind_cols(isplineExtract(gdm_phototrophs)$x %>%
                                         as_tibble() %>%
                                         pivot_longer(everything(), names_to = "predictor", values_to = "x"),
                                       isplineExtract(gdm_phototrophs)$y %>%
                                         as_tibble() %>%
                                         pivot_longer(everything(), names_to = "predictor", values_to = "phototrophs") %>%
                                         dplyr::select(phototrophs))
  
  gdm_splines_consumers <- bind_cols(isplineExtract(gdm_consumers)$x %>%
                                       as_tibble() %>%
                                       pivot_longer(everything(), names_to = "predictor", values_to = "x"),
                                     isplineExtract(gdm_consumers)$y %>%
                                       as_tibble() %>%
                                       pivot_longer(everything(), names_to = "predictor", values_to = "consumers") %>%
                                       dplyr::select(consumers))
  
  gdm_splines_mixotrophs <- bind_cols(isplineExtract(gdm_mixotrophs)$x %>%
                                        as_tibble() %>%
                                        pivot_longer(everything(), names_to = "predictor", values_to = "x"),
                                      isplineExtract(gdm_mixotrophs)$y %>%
                                        as_tibble() %>%
                                        pivot_longer(everything(), names_to = "predictor", values_to = "mixotrophs") %>%
                                        dplyr::select(mixotrophs))
  
  all_splines <- gdm_splines_phototrophs %>%
    full_join(gdm_splines_consumers, by = c("predictor", "x")) %>%
    full_join(gdm_splines_mixotrophs, by = c("predictor", "x"))
  
  predictors_all <- unique(all_splines$predictor)
  predictors_all <- predictors_all[!predictors_all %in% "Geographic"]
  
  #gdm_splines_max_phototrophs <- max(all_splines$phototrophs, na.rm = TRUE)
  #gdm_splines_max_consumers <- max(all_splines$consumers, na.rm = TRUE)
  #gdm_splines_max_mixotrophs <- max(all_splines$mixotrophs, na.rm = TRUE)
  
  #coeff <- round(gdm_splines_max_consumers/gdm_splines_max_phototrophs, 1) + 0.1
  
  all_splines %>%
    ggplot(aes(x = x)) +
    geom_line(aes(y = phototrophs), size = 1, colour = "#FDFF80") +
    geom_line(aes(y = consumers), size = 1, colour = "#76BED0") +
    geom_line(aes(y = mixotrophs), size = 1, colour = "#62BE81") +
    geom_point(aes(x = value, y = -0.05, colour = trophic_status),
               data = metadata %>%
                 dplyr::select(lake_id, all_of(predictors_all)) %>%
                 pivot_longer(!lake_id, names_to = "predictor", values_to = "value") %>%
                 left_join(trophic, by = "lake_id"),
               alpha = 0.7, shape = 39) +
    facet_wrap(~predictor, scales = "free_x", nrow = plot_nrows,
               labeller = labeller(predictor = variable_labels)) +
    scale_y_continuous(name = "Turnover") +
    scale_colour_manual(values = palette_trophic) +
    theme_light() %+replace%
    theme(legend.position = "none", strip.background = element_rect(fill = "black"),
          strip.text.x = element_text(face = "bold"),
          axis.title.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
}

# Define function to plot GDM splines representing phototroph and heterotroph turnover
plotGDMsPhotoCons <- function(gdm_phototrophs, gdm_consumers, plot_nrows = NULL) {
  gdm_splines_phototrophs <- bind_cols(isplineExtract(gdm_phototrophs)$x %>%
                                         as_tibble() %>%
                                         pivot_longer(everything(), names_to = "predictor", values_to = "x"),
                                       isplineExtract(gdm_phototrophs)$y %>%
                                         as_tibble() %>%
                                         pivot_longer(everything(), names_to = "predictor", values_to = "phototrophs") %>%
                                         dplyr::select(phototrophs))
  
  gdm_splines_consumers <- bind_cols(isplineExtract(gdm_consumers)$x %>%
                                       as_tibble() %>%
                                       pivot_longer(everything(), names_to = "predictor", values_to = "x"),
                                     isplineExtract(gdm_consumers)$y %>%
                                       as_tibble() %>%
                                       pivot_longer(everything(), names_to = "predictor", values_to = "consumers") %>%
                                       dplyr::select(consumers))
  
  all_splines <- gdm_splines_phototrophs %>%
    full_join(gdm_splines_consumers, by = c("predictor", "x"))
  
  predictors_all <- unique(all_splines$predictor)
  predictors_all <- predictors_all[!predictors_all %in% "Geographic"]
  
  phototrophs_splines_max <- max(all_splines$phototrophs, na.rm = TRUE)
  consumers_splines_max <- max(all_splines$consumers, na.rm = TRUE)
  
  coeff <- consumers_splines_max/phototrophs_splines_max
  
  all_splines %>%
    ggplot(aes(x = x)) +
    geom_line(aes(y = phototrophs), size = 1, colour = "#FDFF80") +
    geom_line(aes(y = consumers / coeff), size = 1, colour = "#76BED0") +
    geom_point(aes(x = value, y = -0.05, colour = trophic_status),
               data = metadata %>%
                 dplyr::select(lake_id, all_of(predictors_all)) %>%
                 pivot_longer(!lake_id, names_to = "predictor", values_to = "value") %>%
                 left_join(trophic, by = "lake_id"),
               alpha = 0.7, shape = 39) +
    facet_wrap(~predictor, scales = "free_x", nrow = plot_nrows,
               labeller = labeller(predictor = variable_labels)) +
    scale_y_continuous(name = "Turnover") +
    scale_colour_manual(values = palette_trophic) +
    scale_y_continuous(name = "Phototroph turnover (yellow)",
                       sec.axis = sec_axis(~.*coeff, name = "Heterotroph turnover")) +
    #ggtitle(str_to_sentence(explanatory_vars)) +
    theme_bw() %+replace%
    theme(legend.position = "none",
          strip.background = element_rect(fill = "black"),
          strip.text = element_text(colour = "white", face = "bold"),
          strip.text.x = element_text(face = "bold"),
          axis.title.y = element_text(angle = 90),
          axis.title.y.right = element_text(colour = "#76BED0", angle = 270),
          axis.title.x = element_blank(),
          #axis.text.y = element_text(colour = "#FDFF80"),
          #axis.ticks.y = element_line(colour = "#FDFF80"),
          axis.text.y.right = element_text(colour = "#76BED0"),
          axis.ticks.y.right = element_line(colour = "#76BED0"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
}

# Define function to plot GDM splines representing phototroph turnover only
plotGDMsPhoto <- function(gdm_phototrophs, plot_nrows = NULL) {
  gdm_splines_phototrophs <- bind_cols(isplineExtract(gdm_phototrophs)$x %>%
                                         as_tibble() %>%
                                         pivot_longer(everything(), names_to = "predictor", values_to = "x"),
                                       isplineExtract(gdm_phototrophs)$y %>%
                                         as_tibble() %>%
                                         pivot_longer(everything(), names_to = "predictor", values_to = "phototrophs") %>%
                                         dplyr::select(phototrophs))
  
  all_splines <- gdm_splines_phototrophs
  
  predictors_all <- unique(all_splines$predictor)
  predictors_all <- predictors_all[!predictors_all %in% "Geographic"]
  
  all_splines %>%
    ggplot(aes(x = x)) +
    geom_line(aes(y = phototrophs), size = 1, colour = "#FDFF80") +
    geom_point(aes(x = value, y = -0.05, colour = trophic_status),
               data = metadata %>%
                 dplyr::select(lake_id, all_of(predictors_all)) %>%
                 pivot_longer(!lake_id, names_to = "predictor", values_to = "value") %>%
                 left_join(trophic, by = "lake_id"),
               alpha = 0.7, shape = 39) +
    facet_wrap(~predictor, scales = "free_x", nrow = plot_nrows,
               labeller = labeller(predictor = variable_labels)) +
    scale_y_continuous(name = "Phototroph turnover") +
    scale_colour_manual(values = palette_trophic) +
    theme_light() %+replace%
    theme(legend.position = "none", strip.background = element_rect(fill = "black"),
          strip.text.x = element_text(face = "bold"),
          axis.title.y = element_text(angle = 90),
          axis.title.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
}

# Define function to plot GDM splines representing heterotroph turnover only
plotGDMsCons <- function(gdm_consumers, plot_nrows = NULL) {
  gdm_splines_consumers <- bind_cols(isplineExtract(gdm_consumers)$x %>%
                                       as_tibble() %>%
                                       pivot_longer(everything(), names_to = "predictor", values_to = "x"),
                                     isplineExtract(gdm_consumers)$y %>%
                                       as_tibble() %>%
                                       pivot_longer(everything(), names_to = "predictor", values_to = "consumers") %>%
                                       dplyr::select(consumers))
  
  all_splines <- gdm_splines_consumers
  
  predictors_all <- unique(all_splines$predictor)
  predictors_all <- predictors_all[!predictors_all %in% "Geographic"]
  
  all_splines %>%
    ggplot(aes(x = x)) +
    geom_line(aes(y = consumers), size = 1, colour = "#76BED0") +
    geom_point(aes(x = value, y = -0.05, colour = trophic_status),
               data = metadata %>%
                 dplyr::select(lake_id, all_of(predictors_all)) %>%
                 pivot_longer(!lake_id, names_to = "predictor", values_to = "value") %>%
                 left_join(trophic, by = "lake_id"),
               alpha = 0.7, shape = 39) +
    facet_wrap(~predictor, scales = "free_x", nrow = plot_nrows,
               labeller = labeller(predictor = variable_labels)) +
    scale_y_continuous(name = "Heterotroph turnover") +
    scale_colour_manual(values = palette_trophic) +
    theme_light() %+replace%
    theme(legend.position = "none", strip.background = element_rect(fill = "black"),
          strip.text.x = element_text(face = "bold"),
          axis.title.y = element_text(colour = "#76BED0", angle = 90),
          axis.text.y = element_text(colour = "#76BED0"),
          axis.title.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
}

(gdm_functions_plots_geography <- plotGDMsCons(gdm_signif_consumers_bray_geography, plot_nrows = 1))
#ggsave("figures/lp2017-2019watercolumn18s_gdms_geography_functions.pdf", gdm_functions_plots_geography, "pdf", width = 4, height = 3, units = "in")

(gdm_functions_plots_morphometry <- plotGDMsPhoto(gdm_signif_phototrophs_bray_morphometry, plot_nrows = 1))
#ggsave("figures/lp2017-2019watercolumn18s_gdms_morphometry_functions.pdf", gdm_functions_plots_morphometry, "pdf", width = 6, height = 3, units = "in")

(gdm_functions_plots_landuse_soil <- plotGDMsPhotoCons(gdm_signif_phototrophs_bray_landuse_soil, gdm_signif_consumers_bray_landuse_soil, plot_nrows = 2))
#(gdm_functions_plots_landuse_soil <- plotGDMsPhotoCons(gdm_signif_phototrophs_bray_landuse_soil, gdm_signif_consumers_bray_landuse_soil, plot_nrows = 1))
#ggsave("figures/lp2017-2019watercolumn18s_gdms_landuse_soil_functions.pdf", gdm_functions_plots_landuse_soil, "pdf", width = 6, height = 6, units = "in")

(gdm_functions_plots_climate <- plotGDMsCons(gdm_signif_consumers_bray_climate, plot_nrows = 1))
#ggsave("figures/lp2017-2019watercolumn18s_gdms_climate_functions.pdf", gdm_functions_plots_climate, "pdf", width = 4, height = 3, units = "in")

(gdm_functions_plots_waterchem <- plotGDMsPhotoConsMixo(gdm_signif_phototrophs_bray_waterchem, gdm_signif_consumers_bray_waterchem, gdm_signif_mixotrophs_bray_waterchem, plot_nrows = 3))
#(gdm_functions_plots_waterchem <- plotGDMsPhotoConsMixo(gdm_signif_phototrophs_bray_waterchem, gdm_signif_consumers_bray_waterchem, gdm_signif_mixotrophs_bray_waterchem, plot_nrows = 1))
#ggsave("figures/lp2017-2019watercolumn18s_gdms_waterchem_functions.pdf", gdm_functions_plots_waterchem, "pdf", width = 8, height = 9, units = "in")

# Combine GDM spline plots
gdm_functions_plots_all_layout <- "
AAAAAAAAAAA
BBBBB######
CCC########
DD#########
EE#########
"
(gdm_functions_plots_all <- gdm_functions_plots_waterchem /
    gdm_functions_plots_landuse_soil /
    gdm_functions_plots_morphometry /
    gdm_functions_plots_climate /
    gdm_functions_plots_geography + 
    plot_layout(design = gdm_functions_plots_all_layout, guides = "collect"))
#ggsave("figures/lp2017-2019watercolumn18s_gdms_functions_all.pdf", gdm_functions_plots_all, "pdf", width = 22, height = 12.5, units = "in")
