# Create site-pairs for generalized dissimilarity modelling (GDM)
# https://cran.r-project.org/web/packages/gdm/vignettes/gdmVignette.pdf

# Load libraries
library(tidyverse)
library(vegan)
library(gdm)
library(janitor)


#### Import and format data ####
# Import melted protist ASV table
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes_function.tsv")
length(unique(asv_melt$asv_code))

# Import metadata
source("scripts/05_metadata.R")


#### Compute dissimilarity matrix ####
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

# Create ASV table
asvtable_phototrophs <- phototrophs_melt %>%
  dplyr::select(lake_id, asv_code, relseqs) %>%
  pivot_wider(names_from = asv_code, values_from = relseqs, values_fill = 0) %>%
  arrange(lake_id) %>%
  column_to_rownames("lake_id")

# Compute dissimilarity matrix
phototrophs_bray <- vegdist(asvtable_phototrophs, "bray")

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

# Create ASV table
asvtable_consumers <- consumers_melt %>%
  dplyr::select(lake_id, asv_code, relseqs) %>%
  pivot_wider(names_from = asv_code, values_from = relseqs, values_fill = 0) %>%
  arrange(lake_id) %>%
  column_to_rownames("lake_id")

# Compute dissimilarity matrix
consumers_bray <- vegdist(asvtable_consumers, "bray")

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

# Create ASV table
asvtable_mixotrophs <- mixotrophs_melt %>%
  dplyr::select(lake_id, asv_code, relseqs) %>%
  pivot_wider(names_from = asv_code, values_from = relseqs, values_fill = 0) %>%
  arrange(lake_id) %>%
  column_to_rownames("lake_id")

# Compute dissimilarity matrix
mixotrophs_bray <- vegdist(asvtable_mixotrophs, "bray")


#### Create site-pair table ####
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

lakes_nseqs_phototrophs <- countSeqs(phototrophs_melt)
lakes_nseqs_consumers <- countSeqs(consumers_melt)
lakes_nseqs_mixotrophs <- countSeqs(mixotrophs_melt)

# Combine response (community composition) data and predictor (environmental) data into site-pair table format
# Define function to create GDM site-pair table
createSitePair <- function(dist, lakes_nseqs, envtable, geodist = FALSE) {
  # Format dissimilarity matrix
  dissim <- dist %>%
    as.matrix() %>%
    as_tibble(rownames = "lake_id") %>%
    arrange(lake_id)
  
  # Filter study lakes in explanatory data
  envtable <- latlong %>%
    filter(lake_id %in% lakes) %>%
    left_join(envtable, by = "lake_id") %>%
    dplyr::select(lake_id, where(is.numeric)) %>%
    arrange(lake_id)
  
  # Format site-pairs
  sitepairs <- formatsitepair(bioData = dissim, bioFormat = 3,
                              abundance = TRUE, siteColumn = "lake_id",
                              XColumn = "longitude", YColumn = "latitude",
                              predData = envtable,
                              weightType = "custom", custWeights = lakes_nseqs)
  sitepairs <- na.omit(sitepairs)
  print(paste0("Number of rows in site-pair table: ", nrow(sitepairs)))
  
  # Compute GDM with original set of predictors
  gdm <- gdm(sitepairs, geo = geodist)
  print(paste0("Variance explained by GDM built with original set of predictors: ",
               round(gdm$explained, 2), "%"))
  
  # Summarize GDM coefficients
  gdm_coefficients <- tibble(predictor = gdm$predictors,
                             coefficient1 = gdm$coefficients[seq(1, length(gdm$coefficients), 3)],
                             coefficient2 = gdm$coefficients[seq(2, length(gdm$coefficients), 3)],
                             coefficient3 = gdm$coefficients[seq(3, length(gdm$coefficients), 3)])
  
  # Retain predictors with non-zero coefficients
  gdm_predictors <- gdm_coefficients %>%
    mutate(sum_coefficients = coefficient1 + coefficient2 + coefficient3) %>%
    filter(sum_coefficients > 0) %>%
    filter(predictor != "Geographic") %>%
    pull(predictor)
  
  print(paste0("Retained predictor: ", toString(gdm_predictors)))
  
  predictors <- envtable %>%
    dplyr::select(lake_id, latitude, longitude, all_of(gdm_predictors))
  
  # Create new site-pairs with retained predictors
  sitepairs_new <- formatsitepair(bioData = dissim, bioFormat = 3,
                                  abundance = TRUE, siteColumn = "lake_id",
                                  XColumn = "longitude", YColumn = "latitude",
                                  predData = predictors,
                                  weightType = "custom", custWeights = lakes_nseqs)
  
  sitepairs_new <- na.omit(sitepairs_new)
  
  # Recompute GDM with new site-pairs
  gdm_new <- gdm(sitepairs_new, geo = geodist)
  print(paste0("Variance explained by GDM built with updated set of predictors: ",
               round(gdm_new$explained, 2), "%"))
  
  if (length(gdm_predictors) >= 3) {
    return(sitepairs_new)
  } else {
    print("Returning original site-pair table without removing variables")
    return(sitepairs)
  }
}

# Create site-pair tables
phototrophs_sitepair_geography <- createSitePair(phototrophs_bray, lakes_nseqs_phototrophs,
                                                 geography %>%
                                                   dplyr::select(-latitude, -longitude), geodist = TRUE)
phototrophs_sitepair_morphometry <- createSitePair(phototrophs_bray, lakes_nseqs_phototrophs,
                                                   morphometry, geodist = FALSE)
phototrophs_sitepair_landuse_soil <- createSitePair(phototrophs_bray, lakes_nseqs_phototrophs,
                                                    landuse_soil, geodist = FALSE)
phototrophs_sitepair_climate <- createSitePair(phototrophs_bray, lakes_nseqs_phototrophs,
                                               climate, geodist = FALSE)
phototrophs_sitepair_waterchem <- createSitePair(phototrophs_bray, lakes_nseqs_phototrophs,
                                                 waterchem, geodist = FALSE)

consumers_sitepair_geography <- createSitePair(consumers_bray, lakes_nseqs_consumers,
                                               geography %>%
                                                 dplyr::select(-latitude, -longitude), geodist = TRUE)
consumers_sitepair_morphometry <- createSitePair(consumers_bray, lakes_nseqs_consumers,
                                                 morphometry, geodist = FALSE)
consumers_sitepair_landuse_soil <- createSitePair(consumers_bray, lakes_nseqs_consumers,
                                                  landuse_soil, geodist = FALSE)
consumers_sitepair_climate <- createSitePair(consumers_bray, lakes_nseqs_consumers,
                                             climate, geodist = FALSE)
consumers_sitepair_waterchem <- createSitePair(consumers_bray, lakes_nseqs_consumers,
                                               waterchem, geodist = FALSE)

mixotrophs_sitepair_geography <- createSitePair(mixotrophs_bray, lakes_nseqs_mixotrophs,
                                                geography %>%
                                                  dplyr::select(-latitude, -longitude), geodist = TRUE)
mixotrophs_sitepair_morphometry <- createSitePair(mixotrophs_bray, lakes_nseqs_mixotrophs,
                                                  morphometry, geodist = FALSE)
mixotrophs_sitepair_landuse_soil <- createSitePair(mixotrophs_bray, lakes_nseqs_mixotrophs,
                                                   landuse_soil, geodist = FALSE)
mixotrophs_sitepair_climate <- createSitePair(mixotrophs_bray, lakes_nseqs_mixotrophs,
                                              climate, geodist = FALSE)
mixotrophs_sitepair_waterchem <- createSitePair(mixotrophs_bray, lakes_nseqs_mixotrophs,
                                                waterchem, geodist = FALSE)


#### Plot GDM fitted splines/functions ####
# Define function to extract, format, and plot I-splines
plotGDM <- function(sitepairs, envtable, geodist = FALSE) {
  # Filter study lakes in explanatory data
  envtable <- envtable %>%
    filter(lake_id %in% lakes) %>%
    dplyr::select(lake_id, where(is.numeric))
  
  gdm <- gdm(sitepairs, geo = geodist)
  deviance <- round(gdm$explained, 2)
  
  bind_cols(isplineExtract(gdm)$x %>%
              as_tibble() %>%
              pivot_longer(everything(), names_to = "predictor", values_to = "x"),
            isplineExtract(gdm)$y %>%
              as_tibble() %>%
              pivot_longer(everything(), names_to = "predictor", values_to = "y") %>%
              dplyr::select(y)) %>%
    ggplot() +
    facet_wrap(~predictor, scales = "free_x") +
    geom_point(aes(x = x, y = y, colour = predictor)) +
    geom_point(aes(x = value, y = -0.01),
               data = envtable %>%
                 pivot_longer(!lake_id, names_to = "predictor", values_to = "value") %>%
                 filter(predictor %in% gdm$predictors),
               alpha = 0.3) +
    xlab("Environmental gradient") +
    ylab("f(predictor)") +
    ggtitle(paste0(deviance, "% deviance explained")) +
    theme_light() %+replace%
    theme(legend.position = "none", strip.background = element_rect(fill = "black"))
}
# plotGDM(phototrophs_sitepair_geography, geography, geodist = TRUE)
# plotGDM(phototrophs_sitepair_morphometry, morphometry, geodist = FALSE)
# plotGDM(phototrophs_sitepair_landuse_soil, landuse, geodist = FALSE)
# plotGDM(phototrophs_sitepair_climate, climate, geodist = FALSE)
# plotGDM(phototrophs_sitepair_waterchem, waterchem, geodist = FALSE)

# plotGDM(consumers_sitepair_geography, geography, geodist = TRUE)
# plotGDM(consumers_sitepair_morphometry, morphometry, geodist = FALSE)
# plotGDM(consumers_sitepair_landuse_soil, landuse, geodist = FALSE)
# plotGDM(consumers_sitepair_climate, climate, geodist = FALSE)
# plotGDM(consumers_sitepair_waterchem, waterchem, geodist = FALSE)

# plotGDM(mixotrophs_sitepair_geography, geography, geodist = TRUE)
# plotGDM(mixotrophs_sitepair_morphometry, morphometry, geodist = FALSE)
# plotGDM(mixotrophs_sitepair_landuse_soil, landuse, geodist = FALSE)
# plotGDM(mixotrophs_sitepair_climate, climate, geodist = FALSE)
# plotGDM(mixotrophs_sitepair_waterchem, waterchem, geodist = FALSE)


#### Save site-pairs ####
# save(phototrophs_sitepair_geography, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_phototrophs_bray_geography.rda")
# save(phototrophs_sitepair_morphometry, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_phototrophs_bray_morphometry.rda")
# save(phototrophs_sitepair_landuse_soil, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_phototrophs_bray_landuse_soil.rda")
# save(phototrophs_sitepair_climate, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_phototrophs_bray_climate.rda")
# save(phototrophs_sitepair_waterchem, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_phototrophs_bray_waterchem.rda")

# save(consumers_sitepair_geography, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_consumers_bray_geography.rda")
# save(consumers_sitepair_morphometry, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_consumers_bray_morphometry.rda")
# save(consumers_sitepair_landuse_soil, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_consumers_bray_landuse_soil.rda")
# save(consumers_sitepair_climate, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_consumers_bray_climate.rda")
# save(consumers_sitepair_waterchem, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_consumers_bray_waterchem.rda")

# save(mixotrophs_sitepair_geography, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_mixotrophs_bray_geography.rda")
# save(mixotrophs_sitepair_morphometry, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_mixotrophs_bray_morphometry.rda")
# save(mixotrophs_sitepair_landuse_soil, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_mixotrophs_bray_landuse_soil.rda")
# save(mixotrophs_sitepair_climate, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_mixotrophs_bray_climate.rda")
# save(mixotrophs_sitepair_waterchem, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_mixotrophs_bray_waterchem.rda")
