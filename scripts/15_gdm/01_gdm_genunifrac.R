# Create site-pairs for generalized dissimilarity modelling (GDM)
# https://cran.r-project.org/web/packages/gdm/vignettes/gdmVignette.pdf

# Load libraries
library(tidyverse)
library(vegan)
library(gdm)
library(janitor)


#### Import and format data ####
# Import (generalized alpha = 0.5) UniFrac distance matrix
load("output/dissim/lp2017-2019watercolumn18s_protists_unifrac_gen0pt5alpha.rda")

# Import melted protist ASV table
asv_melt <- read_tsv("output/melted/lp2017-2019watercolumn18s_melt_protists_366lakes.tsv")
length(unique(asv_melt$asv_code))

# Import metadata
source("scripts/05_metadata.R")


#### Create site-pair table ####
# Calculate n sequences per lake to weight site-pairs (to account for biases in sampling effort)
lakes_nseqs <- asv_melt %>%
  group_by(lake_id) %>%
  dplyr::summarize(weights = sum(nseqs)) %>%
  arrange(lake_id) %>%
  as.data.frame()

# Combine response (community composition) data and predictor (environmental) data into site-pair table format
# Define function to create GDM site-pair table
createSitePair <- function(dist, envtable, geodist = FALSE) {
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
sitepair_geography_phylo <- createSitePair(unifrac_gen0.5alpha_dist, geography %>%
                                       dplyr::select(-latitude, -longitude), geodist = TRUE)
sitepair_morphometry_phylo <- createSitePair(unifrac_gen0.5alpha_dist, morphometry, geodist = FALSE)
sitepair_landuse_soil_phylo <- createSitePair(unifrac_gen0.5alpha_dist, landuse_soil, geodist = FALSE)
sitepair_climate_phylo <- createSitePair(unifrac_gen0.5alpha_dist, climate, geodist = FALSE)
sitepair_waterchem_phylo <- createSitePair(unifrac_gen0.5alpha_dist, waterchem, geodist = FALSE)


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
# plotGDM(sitepair_geography_phylo, geography, geodist = TRUE)
# plotGDM(sitepair_morphometry_phylo, morphometry, geodist = FALSE)
# plotGDM(sitepair_landuse_soil_phylo, landuse, geodist = FALSE)
# plotGDM(sitepair_climate_phylo, climate, geodist = FALSE)
# plotGDM(sitepair_waterchem_phylo, waterchem, geodist = FALSE)


#### Save site-pairs ####
# save(sitepair_geography_phylo, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_protists_genunifrac_geography.rda")
# save(sitepair_morphometry_phylo, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_protists_genunifrac_morphometry.rda")
# save(sitepair_landuse_soil_phylo, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_protists_genunifrac_landuse_soil.rda")
# save(sitepair_climate_phylo, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_protists_genunifrac_climate.rda")
# save(sitepair_waterchem_phylo, file = "output/gdm/sitepairs/lp2017-2019watercolumn18s_sitepair_protists_genunifrac_waterchem.rda")
