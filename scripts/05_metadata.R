# Curate lake environmental metadata

# Load libraries
library(tidyverse)
library(janitor)
library(geoR)


#### Import metadata ####
# Import lakes list
lakes <- read_tsv("output/env/lakes.tsv", col_names = "lake_id") %>%
  pull(lake_id)
length(lakes)  # Number of lakes

# Import environmental data that has already had missing data filled by ecozone medians
# Rename urbanization variable as "built"
metadata <- read_csv("output/env/lakepulse_data_curated_mediancorrected_nonsaline_2022-04-28.csv", col_names = TRUE) %>%
  rename(area_built = area_urban, fraction_built = fraction_urban)

# Import soil grid properties data and extract data for surfacemost soil depth
# https://www.isric.org/explore/soilgrids/faq-soilgrids
# Name / Description / Mapped units / Conversion factor / Conventional units
# bdod	Bulk density of the fine earth fraction	cg/cm³	100	kg/dm³
# cec	Cation Exchange Capacity of the soil	mmol(c)/kg	10	cmol(c)/kg
# cfvo	Volumetric fraction of coarse fragments (> 2 mm)	cm3/dm3 (vol‰)	10	cm3/100cm3 (vol%)
# clay	Proportion of clay particles (< 0.002 mm) in the fine earth fraction	g/kg	10	g/100g (%)
# nitrogen	Total nitrogen (N)	cg/kg	100	g/kg
# phh2o	Soil pH	pHx10	10	pH
# sand	Proportion of sand particles (> 0.05 mm) in the fine earth fraction	g/kg	10	g/100g (%)
# silt	Proportion of silt particles (≥ 0.002 mm and ≤ 0.05 mm) in the fine earth fraction	g/kg	10	g/100g (%)
# soc	Soil organic carbon content in the fine earth fraction	dg/kg	10	g/kg
# ocd	Organic carbon density	hg/dm³	10	kg/dm³
# ocs	Organic carbon stocks	t/ha	10	kg/m²
soil_grid <- read_csv("data/environmental/soil_grid_profile_mean_20210201.csv", col_names = TRUE) %>%
  rename(lake_id = lakepulse_id) %>%
  dplyr::select(lake_id, contains("0_5"))

metadata <- metadata %>%
  left_join(soil_grid, by = "lake_id")


#### Lake selection and identification ####
# Filter lakes
metadata <- metadata %>%
  filter(lake_id %in% lakes)
nrow(metadata) == length(lakes)  # Should evaluate to TRUE

# Order ecozones by longitude (west to east)
ecozone_longitude_order <- metadata %>%
  dplyr::select(lake_id, ecozone, latitude, longitude) %>%
  group_by(ecozone) %>%
  summarise_if(is.numeric, median, na.rm = TRUE) %>%
  arrange(longitude) %>%
  pull(ecozone)


#### Categorize lakes by trophic status ####
# Categorize trophic states by total phosphorus trigger ranges for Canadian lakes and rivers:
# - Ultraoligotrophic (TP <4 ug/L)
# - Oligotrophic (4-10 ug/L)
# - Mesotrophic (10-20 ug/L)
# - Mesoeutrophic (20-35 ug/L)
# - Eutrophic (35-100 ug/L)
# - Hypereutrophic (>100 ug/L)
metadata <- metadata %>%
  mutate(trophic_status = case_when(tp_tube < 4 ~ "ultraoligotrophic",
                                    tp_tube >= 4 & tp_tube < 10 ~ "oligotrophic",
                                    tp_tube >= 10 & tp_tube < 20 ~ "mesotrophic",
                                    tp_tube >= 20 & tp_tube < 35 ~ "mesoeutrophic",
                                    tp_tube >= 35 & tp_tube < 100 ~ "eutrophic",
                                    tp_tube >= 100 ~ "hypereutrophic"))


#### Plot data distributions ####
ecozones <- metadata %>%
  dplyr::select(lake_id, ecozone)

# Define ecozone palette
palette_ecozone_spectral <- c("#E44A33", "#F17C4A", "#FDAE61",
                              "#FEC980", "#FFE4A0", "#FFFFBF",
                              "#E3F3CD", "#C7E6DB", "#ABD9E9",
                              "#81BAD8", "#569AC7", "#2C7BB6")
names(palette_ecozone_spectral) <- c("Taiga Cordillera", "Boreal Cordillera", "Montane Cordillera",
                                     "Pacific Maritime", "Taiga Plains", "Semi-Arid Plateaux",
                                     "Boreal Plains", "Prairies", "Mixedwood Plains",
                                     "Boreal Shield", "Atlantic Highlands", "Atlantic Maritime")

# Curate facet labels
variable_labels <- c("Geographic distance",
                     
                     "Latitude (\u00B0N)", "Longitude (\u00B0W)", "Latitude (\u00B0N)", "Longitude (\u00B0W)",
                     "Altitude (m)",
                     
                     "Air temp (\u00B0C)", "Precipitations (m)", "Solar radiation (J/m2)",
                     "Wind speed (m/s)", "Ice disappearance (day)",
                     
                     "Area (km2)", "Circularity", "Watershed slope (\u00B0)",
                     "Volume (mcm)", "Max depth (m)", "Discharge (m3/s)",
                     "Residence time (days)",
                     "Watershed area (km2)", "Lake : watershed area",
                     
                     "Crop agriculture (%)", "Forestry (%)", "Natural landscapes (%)",
                     "Pasture (%)", "Built (%)", "Population (people/km2)",
                     
                     "Fine earth density (100 kg/dm3)", "Cation exchange capacity (cmol(c)/kg)", "Coarse fragments (vol%)",
                     "Clay (%)", "Soil nitrogen (g/kg)", "Soil organic carbon density (kg/dm3)",
                     "Soil pH", "Sand (%)", "Silt (%)",
                     "Fine earth organic carbon (g/kg)",
                     
                     "Surface temp (\u00B0C)",
                     "Chlorophyll-a (\u03BCg/L)", "DIC (mg/L)", "DOC (mg/L)",
                     "pH", "Colour (mg/L Pt)",
                     "Total nitrogen (mg/L)", "Total phosphorus (\u03BCg/L)", "Calcium (mg/L)",
                     "Chloride (mg/L)", "Magnesium (mg/L)", "Potassium (mg/L)",
                     "Sodium (mg/L)", "Sulfate (mg/L)")

names(variable_labels) <- c("Geographic",
                            
                            "latitude", "longitude", "latitude_y", "longitude_x",
                            "altitude",
                            
                            "temperature_mean_7d", "precipitation_total_7d", "solar_radiation_net_7d",
                            "windspeed_mean_7d", "ice_disappearance_julianday",
                            
                            "area", "circularity", "slope_100m",
                            "volume", "max_depth", "discharge",
                            "residence_time",
                            "watershed_area", "lakewatershed_area_ratio",
                            
                            "fraction_agriculture", "fraction_forestry", "fraction_natlandscapes",
                            "fraction_pasture", "fraction_built", "population_density",
                            
                            "bdod_mean_0_5", "cec_mean_0_5", "cfvo_mean_0_5",
                            "clay_mean_0_5", "nitrogen_mean_0_5", "ocd_mean_0_5",
                            "phh2o_mean_0_5", "sand_mean_0_5", "silt_mean_0_5",
                            "soc_mean_0_5",
                            
                            "rbr_temperature_mean_tube",
                            "chla", "dic", "doc",
                            "ph_epilimnion", "colour",
                            "tn_tube", "tp_tube", "calcium",
                            "chloride", "magnesium", "potassium",
                            "sodium", "sulfate")

# Curate variable labels without units
variables <- c("Geographic distance",
               
               "Latitude", "Longitude", "Latitude", "Longitude",
               "Altitude",
               
               "Air temp", "Precipitations", "Solar radiation",
               "Wind speed", "Ice disappearance",
               
               "Area", "Circularity", "Watershed slope",
               "Volume", "Max depth", "Discharge",
               "Residence time",
               "Watershed area", "Lake:watershed area",
               
               "Crop agriculture", "Forestry", "Natural landscapes",
               "Pasture", "Built", "Population",
               
               "Fine earth density", "Cation exchange capacity", "Coarse fragments",
               "Clay", "Soil N", "Soil organic C density",
               "Soil pH", "Sand", "Silt",
               "Fine earth organic C",
               
               "Surface temp",
               "Chlorophyll-a", "DIC", "DOC",
               "pH", "Colour",
               "TN", "TP", "Calcium",
               "Chloride", "Magnesium", "Potassium",
               "Sodium", "Sulfate")


names(variables) <- c("Geographic",
                      
                      "latitude", "longitude", "latitude_y", "longitude_x",
                      "altitude",
                      
                      "temperature_mean_7d", "precipitation_total_7d", "solar_radiation_net_7d",
                      "windspeed_mean_7d", "ice_disappearance_julianday",
                      
                      "area", "circularity", "slope_100m",
                      "volume", "max_depth", "discharge",
                      "residence_time",
                      "watershed_area", "lakewatershed_area_ratio",
                      
                      "fraction_agriculture", "fraction_forestry", "fraction_natlandscapes",
                      "fraction_pasture", "fraction_built", "population_density",
                      
                      "bdod_mean_0_5", "cec_mean_0_5", "cfvo_mean_0_5",
                      "clay_mean_0_5", "nitrogen_mean_0_5", "ocd_mean_0_5",
                      "phh2o_mean_0_5", "sand_mean_0_5", "silt_mean_0_5",
                      "soc_mean_0_5",
                      
                      "rbr_temperature_mean_tube",
                      "chla", "dic", "doc",
                      "ph_epilimnion", "colour",
                      "tn_tube", "tp_tube", "calcium",
                      "chloride", "magnesium", "potassium",
                      "sodium", "sulfate")

variable_labels_nounits <- tibble(variable = names(variables),
                                  variable_curated = unname(variables))

# Define function to plot data distributions
plotDistrib <- function(data, plot_ncols = 8) {
  ecozone_medians <- data %>%
    pivot_longer(!lake_id, names_to = "variable", values_to = "value") %>%
    left_join(ecozones, by = "lake_id") %>%
    group_by(ecozone, variable) %>%
    dplyr::summarize(median = median(value))
  
  data_long <- data %>%
    pivot_longer(!lake_id, names_to = "variable", values_to = "value") %>%
    left_join(ecozones, by = "lake_id")
  
  ggplot() +
    geom_histogram(aes(x = value, y = ..count.., fill = factor(ecozone, levels = ecozone_longitude_order)),
                   data = data_long, bins = 30) +
    geom_vline(aes(xintercept  = median, colour = factor(ecozone, levels = ecozone_longitude_order)),
               data = ecozone_medians, linetype = "dashed") +
    facet_wrap(~variable, scales = "free",
               labeller = labeller(variable = variable_labels),
               ncol = plot_ncols) +
    scale_fill_manual(values = palette_ecozone_spectral, name = "Ecozone") +
    scale_colour_manual(values = palette_ecozone_spectral, name = "Ecozone") +
    labs(y = "Number of lakes") +
    theme_light() %+replace%
    theme(axis.title.x = element_blank(),
          strip.background = element_rect(fill = "black"), strip.text.x = element_text(face = "bold"))
}


#### Evaluate collinearity in metadata ####
# Define function to plot correlations
plotCorr <- function(data) {
  # Define function to return upper triangle of distance matrix
  halveDist <- function(mat) {
    for (i in 1:ncol(mat)) {
      for (j in 1:i) {
        mat[i,j] <- NA
      }
    }
    
    return(mat)
  }
  
  data %>%
    dplyr::select(where(is.numeric)) %>%
    dplyr::select(order(colnames(.))) %>%
    cor() %>%
    halveDist() %>%
    as_tibble(rownames = "var1") %>%
    pivot_longer(!var1, names_to = "var2", values_to = "correlation") %>%
    filter(!is.na(correlation)) %>%
    ggplot() +
    geom_tile(aes(x = var1, y = var2, fill = correlation)) +
    geom_text(aes(x = var1, y = var2, label = round(correlation, 2))) +
    scale_fill_gradient2(low = "red", mid = "white", high = "dodgerblue",
                         midpoint = 0, na.value = "white", limits = c(-1,1)) +
    scale_x_discrete(position = "top") +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0),
          panel.background = element_blank(), axis.ticks = element_blank())
}


#### Format environmental data ####
# Format geography data ----
# Format latitude-longitude coordinates
latlong <- metadata %>%
  dplyr::select(lake_id, latitude, longitude)

# Format geography data
geography <- metadata %>%
  mutate(latitude_y = latitude,
         longitude_x = longitude) %>%
  dplyr::select(lake_id, latitude, longitude,
                latitude_y, longitude_x, altitude)

# geography %>%
#   dplyr::select(lake_id, latitude, longitude, altitude) %>%
#   plotDistrib()
# geography %>%
#   dplyr::select(lake_id, latitude, longitude, altitude) %>%
#   plotCorr()

# Format morphometry data ----
morphometry <- metadata %>%
  dplyr::select(lake_id, area, circularity, shoreline_length, slope_100m,
                volume, max_depth, average_depth, discharge, residence_time,
                watershed_area) %>%
  mutate(lakewatershed_area_ratio = area/watershed_area)

morphometry$shoreline_length <- NULL  # Remove shoreline length because highly correlated with area
morphometry$average_depth <- NULL  # Remove average depth because highly correlated with max depth, slope, and volume
# plotDistrib(morphometry)
# plotCorr(morphometry)


# Format land use data ----
# Select land use variables
landuse_fractions <- metadata %>%
  dplyr::select(lake_id, fraction_agriculture, fraction_pasture, fraction_forestry,
         fraction_built, fraction_natlandscapes, fraction_water)

# Recalculate land use fractions for terrestrial area of watershed (i.e. subtract water)
# Combine agriculture and pasture fractions
# Convert fractions to percentages
landuse_terrestrial <- landuse_fractions %>%
  mutate(fraction_all = fraction_agriculture +
           fraction_pasture +
           fraction_forestry +
           fraction_built +
           fraction_natlandscapes +
           fraction_water) %>%
  #mutate(fraction_agri_pasture = fraction_agriculture + fraction_pasture) %>%
  #dplyr::select(-fraction_agriculture, -fraction_pasture) %>%
  mutate(fraction_terrestrial = fraction_all - fraction_water) %>%
  dplyr::select(!fraction_water) %>%
  pivot_longer(!c(lake_id, fraction_all, fraction_terrestrial),
               names_to = "landuse", values_to = "fraction") %>%
  mutate(fraction_landuse_terrestrial = fraction/fraction_terrestrial) %>%
  dplyr::select(-fraction_all, -fraction_terrestrial, -fraction) %>%
  mutate(pct_landuse_terrestrial = fraction_landuse_terrestrial * 100) %>%
  dplyr::select(lake_id, landuse, pct_landuse_terrestrial) %>%
  pivot_wider(names_from = landuse, values_from = pct_landuse_terrestrial)
sum(landuse_terrestrial[,-1])/100 == length(unique(metadata$lake_id))  # Should evaluate to TRUE

metadata <- metadata %>%
  dplyr::select(-names(landuse_terrestrial)[which(grepl(pattern = "^fraction_", names(landuse_terrestrial)))]) %>%
  left_join(landuse_terrestrial, by = "lake_id")

landuse <- metadata %>%
  dplyr::select(lake_id,
                watershed_area,
                population,
                fraction_agriculture, fraction_forestry, fraction_natlandscapes,
                fraction_pasture, fraction_built) %>%
  mutate(population_density = population/watershed_area) %>%
  dplyr::select(-population, -watershed_area)

# plotDistrib(landuse)
# plotCorr(landuse)

# Format soil data ----
soil <- metadata %>%
  dplyr::select(lake_id,
                ends_with("_0_5"))

# plotDistrib(soil)
# plotCorr(soil)

# Combine land use and soil data
landuse_soil <- landuse %>%
  left_join(soil, by = "lake_id")

# Format climate/weather data ----
climate <- metadata %>%
  dplyr::select(lake_id,
                temperature_mean_7d, temperature_mean_30d,
                precipitation_total_7d, precipitation_total_30d,
                solar_radiation_net_7d, solar_radiation_net_30d,
                windspeed_mean_7d, windspeed_mean_30d,
                ice_disappearance_julianday)

plotCorr(climate)  # Remove 30-day weather data
climate <- climate %>%
  dplyr::select(-contains("30d"))

# plotDistrib(climate)
# plotCorr(climate)

# Format water chemistry data ----
waterchem <- metadata %>%
  dplyr::select(lake_id,
                rbr_temperature_mean_tube,
                chla, dic, doc, ph_epilimnion, colour,
                tn_tube, tp_tube,
                calcium, chloride, magnesium, potassium, sodium, sulfate)

# plotDistrib(waterchem)
# plotCorr(waterchem)


#### Combine all environmental data ####
env_all <- latlong %>%
  left_join(geography, by = c("lake_id", "latitude", "longitude")) %>%
  left_join(morphometry, by = "lake_id") %>%
  left_join(landuse, by = "lake_id") %>%
  left_join(soil_grid, by = "lake_id") %>%
  left_join(climate, by = "lake_id") %>%
  left_join(waterchem, by = "lake_id")

# plotCorr(env_all %>%
#            dplyr::select(-latitude, -longitude))

# Write curated metadata to file
# metadata %>%
#   write_tsv("output/metadata.tsv", col_names = TRUE)
