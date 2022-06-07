# Replace missing LakePulse environmental data with ecozone medians

# Load libraries
library(tidyverse)
library(janitor)


#### Import and format data ####
# Import most current environmental dataset
metadata_date <- "2022-04-28"
metadata <- read_csv(paste0("data/environmental/lakepulse_data_curated_", metadata_date, ".csv"), col_names = TRUE) %>%
  rename(lake_id = lakepulse_id)

# Import list of saline lakes
(lakes_saline <- read_tsv("output/env/lakes_saline.tsv", col_names = "lake_id") %>%
    pull(lake_id))
length(lakes_saline)  # Number of saline lakes


#### Remove saline lakes from dataset ####
# Filter out saline lakes
metadata <- metadata %>%
  filter(!lake_id %in% lakes_saline)


#### Combine ecozones ####
# Combine Boreal Cordillera and Taiga Cordillera for the purpose of calculating ecozone medians
metadata <- metadata %>%
  mutate(ecozone_combined = case_when(ecozone == "Boreal Cordillera" | ecozone == "Taiga Cordillera" ~ "Boreal Cordillera/Taiga Cordillera",
                                      TRUE ~ ecozone))


#### Fill in missing chlorophyll-a AM data ####
metadata <- metadata %>%
  mutate(chla = case_when(is.na(chla_am) & !is.na(chla_pm) ~ chla_pm,
                          TRUE ~ chla_am)) %>%
  select(-chla_am, -chla_pm, -chla_day)


#### Calculate variable medians by ecozone ####
# No NAs in basic info, weather (ERA5 Land), or land use data so no need to correct medians
metadata_medians <- metadata %>%
  select(lake_id, ecozone_combined,
         volume, average_depth, discharge, residence_time, slope_100m,
         contains("rbr"),
         contains("do"),
         contains("ph"),
         chla, doc, dic, tss, colour,
         contains("spm"), contains("tp"), contains("tn"), contains("srp"),
         calcium, chloride, magnesium, potassium, sodium, sulfate,
         contains("secchi"),
         contains("bruntvaisala"), centerbuoyancy) %>%
  group_by(ecozone_combined) %>%
  summarise_if(is.numeric, median, na.rm = TRUE)


#### Replace NA values by ecozone medians ####
# Create variable summarizing which variables were median-corrected
metadata$median_corrected <- NA
for (i in 1:nrow(metadata)) {
  for (j in 1:ncol(metadata)) {
    if (colnames(metadata[i, j]) %in% colnames(metadata_medians[,-1]) & is.na(metadata[i, j])) {
      metadata[i, j] <- metadata_medians[which(metadata_medians$ecozone_combined == metadata$ecozone_combined[i]), colnames(metadata)[j]]
      
      var <- colnames(metadata[i, j])
      metadata$median_corrected[i] <- paste(metadata$median_corrected[i], var, sep = ", ")
    }
  }
}
metadata <- metadata %>%
  mutate(median_corrected = str_remove(median_corrected, "^NA, "))


#### Write median-corrected environmental data to file ####
# Include date stamp
# Uncomment to write file:
# metadata %>%
#   select(-ecozone_combined) %>%
#   arrange(lake_id) %>%
#   write_csv(paste0("output/env/lakepulse_data_curated_mediancorrected_nonsaline_", metadata_date, ".csv"), col_names = TRUE)
