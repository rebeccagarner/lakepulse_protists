# Identify saline lakes

# Load libraries
library(tidyverse)
library(janitor)


#### Import and format data ####
# Import environmental data that has already had missing data filled by ecozone medians
# Rename urbanization variable as "built"
metadata <- read_csv("data/environmental/lakepulse_data_curated_2022-04-28.csv", col_names = TRUE) %>%
  rename(lake_id = lakepulse_id) %>%
  rename(area_built = area_urban, fraction_built = fraction_urban)


#### Categorize lakes by trophic status ####
# Categorize trophic states by total phosphorus trigger ranges for Canadian lakes and rivers:
# - Ultraoligotrophic (TP <4 μg/L)
# - Oligotrophic (4-10 μg/L)
# - Mesotrophic (10-20 μg/L)
# - Mesoeutrophic (20-35 μg/L)
# - Eutrophic (35-100 μg/L)
# - Hypereutrophic (>100 μg/L)
metadata <- metadata %>%
  mutate(trophic_status = case_when(tp_tube < 4 ~ "ultraoligotrophic",
                                    tp_tube >= 4 & tp_tube < 10 ~ "oligotrophic",
                                    tp_tube >= 10 & tp_tube < 20 ~ "mesotrophic",
                                    tp_tube >= 20 & tp_tube < 35 ~ "mesoeutrophic",
                                    tp_tube >= 35 & tp_tube < 100 ~ "eutrophic",
                                    tp_tube >= 100 ~ "hypereutrophic"))


#### Categorize lakes by salt status ####
# Extract salinity variables
salinity_vars <- metadata %>%
  dplyr::select(lake_id, lake_name, trophic_status, tp_tube,
                rbr_conductivity_tube, rbr_specific_conductance_tube,
                rbr_salinity_tube,
                calcium, chloride, magnesium, potassium, sodium, sulfate)

# Visualize relationship between conductivity and specific conductance
salinity_vars %>%
  ggplot(aes(x = rbr_conductivity_tube, y = rbr_specific_conductance_tube, colour = rbr_salinity_tube)) +
  geom_point(alpha = 0.7) +
  scale_colour_gradient2(low = "lightblue", mid = "lightgrey", high = "red", midpoint = 10) +
  theme_classic()

# Calculate total ions (sum major anion and cation concentrations in mg/L)
# Categorize lakes by specific conductivity:
# - Fresh (less than 80 μs cm-1),
# - Oligosaline (800–8,000 μs cm-1),
# - Esosaline (8,000–30,000 μs cm-1),
# - Polysaline (30,000–45,000 μs cm-1),
# - Eusaline (45,000–60,000 μs cm-1), and
# - Hypersaline (greater than 60,000 μs cm-1) 
salinity_classes <- salinity_vars %>%
  mutate(salt_status_conductivity = case_when((rbr_conductivity_tube * 1000) < 800 ~ "fresh",
                                              (rbr_conductivity_tube * 1000) >= 800 & (rbr_conductivity_tube * 1000) < 8000 ~ "oligosaline",
                                              (rbr_conductivity_tube * 1000) >= 8000 & (rbr_conductivity_tube * 1000) < 30000 ~ "esosaline",
                                              (rbr_conductivity_tube * 1000) >= 30000 & (rbr_conductivity_tube * 1000) < 45000 ~ "polysaline",
                                              (rbr_conductivity_tube * 1000) >= 45000 & (rbr_conductivity_tube * 1000) < 60000 ~ "eusaline",
                                              (rbr_conductivity_tube * 1000) >= 60000 ~ "hypersaline")) %>%
  mutate(ions_total = calcium + chloride + magnesium + potassium + sodium + sulfate)

# Plot histogram of salinity frequencies
salinity_classes %>%
  pivot_longer(!c(lake_id, lake_name, trophic_status, salt_status_conductivity),
               names_to = "var", values_to = "value") %>%
  filter(var %in% c("rbr_conductivity_tube", "ions_total")) %>%
  ggplot(aes(x = value, y = ..count..,
             fill = factor(salt_status_conductivity, levels = c("fresh", "oligosaline", "esosaline", "polysaline")))) +
  facet_wrap(~var, scales = "free") +
  geom_histogram(bins = 20) +
  scale_fill_manual(values = c("fresh" = "lightblue", "oligosaline" = "yellow", "esosaline" = "orange", "polysaline" = "pink"),
                    na.value = "black") +
  labs(fill = "Salt status") +
  theme_linedraw()

# Plot relationship between total ions and conductivity
salinity_classes %>%
  ggplot(aes(x = rbr_conductivity_tube,
             y = ions_total,
             colour = factor(salt_status_conductivity, levels = c("fresh", "oligosaline", "esosaline", "polysaline")))) +
  geom_point() +
  geom_hline(yintercept = 4000, linetype = "dashed", alpha = 0.3) +
  scale_colour_manual(values = c("fresh" = "lightblue", "oligosaline" = "yellow", "esosaline" = "orange", "polysaline" = "pink"),
                      na.value = "black") +
  labs(colour = "Salt status") +
  theme_classic()


#### Identify saline lakes ####
# Identify saline lakes
lakes_saline <- salinity_classes %>%
  filter((rbr_conductivity_tube * 1000) >= 8000 | ions_total >= 4000) %>%
  arrange(lake_id) %>%
  dplyr::select(lake_id)

# Write list of saline lakes to file
# lakes_saline %>%
#   write_tsv("output/env/lakes_saline.tsv", col_names = FALSE)
