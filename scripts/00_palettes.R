# Define palettes for figures

# Ecozones ----
palette_ecozone_letters <- c(72, 65, 66, 98, 86, 77, 109, 88, 80, 83, 84, 116)
names(palette_ecozone_letters) <- c("Atlantic Highlands", "Atlantic Maritime", "Boreal Cordillera",
                                    "Boreal Plains", "Boreal Shield", "Mixedwood Plains",
                                    "Montane Cordillera", "Pacific Maritime", "Prairies",
                                    "Semi-Arid Plateaux", "Taiga Cordillera", "Taiga Plains")

palette_ecozone_spectral <- c("#E44A33", "#F17C4A", "#FDAE61",
                              "#FEC980", "#FFE4A0", "#FFFFBF",
                              "#E3F3CD", "#C7E6DB", "#ABD9E9",
                              "#81BAD8", "#569AC7", "#2C7BB6")
names(palette_ecozone_spectral) <- c("Taiga Cordillera", "Boreal Cordillera", "Montane Cordillera",
                                     "Pacific Maritime", "Taiga Plains", "Semi-Arid Plateaux",
                                     "Boreal Plains", "Prairies", "Mixedwood Plains",
                                     "Boreal Shield", "Atlantic Highlands", "Atlantic Maritime")


# Environmental factors ----
# Trophic state ----
palette_trophic <- c("#016DCB", "#8DDDFD", "#00FFD0", "#40F99B", "#EEFA41", "#FA9141")
names(palette_trophic) <- c("ultraoligotrophic", "oligotrophic", "mesotrophic", "mesoeutrophic", "eutrophic", "hypereutrophic")

# Explanatory variable categories
palette_environment <- c("#984447", "#ADD9F4", "#476C9B", "#F7A1C4", "#7BD389", "#FDE74C")
names(palette_environment) <- c("Watershed", "Lake morphometry/physical", "Water chemistry",
                                "Spatial MEMs", "Geography", "Climate")


# Taxonomy ----
# Opisthokonts, land plants, and putative protists
palette_taxgroup <- c("#FF2D56", "#D81159", "#00B26E", "yellow")
names(palette_taxgroup) <- c("metazoa", "fungi", "plants", "protists")

# Eukaryotic supergroups
palette_supergroup <- c("#75DAFF", "#FF2D56", "#3DFFE5", "#FF7CEB", "#FF7CEB", "#00E288",
                        "#F6FF00", "#8A84E2", "#E8E8E8", "#CC7D5D")
names(palette_supergroup) <- c("Stramenopiles", "Opisthokonta", "Hacrobia", "Alveolata", "Protalveolata", "Archaeplastida",
                               "Rhizaria", "Excavata", "Amoebozoa", "Apusozoa")

# Eukaryotic divisions
palette_division <- c("#75DAFF", "#276ED8", "#54B2FF", "#C9F7FF", "#1D4692",
                      "#FF2D56", "#D81159", "#FF8787", "#DD6E6E",
                      "#3DFFE5", "#0099AA", "#9BFFEB", "#31E2DF", "#9EE2E1",
                      "#FF7CEB", "#EA00C7", "#FF007B", "#CE00B6", "#FFDAF5",
                      "#00E288", "#00B26E", "#8BE271", "#D2FF91",
                      "#F6FF00", "#FFE100",
                      "#8A84E2", "#AFAFDC", "#6D26E1",
                      "#E8E8E8", "#B7B7B7", "#808080",
                      "#CC7D5D", "#9B6E5D")
names(palette_division) <- c("Ochrophyta", "Opalozoa", "Pseudofungi", "Sagenista", "Stramenopiles_X",
                             "Metazoa", "Fungi", "Choanoflagellida", "Mesomycetozoa",
                             "Cryptophyta", "Katablepharidophyta", "Centroheliozoa", "Haptophyta", "Telonemia",
                             "Dinoflagellata", "Ciliophora", "Apicomplexa", "Perkinsea", "Protalveolata_X",
                             "Chlorophyta", "Streptophyta", "Rhodophyta", "Glaucophyta",
                             "Cercozoa", "Foraminifera",
                             "Discoba", "Metamonada", "Malawimonadidae",
                             "Lobosa", "Conosa", "Breviatea",
                             "Hilomonadea", "Apusomonadidae")


# Trophic functional groups ----
palette_function <- c("#F55D3E", "#B78BD4", "#76BED0", "#56E4C0", "#62BE81",  # orange soda, lavender floral, dark sky blue, medium aquamarine, emerald
                      "#EDD55C", "#FDFF80", "#AB6B42")  # arylide yellow, laser lemon, brown sugar
names(palette_function) <- c("parasites", "consumers, parasites", "consumers", "consumers, mixotrophs", "mixotrophs",
                             "phototrophs, mixotrophs", "phototrophs", "phototrophs, consumers, mixotrophs")
