library(tidyverse)
library(patchwork)
library(dplyr)
library(stringr)

# set working directory
setwd("/Users/natal/R/divergence")

# Load data ---------------------------------------------------------------

# read csv files
species <- read_csv("species_index_table.csv")
pairs <- read_csv("pair_index_table.csv")
divergence <- read_csv("divergence_table.csv")

convert_to_bases <- function(value) {
  value <- gsub(",", "", as.character(value)) # remove commas
  output <- rep(NA_real_, length(value))
  # convert Gb
  gb_mask <- str_detect(value, "Gb")
  output[gb_mask] <- suppressWarnings(as.numeric(str_replace(value[gb_mask], "Gb", ""))) * 1e9
  # convert Mb
  mb_mask <- str_detect(value, "Mb")
  output[mb_mask] <- suppressWarnings(as.numeric(str_replace(value[mb_mask], "Mb", ""))) * 1e6
  # convert kb
  kb_mask <- str_detect(value, "kb")
  output[kb_mask] <- suppressWarnings(as.numeric(str_replace(value[kb_mask], "kb", ""))) * 1e3
  # assume raw base values if numeric
  base_mask <- !gb_mask & !mb_mask & !kb_mask
  output[base_mask] <- suppressWarnings(as.numeric(value[base_mask]))
  return(output)
}


# convert N50 values
species <- species %>%
  mutate(
    N50_scaffold = convert_to_bases(N50_scaffold),
    N50_contig = convert_to_bases(N50_contig)
  )

# join species metadata for both species in a pair
df <- pairs %>%
  left_join(species, by = c("species_index1" = "Species_index")) %>%
  rename_with(~ str_replace(., "$", "1"), c("Species", "Accession", "Type", "Class", "Order", "Family",
                                            "Genome_size", "Body_size_kg", "Generation_time_yrs",
                                            "Repeat_percent", "N50_scaffold", "N50_contig",
                                            "Coverage", "Quality")) %>%
  left_join(species, by = c("species_index2" = "Species_index")) %>%
  rename_with(~ str_replace(., "$", "2"), c("Species", "Accession", "Type", "Class", "Order", "Family",
                                            "Genome_size", "Body_size_kg", "Generation_time_yrs",
                                            "Repeat_percent", "N50_scaffold", "N50_contig",
                                            "Coverage", "Quality")) %>%
  left_join(divergence, by = "pair") %>%
  # compute averaged values
  mutate(
    Genome_size_avg = if_else(is.na(Genome_size1) & is.na(Genome_size2), NA_real_, rowMeans(select(., Genome_size1, Genome_size2), na.rm = TRUE)),
    Body_size_kg_avg = if_else(is.na(Body_size_kg1) & is.na(Body_size_kg2), NA_real_, rowMeans(select(., Body_size_kg1, Body_size_kg2), na.rm = TRUE)),
    Generation_time_yrs_avg = if_else(is.na(Generation_time_yrs1) & is.na(Generation_time_yrs2), NA_real_, rowMeans(select(., Generation_time_yrs1, Generation_time_yrs2), na.rm = TRUE)),
    Repeat_percent_avg = if_else(is.na(Repeat_percent1) & is.na(Repeat_percent2), NA_real_, rowMeans(select(., Repeat_percent1, Repeat_percent2), na.rm = TRUE)),
    N50_scaffold_avg = if_else(is.na(N50_scaffold1) & is.na(N50_scaffold2), NA_real_, rowMeans(select(., N50_scaffold1, N50_scaffold2), na.rm = TRUE)),
    N50_contig_avg = if_else(is.na(N50_contig1) & is.na(N50_contig2), NA_real_, rowMeans(select(., N50_contig1, N50_contig2), na.rm = TRUE)),
    Coverage_avg = if_else(is.na(Coverage1) & is.na(Coverage2), NA_real_, rowMeans(select(., Coverage1, Coverage2), na.rm = TRUE))
  ) %>%
  # convert Timetree_div to numeric
  mutate(Timetree_div = if_else(Timetree_div == "Unknown", NA_character_, Timetree_div),
         Timetree_div = parse_number(Timetree_div)) %>%
  # deal with any other 0/NA/Unknown values
  mutate(across(c(raw, k3p, jc, COI_divergence, Timetree_div), ~ if_else(. == 0 | . == "Unknown" | is.na(.), NA_real_, as.numeric(.))))

# optional filter: Keep only pairs where both species have "Chromosome" quality
filter_chromosomes_only <- FALSE  # Set to TRUE to apply the filter
if (filter_chromosomes_only) {
  df <- df %>% filter(Quality1 == "Chromosome" & Quality2 == "Chromosome")
}

#exclude k2p values equal to zero
df <- df %>% filter(!is.na(k2p) & k2p != 0)

# view working df
glimpse(df)
#View(df)

# Histogram of k2p values -------------------------------------------------
histogram_plot <- ggplot(df, aes(x = k2p)) +
  geom_histogram(binwidth = 0.0001, fill = "skyblue", color = "black") +
  labs(
    title = "Histogram of k2p Divergence Values",
    x = "k2p Divergence",
    y = "Count"
  )

print(histogram_plot)

# Scatterplot of k2p by taxonomic order -----------------------------------

# filter out rows with NA in Order1
df_filtered <- df %>% 
  filter(!is.na(Order1))

# reorder Order1 factor so orders are grouped by class
order_levels <- df_filtered %>% 
  distinct(Order1, Class1) %>% 
  arrange(Class1, Order1) %>% 
  pull(Order1)

# update Order1 column to be a factor with the new levels
df_filtered <- df_filtered %>% 
  mutate(Order1 = factor(Order1, levels = order_levels))

# k2p Divergence vs Order with arranged Orders 
scatter_plot <- ggplot(df_filtered, aes(x = Order1, y = k2p, color = Class1)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "k2p Divergence by Order",
    x = "Order",
    y = "k2p Divergence",
    color = "Taxonomic Class"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(scatter_plot)

# Options to view specific datapoints --------------------------------------

# print Species Names for the 10 Lowest k2p Divergence Values
lowest10 <- df %>%
  arrange(k2p) %>%
  slice_head(n = 10) %>%
  select(pair, Species1, Species2, k2p)

print(lowest10)

# view data for a specific species
species_of_interest <- "Dipodomys_spectabilis"  # replace with your species of interest
divergence_by_species <- df %>%
  filter(Species1 == species_of_interest | Species2 == species_of_interest) %>%
  select(pair, Species1, Species2, raw, k2p, k3p, jc, COI_divergence)

print(divergence_by_species)

# view specific row
#View(df %>% filter(pair == 12))

# Plot k2p and alignment or repeat percent ------------------------------------------

# percent aligned 
k2p_aligned_plot <- ggplot(df_filtered, aes(x = alignment_percent, y = k2p, color = Class1)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title = "alignment percent vs k2p Divergence by Class",
    x = "alignment percent",
    y = "k2p Divergence",
    color = "Taxonomic Class"
  ) +
  theme_minimal() +
  facet_wrap(~ Class1, scales = "free")  # create a separate plot for each class

print(k2p_aligned_plot)

# repeat percent
k2p_repeat_plot <- ggplot(df_filtered, aes(x = Repeat_percent_avg, y = k2p, color = Class1)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title = "repeat percent vs k2p Divergence by Class",
    x = "repeat percent",
    y = "k2p Divergence",
    color = "Taxonomic Class"
  ) +
  theme_minimal() +
  facet_wrap(~ Class1, scales = "free")  # create a separate plot for each class

print(k2p_repeat_plot)
# Plot k2p and body size, generation time, or genome size -----------------

# genome size 
ggplot(df_filtered, aes(x = Genome_size_avg, y = k2p, color = Class1)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title = "Genome_size_avg vs k2p Divergence by Class",
    x = "Genome_size_avg",
    y = "k2p Divergence",
    color = "Taxonomic Class"
  ) +
  theme_minimal() +
  facet_wrap(~ Class1, scales = "free")  

# body size
ggplot(df_filtered, aes(x = Body_size_kg_avg, y = k2p, color = Class1)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title = "Body_size_kg_avg vs k2p Divergence by Class",
    x = "Body_size_kg_avg",
    y = "k2p Divergence",
    color = "Taxonomic Class"
  ) +
  theme_minimal() +
  facet_wrap(~ Class1, scales = "free")  

# generation time
ggplot(df_filtered, aes(x = Generation_time_yrs_avg, y = k2p, color = Class1)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title = "Generation_time_yrs_avg vs k2p Divergence by Class",
    x = "Generation_time_yrs_avg",
    y = "k2p Divergence",
    color = "Taxonomic Class"
  ) +
  theme_minimal() +
  facet_wrap(~ Class1, scales = "free")  

