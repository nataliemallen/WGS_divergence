# load libraries
library(tidyverse)
library(patchwork)

# set working directory
setwd("/Users/natal/R/divergence")

# load data -----------------------------------------------------------

# read csv files
species <- read_csv("species_index_table.csv")
pairs <- read_csv("pair_index_table.csv")
divergence <- read_csv("april_divergence_updated.csv")

# function to standardize N50 values
convert_to_bases <- function(value) {
  value <- gsub(",", "", as.character(value)) 
  output <- rep(NA_real_, length(value))
  gb_mask <- str_detect(value, "Gb")
  output[gb_mask] <- suppressWarnings(as.numeric(str_replace(value[gb_mask], "Gb", ""))) * 1e9
  mb_mask <- str_detect(value, "Mb")
  output[mb_mask] <- suppressWarnings(as.numeric(str_replace(value[mb_mask], "Mb", ""))) * 1e6
  kb_mask <- str_detect(value, "kb")
  output[kb_mask] <- suppressWarnings(as.numeric(str_replace(value[kb_mask], "kb", ""))) * 1e3
  base_mask <- !gb_mask & !mb_mask & !kb_mask
  output[base_mask] <- suppressWarnings(as.numeric(value[base_mask]))
  return(output)
}

# convert N50 
species <- species %>%
  mutate(
    N50_scaffold = convert_to_bases(N50_scaffold),
    N50_contig = convert_to_bases(N50_contig)
  )

# flag any rows where N50_scaffold < N50_contig 
species <- species %>%
  mutate(flag_N50_issue = if_else(N50_scaffold < N50_contig, TRUE, FALSE, missing = FALSE))

# join species metadata for each pair
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
  # averaged values
  mutate(
    Genome_size_avg = if_else(is.na(Genome_size1) & is.na(Genome_size2), NA_real_, rowMeans(select(., Genome_size1, Genome_size2), na.rm = TRUE)),
    Body_size_kg_avg = if_else(is.na(Body_size_kg1) & is.na(Body_size_kg2), NA_real_, rowMeans(select(., Body_size_kg1, Body_size_kg2), na.rm = TRUE)),
    Generation_time_yrs_avg = if_else(is.na(Generation_time_yrs1) & is.na(Generation_time_yrs2), NA_real_, rowMeans(select(., Generation_time_yrs1, Generation_time_yrs2), na.rm = TRUE)),
    Repeat_percent_avg = if_else(is.na(Repeat_percent1) & is.na(Repeat_percent2), NA_real_, rowMeans(select(., Repeat_percent1, Repeat_percent2), na.rm = TRUE)),
    N50_scaffold_avg = if_else(is.na(N50_scaffold1) & is.na(N50_scaffold2), NA_real_, rowMeans(select(., N50_scaffold1, N50_scaffold2), na.rm = TRUE)),
    N50_contig_avg = if_else(is.na(N50_contig1) & is.na(N50_contig2), NA_real_, rowMeans(select(., N50_contig1, N50_contig2), na.rm = TRUE)),
    Coverage_avg = if_else(is.na(Coverage1) & is.na(Coverage2), NA_real_, rowMeans(select(., Coverage1, Coverage2), na.rm = TRUE))
  )



#temporarily exclude k2p values equal to zero
df <- df %>% filter(!is.na(k2p) & k2p != 0)

# optional
## exclude reptiles
# df <- df %>%
#   filter(Class1 != "Reptilia" & Class2 != "Reptilia")
# ## combine fish
# df <- df %>%
#   mutate(Class1 = case_when(
#     Class1 %in% c("Actinopterygii", "Osteichthyes") ~ "Fish",
#     TRUE ~ Class1
#   ))

glimpse(df)
#View(df)

# histograms ---------------------------------------------------------------
# k20
k2p_histogram <- ggplot(df, aes(x = k2p)) +
  geom_histogram(binwidth = 0.0001, fill = "skyblue", color = "black") +
  labs(
    title = "Histogram of k2p Divergence Values",
    x = "k2p Divergence",
    y = "Count"
  )

print(k2p_histogram)

# N50
n50_histogram <- ggplot(df, aes(x = N50_scaffold_avg)) +
  geom_histogram(binwidth = 1000000, fill = "skyblue", color = "black") +
  labs(
    title = "Histogram of scaffold N50",
    x = "Avg scaffold N50",
    y = "Count"
  )

print(n50_histogram)

# view specific datapoints ------------------------------------------------

# print species names for lowest k2p values
lowest10 <- df %>%
  arrange(k2p) %>%        
  slice_head(n = 10) %>%  
  select(pair, Species1, Species2, k2p)

print(lowest10)

# see divergence for a specific species
species_of_interest <- "Dipodomys_spectabilis"  
divergence_by_species <- df %>%
  filter(Species1 == species_of_interest | Species2 == species_of_interest) %>%
  select(pair, Species1, Species2, raw, k2p, k3p, jc)

print(divergence_by_species)

# k2p boxplots ------------------------------------------------

ggplot(df, aes(x = Class1, y = k2p, fill = Class1)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(alpha = 0.4, color = "black", width = 0.2) + 
  labs(x = "Class", y = "Divergence (k2p)", title = "Divergence by Class") +
  theme_minimal() +
  theme(text = element_text(size = 15)) +
  scale_fill_viridis_d(option = "D") 

ggplot(df, aes(x = Order1, y = k2p, fill = Class1)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.4, color = "black", width = 0.2) + 
  theme_minimal() +
  labs(x = "Order", y = "k2p divergence", title = "k2p divergence by Order (colored by Class)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# k2p Scatter Plots by Order -----------------------------------------------

# create order levels for plot arrangement
df_filtered <- df %>% 
  filter(!is.na(Order1))

order_levels <- df_filtered %>%
  distinct(Order1, Class1) %>%
  arrange(Class1, Order1) %>%
  distinct(Order1, .keep_all = TRUE) %>%
  pull(Order1) %>%
  as.character()

order_levels <- as.character(order_levels)

df_filtered <- df_filtered %>%
  mutate(Order1 = factor(Order1, levels = order_levels))

levels(df_filtered$Order1)

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

# divergence variable plots -----------------------------------------------

df_filtered_divergence <- df %>% 
  filter(!is.na(k2p))

# k2p and genome size
ggplot(df_filtered_divergence, aes(x = Genome_size_avg, y = k2p, color = Class1)) +
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

# k2p and body size
ggplot(df_filtered_divergence, aes(x = Body_size_kg_avg, y = k2p, color = Class1)) +
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

# k2p and generation time
ggplot(df_filtered_divergence, aes(x = Generation_time_yrs_avg, y = k2p, color = Class1)) +
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

# k2p and repeat content
ggplot(df_filtered_divergence, aes(x = Repeat_percent_avg, y = k2p, color = Class1)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title = "Repeat_percent_avg vs k2p Divergence by Class",
    x = "Repeat_percent_avg",
    y = "k2p Divergence",
    color = "Taxonomic Class"
  ) +
  theme_minimal() +
  facet_wrap(~ Class1, scales = "free")  

# log transformation for body size
ggplot(df_filtered_divergence, aes(x = Body_size_kg_avg, y = k2p, color = Class1)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_x_log10() +  
  labs(
    title = "Body Size vs k2p Divergence by Class",
    x = "Body Size (kg, log scale)",
    y = "k2p Divergence",
    color = "Taxonomic Class"
  ) +
  theme_minimal() +
  facet_wrap(~ Class1, scales = "free")

# extract pairs for timetree ----------------------------------------------

## only need to run this once

# create the new dataframe
df_output <- df %>%
  mutate(
    Genus = sub("_.*", "", Species1),  
    Species1 = sub(".*_", "", Species1), 
    Species2 = sub(".*_", "", Species2)   
  ) %>%
  select(Genus, Species1, Species2) %>%  
  mutate(
    TTID1 = NA,  
    TTID2 = NA,  
    DIVERGENCE = NA  
  )

write_csv(df_output, "species_pairs_timetree_april.csv")

head(df_output)

View(df)

# add timetree values to df -----------------------------------------------

# read the timetree divergence data
timedivergence <- read_csv("species_pairs_timetree_april_updated.csv")

# create species names with genus to match df
timedivergence <- timedivergence %>%
  mutate(Species1_full = paste(Genus, Species1, sep = "_"),  
         Species2_full = paste(Genus, Species2, sep = "_")) 

# create timetree column in dr
df <- df %>%
  left_join(timedivergence, by = c("Species1" = "Species1_full", "Species2" = "Species2_full")) %>%
  mutate(Timetree_div = if_else(is.na(DIVERGENCE) | DIVERGENCE == "Unknown", "Unknown", DIVERGENCE)) %>%
  select(-DIVERGENCE)  

df$Timetree_div <- as.numeric(df$Timetree_div)

head(df)


# k2p and timetree plots --------------------------------------------------------

df_filtered_divergence <- df %>% 
  filter(!is.na(k2p) & !is.na(Timetree_div))

df_filtered_divergence <- df_filtered_divergence %>%
  filter(Timetree_div != "Unknown")

# k2p genome divergence vs timetree divergence, colored by Class
time_vs_k2p_plot <- ggplot(df_filtered_divergence, aes(x = Timetree_div, y = k2p, color = Class1)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "Timetree Divergence vs k2p Divergence by Class",
    x = "Timetree Divergence",
    y = "k2p Divergence",
    color = "Taxonomic Class"
  ) +
  theme_minimal()

print(time_vs_k2p_plot)

# with trend lines
time_vs_k2p_plot <- ggplot(df_filtered_divergence, aes(x = Timetree_div, y = k2p, color = Class1)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) + 
  labs(
    title = "Timetree Divergence vs k2p Divergence by Class",
    x = "Timetree Divergence",
    y = "k2p Divergence",
    color = "Taxonomic Class"
  ) +
  theme_minimal()

print(time_vs_k2p_plot)

# unified trend line
time_vs_k2p_plot <- ggplot(df_filtered_divergence, aes(x = Timetree_div, y = k2p, color = Class1)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +  
  labs(
    title = "Timetree Divergence vs k2p Divergence by Class",
    x = "Timetree Divergence",
    y = "k2p Divergence",
    color = "Taxonomic Class"
  ) +
  theme_minimal()

print(time_vs_k2p_plot)

# with panels
time_vs_k2p_plot <- ggplot(df_filtered_divergence, aes(x = Timetree_div, y = k2p, color = Class1)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    title = "Timetree Divergence vs k2p Divergence by Class",
    x = "Timetree Divergence",
    y = "k2p Divergence",
    color = "Taxonomic Class"
  ) +
  theme_minimal() +
  facet_wrap(~ Class1, scales = "free")  

print(time_vs_k2p_plot)
