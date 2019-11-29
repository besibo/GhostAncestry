# ---------------------------------------------------------------------------------
# --------------------------------- All scaffolds ---------------------------------
# ---------------------------------------------------------------------------------

# Here, we try to extend what we did for one specific scaffold to the whole list of scaffold of more than 1Mbase

# ----------------------------------------------------------------------------------
# ------------------------------------ Preamble ------------------------------------
# ----------------------------------------------------------------------------------
# Here, I load all necessary packages, functions, and data
# And I specify several parameters about the scaffold of interest and out files

# -------- Load useful functions --------
# Packages
library(tidyverse)
library(purrr)
library(furrr)

# Custom functions I created that are needed to generate the outfiles
source("Custom_functions.R")
# ---------------------------------------



# ----- Import necessary data files -----
## list of autosomes of length > 1Mbp
autosomes <- read_table2("2_All scaffolds/_infiles/autosomes_sup1Mbp_length.txt", col_names = FALSE)
colnames(autosomes) <- c("chr", "size")

## haplo file
haplo <- read_table2("2_All scaffolds/_infiles/NWAc_ghost_allscaffolds.haplo")

## .bed file
bed <- read_table2("2_All scaffolds/_infiles/mask_sort_3_merged.bed", col_names = FALSE)
colnames(bed) <- c("chr", "start", "end")

## Remove scaffold shorter than 1MB in the bed file
bed <- autosomes %>% 
  left_join(bed) %>% 
  select(-size)
# ---------------------------------------



# -------------- Parameters -------------
## Specify the width of the window (bp)
window_width <- 1000

## Specify the width of the window to compute values of mu by region
mu_width <- 100000 # Window size for the mutation rate computations

## Add the lower bound of the last window for all scaffold
## and add a table with the lower bound of all windows for all scaffold
autosomes <- autosomes %>% 
  mutate(low_bound_last_window = size - (size %% window_width),
         window_inf = map(low_bound_last_window, ~seq(0, ., by = window_width)))

## Unnest this table to get more than 2 million row
autosomes_unnested  <- autosomes %>%
  unnest(window_inf)


# Join the autosomes and bed tables to keep "bed information" 
# only for autosomes of length greater than 1Mbp
#bed_1M <- inner_join(autosomes, bed)

# We have 116 scaffolds of more than 1Mbp
# length(unique(bed_1M$chr))
# ---------------------------------------



# ---------- Output file names ----------
obs_file_name      <- "2_All scaffolds/_outfiles/All_Scaffolds_obs.txt"
weight_file_name   <- "2_All scaffolds/_outfiles/All_Scaffolds_weight.txt"
mutation_file_name <- "2_All scaffolds/_outfiles/All_Scaffolds_mut.txt"
# ---------------------------------------


# ------- Use multiple processors -------
plan(multiprocess)
# ---------------------------------------


# ----------------------------------------------------------------------------------
# -------------------------------- Observation file --------------------------------
# ----------------------------------------------------------------------------------

# 1. Filter positions for which ind0 has a private mutation
res <- haplo %>% 
  unite("other", ind1:ind21, sep = "", 
        remove = FALSE) %>%                      # Create a string (called "other") by concatenating the nucleotide of all individuals except ind0
  filter(!str_detect(other, ind0)) %>%           # Remove all positions for which the base of ind0 is detected in "other"
  mutate(n_base = n_distinct_bases(other)) %>%   # Count the number of distinct bases in "other" (excluding Ns)
  filter(n_base < 2) %>%                         # Keep the positions for which only one base (excluding Ns) was found in "other"
  select(-other, -n_base)                        # Remove temporary columns

# 2. Create a new table in which the rows contain:
#   - the bottom size of each window
#   - the number of private mutations in each window
#   - the list of private mutations in each window
pos_list <- res %>% 
  select(chr, pos) %>% 
  mutate(window_inf = pos - (pos %% window_width)) %>% 
  nest(data = c(pos)) %>% 
  mutate(n = future_map_int(data, ~nrow(.)),
         list = future_map_chr(data, ~str_flatten(pull(.), " "))) %>% 
  select(-data)

# 3. Since pos_list doesn't contain all possible windows (but only those for which
# private mutations are observed), we merge it with all_windows_lb to include the 
# missing windows and we tidy the final output by replacing NAs with appropriate values
pos_list_full <- autosomes_unnested %>% 
  left_join(pos_list) %>% 
  mutate(n = replace_na(n, 0),
         list = replace_na(list, ""))

# 4. Finally, write this table on disk and make sure
# no scientific notation will popup in the outfile
pos_list_full %>% 
  select(chr, window_inf, n, list) %>% 
  mutate(window_inf = format(window_inf, scientific = FALSE)) %>% 
  write.table(obs_file_name, quote = FALSE, row.names = FALSE, col.names = FALSE)


# ---------------------------------------------------------------------------------
# ---------------------------------- Weight file ----------------------------------
# ---------------------------------------------------------------------------------

# 1. Filter data for scaffolds of more than 1Mb and
# add 2 new columns to the bed file indicating the bottom values of the
# windows for which the uncallable regions start (div_start) and end (div_end) 
bed_lb <- bed %>% 
  mutate(window_start = start - (start %% window_width),
         window_end   = end   - (end   %% window_width))

# 2. Apply the add_lines() function to each row of bed_one_scaf in order
# to separate into multiple rows the uncallable regions spanning
# multiple windows. 
# Warning: this takes a few minutes on my computer (3-6).
bed_windows <- bed_lb %>% 
  mutate(new = future_pmap(., add_lines)) %>% 
  select(new) %>% 
  unnest_legacy(new)

# 3. compute the fraction of callable bases for each window in bed_windows:
# - first, tag each line of the table with the lower bound of the appropriate window 
#   (several rows can indeed fall inside the same window).
# - then, compute the total number of uncallable bases for each window
# - finally, transform these numbers as a fraction of callable bases
bed_clean <- bed_windows %>% 
  mutate(window_inf = start - (start %% window_width),
         uncall = end - start) %>% 
  group_by(chr, window_inf) %>% 
  summarise(uncall = sum(uncall)) %>% 
  mutate(callable = (window_width - uncall) / window_width) %>% 
  select(-uncall)

# 4. Manually add the missing windows and assign them a callable rate of 1
# This step is necessary because some windows may still be missing in bed_clean.
# indeed, if a window is 100% callable it won't appear in the original .bed file. 
bed_clean <- autosomes_unnested %>% 
  left_join(bed_clean) %>% 
  mutate(callable = replace_na(callable, 1)) %>% 
  select(chr, window_inf, callable)

# 5. Finally, write this table on disk and make sure
# no scientific notation will popup in the outfile
bed_clean %>% 
  mutate(window_inf = format(window_inf, scientific = FALSE)) %>% 
  write.table(weight_file_name, quote = FALSE, row.names = FALSE, col.names = FALSE)



# ---------------------------------------------------------------------------------
# --------------------------------- Mutation file ---------------------------------
# ---------------------------------------------------------------------------------

# 1. Compute mu_glob: 
#   - count the number of mutated sites in haplo for each scaffold
#   - add the length of the scaffold found in "autosomes"
#   - divide the count by the length
mu_scaffold <- haplo %>% 
  count(chr) %>% 
  left_join(autosomes) %>% 
  mutate(mu_glob = n / size,
         last_region = size %% mu_width)

# 2. Compute mu_reg
#   - For each region of 100,000 bp (or any other length), count the number 
#     of mutated positions and divide by the region size
mu_window <- haplo %>% 
  select(chr, pos) %>% 
  mutate(lb_win_mu = pos - (pos %% mu_width)) %>% 
  group_by(chr, lb_win_mu) %>% 
  summarise(n_mut = n(),
            mu_reg = n_mut / mu_width)

#    - Replace the last value: we have to divide by the true size of the last region
#      since it is (almost) always smaller than mu_width
mu_window <- mu_window %>% 
  nest(data = c(-chr)) %>% 
  left_join(mu_scaffold) %>% 
  mutate(data2 = map2(.x = data, .y = last_region, ~replace_last(.x, .y, 3))) %>% 
  select(chr, data2) %>% 
  unnest_legacy()
  
# 3. Compute mu_relat by dividing each mu_reg by mu_glob
mu_relat <- mu_window %>% 
  left_join(mu_scaffold, by = "chr") %>% 
  select(-n, -size) %>% 
  mutate(mu_relat = mu_reg / mu_glob)

# 4. create the final table by duplicating mu_relat for each 1000bp window along the scaffold
mu_table <- autosomes_unnested %>% 
  mutate(lb_win_mu = window_inf - (window_inf %% mu_width)) %>% 
  left_join(mu_relat %>% select(-window_inf)) %>% 
  select(chr, window_inf, mu_relat)

# 5. Write that table on disk and make sure
# no scientific notation will popup in the outfile
mu_table %>% 
  mutate(window_inf = format(window_inf, scientific = FALSE)) %>% 
  write.table(mutation_file_name, quote = FALSE, row.names = FALSE, col.names = FALSE)

