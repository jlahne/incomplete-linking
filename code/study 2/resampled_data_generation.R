# This script imports the cleaned data from the Cocoa Wheel of Excellence
# incomplete-block design study and then uses pairwise resampling and classic
# bootstrapping to generate measures for estimating stability of the designs.


# Load resources ----------------------------------------------------------

library(here)
source(here("code/helper_functions.R"))

clean_incomplete_linking <-
  read_csv(here("data/study 2/clean_incomplete_linking_table.csv"))

clean_incomplete_sorting <-
  read_csv(here("data/study 2/clean_incomplete_sorting_table.csv"))

clean_complete_sorting <-
  read_csv(here("data/study 2/clean_complete_sorting_table.csv"))

ibd_design <-
  read_csv(here("data/study 2/ibd_design.csv")) %>%
  as.matrix()

# Examine cooccurrence ----------------------------------------------------

### Incomplete Linking:

# Total sample sets (of 3 blocks) = 30

clean_incomplete_linking %>%
  separate(col = sample_set, into = c("design", "rep", "subject")) %>%
  count(subject)

# Presentation count minimum = 21, maximum = 24
clean_incomplete_linking %>%
  count(from, to) %>%
  filter(from == to) %>%
  pull(n) %>%
  summary()

# Minimum co-occurrence = 3, maximum = 8

clean_incomplete_linking %>%
  count(from, to) %>%
  filter(from != to) %>%
  pull(n) %>%
  summary()

### Incomplete sorting:

# Total sample sets (of 3 blocks) = 58

clean_incomplete_sorting %>%
  separate(sample_set, into = c("subject", "rep")) %>%
  count(subject)

# Presentation count minimum = 41, maximum = 48

clean_incomplete_sorting %>%
  count(from, to) %>%
  filter(from == to) %>%
  pull(n) %>%
  summary()

# Minimum co-occurrence = 6, maximum = 16

clean_incomplete_sorting %>%
  count(from, to) %>%
  filter(from != to) %>%
  pull(n) %>%
  summary()

# Minimum/maximum co-occurrence = 53

clean_complete_sorting %>%
  count(from, to) %>%
  pull(n) %>%
  summary()

# To make the experiment more "fair" we will pull 30 random sample sets from the
# complete sorting study

set.seed(1234) # for reproducibility

clean_complete_sorting <-
  clean_complete_sorting %>%
  nest(data = -sample_set) %>%
  slice_sample(n = 30) %>%
  unnest(everything())

# And we will pull the first 30 subjects from the incomplete sorting study
# This results in a presentation count minimum = 21, maximum = 24
# And co-occurrence min = 3, max = 24

clean_incomplete_sorting <-
  clean_incomplete_sorting %>%
  nest(data = -sample_set) %>%
  separate(sample_set, into = c("judge", "block"), sep = "_") %>%
  mutate(judge = as.numeric(judge), block = as.numeric(block)) %>%
  arrange(judge, block) %>%
  slice_head(n = 90) %>%
  unite(judge, block, col = "sample_set") %>%
  unnest(cols = data)

# Resample data -----------------------------------------------------------

tic()
pb <- progress::progress_bar$new(total = 4)

num_boots <- 1000 # set to 10-100 for testing

incomplete_linking_resampled_pairwise <-
  expand_grid(boot_id = 1:num_boots,
              ibd_row = 1:93) %>%
  mutate(simulated_graph = map(.x = ibd_row,
                               .f = ~simulate_graph(
                                 original_results = clean_incomplete_linking,
                                 selected_samples = ibd_design[.x, ],
                                 as_sorting = FALSE
                               )))

pb$tick()

# incomplete sorting

incomplete_sorting_resampled_pairwise <-
  expand_grid(boot_id = 1:num_boots,
              ibd_row = 1:93) %>%
  mutate(simulated_graph = map(.x = ibd_row,
                               .f = ~simulate_graph(
                                 original_results = clean_incomplete_sorting,
                                 selected_samples = ibd_design[.x, ],
                                 as_sorting = TRUE
                               )))

pb$tick()

# complete sorting - pairwise

complete_sorting_resampled_pairwise <-
  expand_grid(boot_id = 1:num_boots,
              ibd_row = 1:31) %>%
  mutate(simulated_graph = map(.x = ibd_row,
                               .f = ~simulate_graph(
                                 original_results = clean_incomplete_linking,
                                 selected_samples = 1:62,
                                 as_sorting = TRUE
                               )))

pb$tick()

# complete sorting - classic

complete_sorting_resampled_classic <-
  tibble(boot_id = 1:num_boots) %>%
  mutate(sorting_rc_groups = map(.x = boot_id,
                                 .f = ~nest(clean_complete_sorting, data = -sample_set) %>%
                                   slice_sample(prop = 1, replace = TRUE) %>%
                                   mutate(dissimilarity = map(.x = data, .f = ~get_dissimilarity_table(.x))) %>%
                                   unnest(dissimilarity) %>%
                                   get_dissimilarity_matrix() %>%
                                   get_recursive_partition_groups()))

pb$tick()

save(incomplete_linking_resampled_pairwise,
     incomplete_sorting_resampled_pairwise,
     complete_sorting_resampled_classic,
     complete_sorting_resampled_pairwise,
     file = here("outputs/study 2/chocolate wheel resampled data.RData"))

toc()
