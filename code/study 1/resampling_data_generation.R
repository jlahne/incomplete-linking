# This file does the resampling data generation for Study 1 in the manuscript

# Setup -------------------------------------------------------------------

# Import the helper functions (also loads packages)
library(here)
source(here("code/helper_functions.R"))

# Load the data from this and previous studies
link_table_incomplete_linking <- 
  read_csv(here("data/study 1/study_1_edgelist_incomplete_linking.csv"))
link_table_complete_linking <- 
  read_csv(here("data/study 1/study_1_edgelist_complete_linking.csv"))
link_table_complete_sorting <- 
  read_csv(here("data/study 1/study_1_edgelist_complete_sorting.csv"))
sample_ids <- 
  read_csv(here("data/study 1/study_1_samples.csv"))
ibd_design <-
  bind_rows(
    read_csv(here("data/study 1/Choc Block 1 Design.csv"), col_names = FALSE),
    read_csv(here("data/study 1/Choc Block 2 Design.csv"), col_names = FALSE)
  ) %>%
  as.matrix()


num_boots <- 1000 # set this to be ~100 for testing

set.seed(1234) # for replicability

# We pull a random sample of the complete linking and sorting sample sets in
# order to make the co-occurrence frequency (concurrency?) more similar to the
# incomplete block design.  This is still large enough for a typical
# sorting/linking study.
used_samples_complete_sort <- sample(x = 1:62, size = 20)
used_samples_complete_link <- sample(x = 1:63, size = 20)
link_table_complete_linking <- 
  link_table_complete_linking %>%
  filter(sample_set %in% paste0("J", used_samples_complete_link))
link_table_complete_sorting <- 
  link_table_complete_sorting %>%
  filter(sample_set %in% paste0("J", used_samples_complete_sort))


# Pairwise resampling -----------------------------------------------------

pairwise_resample_incomplete_linking <- 
  expand_grid(boot_id = 1:num_boots,
              ibd_row = 1:40) %>%
  mutate(simulated_graph = map(.x = ibd_row,
                               .f = ~simulate_graph(link_table_incomplete_linking,
                                                    selected_samples = ibd_design[.x, ],
                                                    as_sorting = FALSE)))

pairwise_resample_complete_linking <- 
  expand_grid(boot_id = 1:num_boots, 
              sample_set = 1:20) %>%
  mutate(simulated_graph = map(.x = sample_set,
                               .f = ~simulate_graph(link_table_complete_linking,
                                                    selected_samples = 1:10,
                                                    as_sorting = FALSE)))

pairwise_resample_complete_sorting <- 
  expand_grid(boot_id = 1:num_boots, 
              sample_set = 1:20) %>%
  mutate(simulated_graph = map(.x = sample_set,
                               .f = ~simulate_graph(link_table_complete_sorting,
                                                    selected_samples = 1:10,
                                                    as_sorting = TRUE)))

# Convert the graph results to partition groups

resampled_groups_incomplete_linking <-
  pairwise_resample_incomplete_linking %>%
  transmute(boot_id, 
            ibd_row,
            dissimilarity = map(.x = simulated_graph, 
                                .f = ~get_dissimilarity_from_graph(.x))) %>%
  unnest(dissimilarity) %>%
  nest(data = -boot_id) %>%
  transmute(boot_id,
            recursive_partition = map(.x = data,
                                      .f = ~get_dissimilarity_matrix(.x) %>%
                                        get_recursive_partition_groups()))

resampled_groups_complete_linking <-
  pairwise_resample_complete_linking %>%
  transmute(boot_id, 
            sample_set,
            dissimilarity = map(.x = simulated_graph, 
                                .f = ~get_dissimilarity_from_graph(.x))) %>%
  unnest(dissimilarity) %>%
  nest(data = -boot_id) %>%
  transmute(boot_id,
            recursive_partition = map(.x = data,
                                      .f = ~get_dissimilarity_matrix(.x) %>%
                                        get_recursive_partition_groups()))

resampled_groups_complete_sorting <-
  pairwise_resample_complete_sorting %>%
  transmute(boot_id, 
            sample_set,
            dissimilarity = map(.x = simulated_graph, 
                                .f = ~get_dissimilarity_from_graph(.x))) %>%
  unnest(dissimilarity) %>%
  nest(data = -boot_id) %>%
  transmute(boot_id,
            recursive_partition = map(.x = data,
                                      .f = ~get_dissimilarity_matrix(.x) %>%
                                        get_recursive_partition_groups()))

# Classic bootstrapping ---------------------------------------------------

# we can save calculation time by first converting the links to dissimilarities,
# and then bootstrapping (this is not possible in the pairwise resampling
# because that relies on simulation to generate NEW dissimilarities, not
# resampling the observed ones)

dissimilarity_table_complete_linking <-
  link_table_complete_linking %>%
  nest(data = -sample_set) %>%
  transmute(sample_set,
            graph_dissimilarity = map(.x = data,
                                      .f = ~get_dissimilarity_table(.x)))

dissimilarity_table_complete_sorting <-
  link_table_complete_sorting %>%
  nest(data = -sample_set) %>%
  transmute(sample_set,
            graph_dissimilarity = map(.x = data,
                                      .f = ~get_dissimilarity_table(.x)))

# then we do the classic bootstrap: resample sample_ids with replacement

bootstrapped_groups_complete_linking <- 
  dissimilarity_table_complete_linking %>%
  slice_sample(prop = num_boots, replace = TRUE) %>%
  transmute(boot_id = rep(1:num_boots, each = 20),
            graph_dissimilarity) %>%
  unnest(graph_dissimilarity) %>%
  nest(data = -boot_id) %>%
  transmute(boot_id,
            recursive_partition = map(.x = data,
                                      .f = ~get_dissimilarity_matrix(.x) %>% 
                                        get_recursive_partition_groups()))

bootstrapped_groups_complete_sorting <-
  dissimilarity_table_complete_sorting %>%
  slice_sample(prop = num_boots, replace = TRUE) %>%
  transmute(boot_id = rep(1:num_boots, each = 20),
            graph_dissimilarity) %>%
  unnest(graph_dissimilarity) %>%
  nest(data = -boot_id) %>%
  transmute(boot_id,
            recursive_partition = map(.x = data,
                                      .f = ~get_dissimilarity_matrix(.x) %>% 
                                        get_recursive_partition_groups()))


# Save the results for further analysis -----------------------------------

save(resampled_groups_complete_linking,
     resampled_groups_complete_sorting,
     resampled_groups_incomplete_linking,
     bootstrapped_groups_complete_linking,
     bootstrapped_groups_complete_sorting,
     pairwise_resample_complete_linking,
     pairwise_resample_complete_sorting,
     pairwise_resample_incomplete_linking,
     file = here("outputs/study 1/resampled_data_files.RData"))


























