# This file analyzes both the results and the resampled results for study 1

# Setup -------------------------------------------------------------------

# Import the helper functions (also loads packages)
library(here)
source(here("code/helper_functions.R"))
library(flextable)

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

# Load the resampled data generated from the "resampling data generation.R"
# script
load(here("outputs/study 1/resampled_data_files.RData"))

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

# Cartoons of resampling approach -----------------------------------------

plot_sort_cartoon <- function(data = link_table_incomplete_linking,
                              design = sample(x = 1:10, size = 6, replace = FALSE),
                              sort = FALSE, edge_colors = c("black", "grey"),
                              edge_types = c("dashed", "solid"),
                              edge_width = c(1, 0.2)){

  data <-
    data %>%
    select(from, to, presence)

  data %>%
    graph_from_data_frame(directed = FALSE) %>%
    as_tbl_graph() %>%
    graph_join(
      data %>%
        simulate_graph2(selected_samples = design, as_sorting = sort) %>%
        activate(edges) %>%
        mutate(resample = "selected") %>%
        activate(nodes) %>%
        mutate(resample = "selected")
    ) %>%
    activate(edges) %>%
    mutate(presence = as.factor(presence),
           resample = ifelse(is.na(resample), "unselected", resample)) %>%
    activate(nodes) %>%
    mutate(resample = ifelse(is.na(resample), "unselected", resample)) %>%
    ggraph(layout = "circle") +
    geom_edge_fan(aes(linetype = presence, color = resample, width = resample)) +
    geom_node_point(aes(color = resample), size = 10, fill = "white", shape = 21) +
    geom_node_text(aes(label = name, color = resample)) +
    coord_fixed() +
    theme_graph() +
    scale_edge_linetype_manual(values = edge_types) +
    scale_edge_width_manual(values = edge_width) +
    scale_edge_color_manual(values = edge_colors) +
    scale_color_manual(values = edge_colors) +
    theme(legend.position = "none")
}

# The complete design
p_observed_similarities <-
  plot_sort_cartoon(design = 1:10, edge_colors = c("black", "black"),
                    edge_width = c(0.1, 0.1)) +
  labs(caption = "Study 1 (incomplete linking) multigraph of observed responses")

set.seed(123)
p_simulated_similarities_1_1 <-
  plot_sort_cartoon(design = ibd_design[1,]) +
  labs(caption = "Simulated draw 1 (resample) from Study 1, IBD row 1 based on observed responses")

set.seed(123)
p_simulated_similarities_1_1_selected <-
  plot_sort_cartoon(design = ibd_design[1,],
                    edge_colors = c("black", "transparent"),
                    edge_types = c("blank", "solid")) +
  labs(caption = "Resample 1 from Study 1, IBD row 1 with only retained 'similar' judgments")

set.seed(456)
p_simulated_similarities_1_2 <-
  plot_sort_cartoon(design = ibd_design[1,]) +
  labs(caption = "Simulated draw 2 (resample) from Study 1, IBD row 1 based on observed responses")

set.seed(456)
p_simulated_similarities_1_2_selected <-
  plot_sort_cartoon(design = ibd_design[1,],
                    edge_colors = c("black", "transparent"),
                    edge_types = c("blank", "solid")) +
  labs(caption = "Resample 2 from Study 1, IBD row 1 with only retained 'similar' judgments")

set.seed(123)
plot_sort_cartoon(design = 1:10, sort = FALSE)
set.seed(123)
plot_sort_cartoon(design = 1:10, sort = FALSE,
                  edge_colors = c("black", "transparent"),
                  edge_types = c("blank", "solid"))

ggsave(filename = "img/study_1_full_multigraph.png",
       plot = p_observed_similarities,
       height = 3.5, width = 3.5, units = "in", scale = 2.5, dpi = 300)
ggsave(filename = "img/study_1_ibd_row_1_possibilities_1.png",
       plot = p_simulated_similarities_1_1,
       height = 2.5, width = 2.5, units = "in", scale = 2.7, dpi = 300)
ggsave(filename = "img/study_1_ibd_row_1_resample_1.png",
       plot = p_simulated_similarities_1_1_selected,
       height = 2.5, width = 2.5, units = "in", scale = 2.7, dpi = 300)
ggsave(filename = "img/study_1_ibd_row_1_possibilities_2.png",
       plot = p_simulated_similarities_1_2,
       height = 2.5, width = 2.5, units = "in", scale = 2.7, dpi = 300)
ggsave(filename = "img/study_1_ibd_row_1_resample_2.png",
       plot = p_simulated_similarities_1_2_selected,
       height = 2.5, width = 2.5, units = "in", scale = 2.7, dpi = 300)


# Grouping analysis by recursive partitioning -----------------------------

# Incomplete linking

diss_matrix_incomplete_linking <-
  link_table_incomplete_linking %>%
  nest(data = -id) %>%
  transmute(dissimilarity = map(.x = data,
                                .f = ~get_dissimilarity_table(.x))) %>%
  unnest(dissimilarity) %>%
  get_dissimilarity_matrix()

recursive_partition_groups_incomplete_linking <-
  diss_matrix_incomplete_linking %>%
  get_recursive_partition_groups()

additive_tree_incomplete_linking <-
  diss_matrix_incomplete_linking %>%
  dist(method = "max") %>%
  nj() %>%
  as_tbl_graph() %>%
  left_join(recursive_partition_groups_incomplete_linking,
            by = c("name" = "sample")) %>%
  left_join(sample_ids %>% mutate(sample = as.character(sample)),
            by = c("name" = "sample"))

# Complete linking

diss_matrix_complete_linking <-
  link_table_complete_linking %>%
  nest(data = -sample_set) %>%
  transmute(dissimilarity = map(.x = data,
                                .f = ~get_dissimilarity_table(.x))) %>%
  unnest(dissimilarity) %>%
  get_dissimilarity_matrix()

recursive_partition_groups_complete_linking <-
  diss_matrix_complete_linking %>%
  get_recursive_partition_groups()

additive_tree_complete_linking <-
  diss_matrix_complete_linking %>%
  dist(method = "max") %>%
  nj() %>%
  as_tbl_graph() %>%
  left_join(recursive_partition_groups_complete_linking,
            by = c("name" = "sample"))

# Complete sorting

diss_matrix_complete_sorting <-
  link_table_complete_sorting %>%
  nest(data = -sample_set) %>%
  transmute(dissimilarity = map(.x = data,
                                .f = ~get_dissimilarity_table(.x))) %>%
  unnest(dissimilarity) %>%
  get_dissimilarity_matrix()

recursive_partition_groups_complete_sorting <-
  diss_matrix_complete_sorting %>%
  get_recursive_partition_groups()

additive_tree_complete_sorting <-
  diss_matrix_complete_sorting %>%
  dist(method = "max") %>%
  nj() %>%
  as_tbl_graph() %>%
  left_join(recursive_partition_groups_complete_sorting,
            by = c("name" = "sample"))

# And now we plot these

p_tree_incomplete_link <-
  additive_tree_incomplete_linking %>%
  mutate(type = ifelse(str_detect(name, "Node"), "internal", "leaf"),
         name.y = str_replace_all(name.y, "_", " ") %>%
           str_remove_all("\\?")) %>%
  mutate(group = as.factor(group)) %>%
  activate(edges) %>%
  mutate(color = .N()$group[to]) %>%
  ggraph(layout = "unrooted", length = length) +
  geom_edge_link(aes(color = color), show.legend = FALSE) +
  geom_node_text(data = . %>% filter(type == "leaf"),
                 aes(label = name.y,
                     color = group,
                     angle = atan(y / x) * 180 / pi,
                     hjust = ifelse(x > 0, "left", "right")),
                 show.legend = FALSE) +
  coord_equal() +
  theme_graph() +
  expand_limits(x = c(-13, 10), y = c(-15, 12)) +
  labs(caption = "incomplete linking")

p_tree_complete_link <-
  additive_tree_complete_linking %>%
  mutate(type = ifelse(str_detect(name, "Node"), "internal", "leaf"),
         name = str_replace_all(name, "_", " ") %>%
           str_remove_all("\\?")) %>%
  mutate(group = as.factor(group)) %>%
  activate(edges) %>%
  mutate(color = .N()$group[to]) %>%
  ggraph(layout = "unrooted", length = length) +
  geom_edge_link(aes(color = color), show.legend = FALSE) +
  geom_node_text(data = . %>% filter(type == "leaf"),
                 aes(label = name,
                     color = group,
                     angle = atan(y / x) * 180 / pi,
                     hjust = ifelse(x > 0, "left", "right")),
                 show.legend = FALSE) +
  coord_equal() +
  theme_graph() +
  expand_limits(x = c(-30, 20), y = c(-30, 20)) +
  labs(caption = "complete linking")

p_tree_complete_sort <-
  additive_tree_complete_sorting %>%
  mutate(type = ifelse(str_detect(name, "Node"), "internal", "leaf"),
         name = str_replace_all(name, "_", " ") %>%
           str_remove_all("\\?")) %>%
  mutate(group = as.factor(group)) %>%
  activate(edges) %>%
  mutate(color = .N()$group[to]) %>%
  ggraph(layout = "unrooted", length = length) +
  geom_edge_link(aes(color = color), show.legend = FALSE) +
  geom_node_text(data = . %>% filter(type == "leaf"),
                 aes(label = name,
                     color = group,
                     angle = atan(y / x) * 180 / pi,
                     hjust = ifelse(x > 0, "left", "right")),
                 show.legend = FALSE) +
  coord_equal() +
  theme_graph() +
  expand_limits(x = c(-40, 20), y = c(-30, 20)) +
  labs(caption = "complete sorting")

ggsave(filename = "img/study_1_tree_complete_sort.png", plot = p_tree_complete_sort,
       height = 6, width = 6, units = "in", scale = 1.2, dpi = 300)
ggsave(filename = "img/study_1_tree_incomplete_link.png", plot = p_tree_incomplete_link,
       height = 6, width = 6, units = "in", scale = 1.2, dpi = 300)
ggsave(filename = "img/study_1_tree_complete_link.png", plot = p_tree_complete_link,
       height = 6, width = 6, units = "in", scale = 1.2, dpi = 300)

# Graph statistics --------------------------------------------------------

link_table_combined <-
  bind_rows(
    link_table_complete_linking %>%
      mutate(study = "complete linking"),
    link_table_complete_sorting %>%
      mutate(study = "complete sorting"),
    link_table_incomplete_linking %>%
      left_join(sample_ids, by = c("from" = "sample")) %>%
      left_join(sample_ids, by = c("to" = "sample")) %>%
      transmute(sample_set = id,
                from = name.x,
                to = name.y,
                presence,
                study = "incomplete_linking")
  )

graph_statistics <-
  link_table_combined %>%
  nest(data = -c(study, sample_set)) %>%
  mutate(graph = map(.x = data, ~graph_from_links(.x))) %>%
  transmute(study, sample_set,
            nodes = map(.x = graph,
                        .f = ~V(.x)$name),
            transitivity = map(.x = graph,
                               .f = ~transitivity(.x, type = "local", isolate = "zero")),
            degree = map(.x = graph,
                         .f = ~igraph::degree(.x, normalized = TRUE)),
            components = map(.x = graph,
                             .f = ~component_count(.x))) %>%
  unnest(-c(study, sample_set))

ft_graph_stats <-
  graph_statistics %>%
  pivot_longer(transitivity:components) %>%
  group_by(study, name) %>%
  ggdist::median_qi(value) %>%
  mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  transmute(study, name,
            CI = str_c(value, " (", .lower, ", ", .upper, ")")) %>%
  pivot_wider(names_from = name, values_from = CI) %>%
  flextable()

ft_graph_stats %>%
  save_as_docx(path = "tables/study_1_graph_stats.docx")

# Jaccard stability (resampling analysis) ---------------------------------

# Unfortunately this is processor-intensive, and takes about 20-30 minutes to
# run.

# Incomplete linking (pairwise resampling by necessity)

observed_partitions_incomplete_linking <-
  recursive_partition_groups_incomplete_linking %>%
  recursive_partition_pairwise()

jaccard_stability_incomplete_linking <-
  pairwise_resample_incomplete_linking %>%
  transmute(boot_id,
            dissimilarites = map(.x = simulated_graph,
                                 .f = ~get_dissimilarity_from_graph(.x))) %>%
  unnest(everything()) %>%
  nest(data = -boot_id) %>%
  transmute(boot_id,
            jaccard_similarity = map(.x = data,
                                     .f = ~get_dissimilarity_matrix(.x) %>%
                                       get_recursive_partition_groups() %>%
                                       recursive_partition_pairwise() %>%
                                       jaccard_index(observed_partitions_incomplete_linking, .))) %>%
  unnest(everything())

# Complete linking (pairwise)

observed_partitions_complete_linking <-
  recursive_partition_groups_complete_linking %>%
  recursive_partition_pairwise()

jaccard_stability_pairwise_complete_linking <-
  pairwise_resample_complete_linking %>%
  transmute(boot_id,
            dissimilarites = map(.x = simulated_graph,
                                 .f = ~get_dissimilarity_from_graph(.x))) %>%
  unnest(everything()) %>%
  nest(data = -boot_id) %>%
  transmute(boot_id,
            jaccard_similarity = map(.x = data,
                                     .f = ~get_dissimilarity_matrix(.x) %>%
                                       get_recursive_partition_groups() %>%
                                       recursive_partition_pairwise() %>%
                                       jaccard_index(observed_partitions_complete_linking, .))) %>%
  unnest(everything())

# Complete linking (bootstrap)

jaccard_stability_bootstrap_complete_linking <-
  bootstrapped_groups_complete_linking %>%
  transmute(boot_id,
            jaccard_similarity = map(.x = recursive_partition,
                                     .f = ~recursive_partition_pairwise(.x) %>%
                                       jaccard_index(observed_partitions_complete_linking, .))) %>%
  unnest(everything())

# Complete sorting (pairwise)

observed_partitions_complete_sorting <-
  recursive_partition_groups_complete_sorting %>%
  recursive_partition_pairwise()

jaccard_stability_pairwise_complete_sorting <-
  pairwise_resample_complete_sorting %>%
  transmute(boot_id,
            dissimilarites = map(.x = simulated_graph,
                                 .f = ~get_dissimilarity_from_graph(.x))) %>%
  unnest(everything()) %>%
  nest(data = -boot_id) %>%
  transmute(boot_id,
            jaccard_similarity = map(.x = data,
                                     .f = ~get_dissimilarity_matrix(.x) %>%
                                       get_recursive_partition_groups() %>%
                                       recursive_partition_pairwise() %>%
                                       jaccard_index(observed_partitions_complete_sorting, .))) %>%
  unnest(everything())

# Complete sorting (bootstrap)

jaccard_stability_bootstrap_complete_sorting <-
  bootstrapped_groups_complete_sorting %>%
  transmute(boot_id,
            jaccard_similarity = map(.x = recursive_partition,
                                     .f = ~recursive_partition_pairwise(.x) %>%
                                       jaccard_index(observed_partitions_complete_sorting, .))) %>%
  unnest(everything())

# Combine everything

jaccard_stability_statistics <-
  jaccard_stability_incomplete_linking %>%
  group_by(sample) %>%
  summarize(stability = mean(jaccard_index)) %>%
  left_join(recursive_partition_groups_incomplete_linking) %>%
  mutate(sample = as.numeric(sample)) %>%
  left_join(sample_ids) %>%
  transmute(study = "incomplete linking",
            type = "pairwise",
            sample = name,
            stability,
            group) %>%
  bind_rows(
    jaccard_stability_bootstrap_complete_linking %>%
      group_by(sample) %>%
      summarize(stability = mean(jaccard_index)) %>%
      transmute(study = "complete linking",
                type = "bootstrap",
                sample,
                stability) %>%
      left_join(recursive_partition_groups_complete_linking),
    jaccard_stability_bootstrap_complete_sorting %>%
      group_by(sample) %>%
      summarize(stability = mean(jaccard_index)) %>%
      transmute(study = "complete sorting",
                type = "bootstrap",
                sample,
                stability) %>%
      left_join(recursive_partition_groups_complete_sorting),
    jaccard_stability_pairwise_complete_linking %>%
      group_by(sample) %>%
      summarize(stability = mean(jaccard_index)) %>%
      transmute(study = "complete linking",
                type = "pairwise",
                sample,
                stability) %>%
      left_join(recursive_partition_groups_complete_linking),
    jaccard_stability_pairwise_complete_sorting %>%
      group_by(sample) %>%
      summarize(stability = mean(jaccard_index)) %>%
      transmute(study = "complete sorting",
                type = "pairwise",
                sample,
                stability) %>%
      left_join(recursive_partition_groups_complete_sorting)) %>%
  mutate(sample = str_replace_all(sample, "_|\\?", " ") %>%
           str_replace_all("  ", " ") %>%
           str_squish())

# Jaccard stability for individual samples

library(tidytext)

p_jaccard_individual <-
  jaccard_stability_statistics %>%
  unite(study, type, group, col = "study_type", remove = FALSE) %>%
  mutate(sample = factor(sample) %>%
           reorder_within(group, within = study)) %>%
  ggplot(aes(x = sample, y = stability, color = type)) +
  geom_line(aes(group = study_type), size = 1) +
  geom_point(shape = 21, size = 3, fill = "white") +
  facet_wrap(~ study, scales = "free_x") +
  scale_x_reordered() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "top") +
  labs(x = NULL,
       y = "Jaccard stability",
       color = "resampling strategy") +
  scale_color_manual(values = c("darkgreen", "tan")) +
  lims(y = c(0, 1))

ggsave(filename = "img/study_1_jaccard_individual.png",
       plot = p_jaccard_individual,
       units = "in", height = 4, width = 6, scale = 2, dpi = 300)


# Quick and dirty comparison ----------------------------------------------

# Below is GPA merely to check that the configurations from the 3 methods are
# indeed similar, except for the placement of MES48

gpa_res <-
  diss_matrix_incomplete_linking %>%
  cmdscale() %>%
  as_tibble(rownames = "sample") %>%
  left_join(sample_ids %>% mutate(sample = as.character(sample))) %>%
  transmute(name = str_remove_all(name, " "),
            x = V1,
            y = V2) %>%
  left_join(
    diss_matrix_complete_linking %>%
      cmdscale() %>%
      as_tibble(rownames = "name") %>%
      mutate(name = str_remove_all(name, " "))
  ) %>%
  left_join(
    diss_matrix_complete_sorting %>%
      cmdscale() %>%
      as_tibble(rownames = "name") %>%
      mutate(name = str_remove_all(name, " ")),
    by = "name"
  ) %>%
  column_to_rownames("name") %>%
  FactoMineR::GPA(group = c(2, 2, 2))

gpa_res$RV
plot(gpa_res)












