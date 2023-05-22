# This script imports the resampled data from the Cocoa Wheel of Excellence
# study (study 2 in the manuscript) and analyzes it and the original data for
# grouping, stability, etc.

# Load resources ----------------------------------------------------------

library(here)
library(gt)
library(flextable)
source(here("code/helper_functions.R"))
library(tidytext)

clean_incomplete_linking <-
  read_csv(here("data/study 2/clean_incomplete_linking_table.csv"))

clean_incomplete_sorting <-
  read_csv(here("data/study 2/clean_incomplete_sorting_table.csv"))

clean_complete_sorting <-
  read_csv(here("data/study 2/clean_complete_sorting_table.csv"))

ibd_design <-
  read_csv(here("data/study 2/ibd_design.csv")) %>%
  as.matrix()

load(here("outputs/study 2/chocolate wheel resampled data.RData"))

chocolate_wheel <-
  read_csv(here("data/study 2/Chocolate Wheel Samples.csv")) %>%
  mutate_all(~tolower(.))


# Pull 24 sample sets from complete sorting -------------------------------

set.seed(1234) # for reproducibility

clean_complete_sorting <-
  clean_complete_sorting %>%
  nest(data = -sample_set) %>%
  slice_sample(n = 24) %>%
  unnest(everything())

# Grouping analysis by recursive partitioning -----------------------------

# Incomplete linking

diss_matrix_incomplete_linking <-
  clean_incomplete_linking %>%
  nest(data = -sample_set) %>%
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
            by = c("name" = "sample"))

# Incomplete sorting

diss_matrix_incomplete_sorting <-
  clean_incomplete_sorting %>%
  nest(data = -sample_set) %>%
  transmute(dissimilarity = map(.x = data,
                                .f = ~get_dissimilarity_table(.x))) %>%
  unnest(dissimilarity) %>%
  get_dissimilarity_matrix()

recursive_partition_groups_incomplete_sorting <-
  diss_matrix_incomplete_sorting %>%
  get_recursive_partition_groups()

additive_tree_incomplete_sorting <-
  diss_matrix_incomplete_sorting %>%
  dist(method = "max") %>%
  nj() %>%
  as_tbl_graph() %>%
  left_join(recursive_partition_groups_incomplete_sorting,
            by = c("name" = "sample"))

# Complete sorting

diss_matrix_complete_sorting <-
  clean_complete_sorting %>%
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

# Plots

library(paletteer)

this_palette <- paletteer_d("pals::alphabet")

# Additive trees

tci_incomplete_link <-
  diss_matrix_incomplete_linking %>%
  dist(method = "max") %>%
  nj() %>%
  TreeTools::TotalCopheneticIndex()

p_tree_incomplete_link <-
  additive_tree_incomplete_linking %>%
  mutate(type = ifelse(str_detect(name, "Node"), "internal", "leaf"),
         group = as.factor(group)) %>%
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
  expand_limits(x = c(-5, 5), y = c(-5, 5)) +
  scale_color_manual(values = this_palette, aesthetics = c("color", "edge_color")) +
  labs(caption = "incomplete linking") +
  annotate(geom = "text", x = 0, y = -5, hjust = 0,
           label = bquote(Total~Cophenetic~Index==.(tci_incomplete_link)))

tci_incomplete_sort <-
  diss_matrix_incomplete_sorting %>%
  dist(method = "max") %>%
  nj() %>%
  TreeTools::TotalCopheneticIndex()

p_tree_incomplete_sort <-
  additive_tree_incomplete_sorting %>%
  mutate(type = ifelse(str_detect(name, "Node"), "internal", "leaf"),
         group = as.factor(group)) %>%
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
  expand_limits(x = c(-10, 10), y = c(-10, 10)) +
  scale_color_manual(values = this_palette, aesthetics = c("color", "edge_color")) +
  labs(caption = "incomplete sorting") +
  annotate(geom = "text", x = 0.5, y = 10, hjust = 0,
           label = bquote(Total~Cophenetic~Index==.(tci_incomplete_sort)))

tci_complete_sort <-
  diss_matrix_complete_sorting %>%
  dist(method = "max") %>%
  nj() %>%
  TreeTools::TotalCopheneticIndex()

p_tree_complete_sort <-
  additive_tree_complete_sorting %>%
  mutate(type = ifelse(str_detect(name, "Node"), "internal", "leaf"),
         group = as.factor(group)) %>%
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
  expand_limits(x = c(-20, 20), y = c(-20, 20)) +
  scale_color_manual(values = this_palette, aesthetics = c("color", "edge_color")) +
  labs(caption = "complete sorting") +
  annotate(geom = "text", x = 1, y = 17, hjust = 0,
           label = bquote(Total~Cophenetic~Index==.(tci_complete_sort)))

ggsave(filename = "img/study_2_tree_incomplete_linking.png",
       plot = p_tree_incomplete_link,
       width = 3, height = 3, units = "in", scale = 3, dpi = 300)
ggsave(filename = "img/study_2_tree_incomplete_sorting.png",
       plot = p_tree_incomplete_sort,
       width = 3, height = 3, units = "in", scale = 3, dpi = 300)
ggsave(filename = "img/study_2_tree_complete_sorting.png",
       plot = p_tree_complete_sort,
       width = 3, height = 3, units = "in", scale = 3, dpi = 300)

# Graph statistics --------------------------------------------------------

graph_statistics <-
  tibble(study = c("incomplete linking", "incomplete sorting", "complete sorting"),
         data = list(clean_incomplete_linking,
                     clean_incomplete_sorting,
                     clean_complete_sorting)) %>%
  unnest(data) %>%
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
  save_as_docx(path = "tables/study_2_graph_stats.docx")

# Jaccard stability (resampling analysis) ---------------------------------

# Unfortunately this is quite slow; it will take ~30 minutes to run all of this.
# Therefore, if possible do not run the commented code, simply load the saved
# Jaccard stability files.
#
# # Incomplete linking
#
# observed_partitions_incomplete_linking <-
#   recursive_partition_groups_incomplete_linking %>%
#   recursive_partition_pairwise()
#
# tic()
# jaccard_stability_incomplete_linking <-
#   incomplete_linking_resampled_pairwise %>%
#   transmute(boot_id,
#             dissimilarities = map(.x = simulated_graph,
#                                   .f = ~get_dissimilarity_from_graph(.x))) %>%
#   unnest(everything()) %>%
#   nest(data = -boot_id) %>%
#   transmute(boot_id,
#             jaccard_similarity = map(.x = data,
#                                      .f = ~get_dissimilarity_matrix(.x) %>%
#                                        get_recursive_partition_groups() %>%
#                                        recursive_partition_pairwise() %>%
#                                        jaccard_index(observed_partitions_incomplete_linking, .))) %>%
#   unnest(everything())
# toc()
#
# # Incomplete sorting
#
# observed_partitions_incomplete_sorting <-
#   recursive_partition_groups_incomplete_sorting %>%
#   recursive_partition_pairwise()
#
# tic()
# jaccard_stability_incomplete_sorting <-
#   incomplete_sorting_resampled_pairwise %>%
#   transmute(boot_id,
#             dissimilarities = map(.x = simulated_graph,
#                                   .f = ~get_dissimilarity_from_graph(.x))) %>%
#   unnest(everything()) %>%
#   nest(data = -boot_id) %>%
#   transmute(boot_id,
#             jaccard_similarity = map(.x = data,
#                                      .f = ~get_dissimilarity_matrix(.x) %>%
#                                        get_recursive_partition_groups() %>%
#                                        recursive_partition_pairwise() %>%
#                                        jaccard_index(observed_partitions_incomplete_sorting, .))) %>%
#   unnest(everything())
# toc()
#
#
# # Complete sorting - pairwise
#
# observed_partitions_complete_sorting <-
#   recursive_partition_groups_complete_sorting %>%
#   recursive_partition_pairwise()
#
# tic()
# jaccard_stability_complete_sorting_pairwise <-
#   complete_sorting_resampled_pairwise %>%
#   transmute(boot_id,
#             dissimilarities = map(.x = simulated_graph,
#                                   .f = ~get_dissimilarity_from_graph(.x))) %>%
#   unnest(everything()) %>%
#   nest(data = -boot_id) %>%
#   transmute(boot_id,
#             jaccard_similarity = map(.x = data,
#                                      .f = ~get_dissimilarity_matrix(.x) %>%
#                                        get_recursive_partition_groups() %>%
#                                        recursive_partition_pairwise() %>%
#                                        jaccard_index(observed_partitions_complete_sorting, .))) %>%
#   unnest(everything())
# toc()
#
# # Complete sorting - classic
#
# tic()
# jaccard_stability_complete_sorting_classic <-
#   complete_sorting_resampled_classic %>%
#   transmute(boot_id,
#             jaccard_similarity = map(.x = sorting_rc_groups,
#                                      .f = ~recursive_partition_pairwise(.x) %>%
#                                        jaccard_index(observed_partitions_complete_sorting, .))) %>%
#   unnest(everything())
# toc()
#
# save(list = c("jaccard_stability_incomplete_linking",
#               "jaccard_stability_incomplete_sorting",
#               "jaccard_stability_complete_sorting_pairwise",
#               "jaccard_stability_complete_sorting_classic"),
#      file = here("outputs/study 2/chocolate wheel jaccard stability.RData"))

# Start here if the file "chocolate wheel jaccard stability.RData" is available.

load(here("outputs/study 2/chocolate wheel jaccard stability.RData"))

# Get a quantile estimate for the complete sorting

jaccard_stability_statistics <-
  # Start with the results for the complete sort, classic bootstrap
  jaccard_stability_complete_sorting_classic %>%
  group_by(sample) %>%
  ggdist::mean_qi(jaccard_index) %>%
  mutate(study = "complete sorting",
         type = "classic bootstrap") %>%
  left_join(recursive_partition_groups_complete_sorting) %>%
  # Add the rest of the studies
  bind_rows(
    # The pairwise complete sort
    jaccard_stability_complete_sorting_pairwise %>%
      group_by(sample) %>%
      ggdist::mean_qi(jaccard_index) %>%
      mutate(study = "complete sorting",
             type = "pairwise resample") %>%
      left_join(recursive_partition_groups_complete_sorting),
    jaccard_stability_incomplete_linking %>%
      group_by(sample) %>%
      ggdist::mean_qi(jaccard_index) %>%
      mutate(study = "incomplete linking",
             type = "pairwise resample") %>%
      left_join(recursive_partition_groups_incomplete_linking),
    jaccard_stability_incomplete_sorting %>%
      group_by(sample) %>%
      ggdist::mean_qi(jaccard_index) %>%
      mutate(study = "incomplete sorting",
             type = "pairwise resample") %>%
      left_join(recursive_partition_groups_incomplete_sorting)
  )

# Individual sample stability, line plot

horizontal_line_summary_data <-
  jaccard_stability_statistics %>%
  group_by(study, type) %>%
  summarize(mean = mean(jaccard_index))

p_stability_plot_individual <-
  jaccard_stability_statistics %>%
  group_by(sample) %>%
  mutate(average_stability = mean(jaccard_index)) %>%
  ungroup() %>%
  mutate(sample = factor(sample) %>% fct_reorder(-average_stability)) %>%
  unite(study, type, col = "study_type", remove = FALSE) %>%
  ggplot(aes(x = sample, y = jaccard_index, group = study_type)) +
  geom_hline(aes(yintercept = mean, color = study, linetype = type),
             size = 0.2,
             data = horizontal_line_summary_data) +
  geom_line(aes(color = study, linetype = type), size = 0.75, alpha = 0.8) +
  geom_point(aes(color = study), shape = 21, fill = "white", size = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "bottom",
        plot.subtitle = element_text(face = "bold")) +
  lims(y = c(0, 1)) +
  labs(x = NULL,
       y = "Jaccard stability",
       linetype = "resampling strategy",
       subtitle = "Individual sample stability") +
  scale_color_manual(values = this_palette[c(6, 11, 12)])

# Groupwise sample stability, line plot

p_stability_plot_group <-
  jaccard_stability_statistics %>%
  group_by(study, type, group) %>%
  summarize(group_stability = mean(jaccard_index)) %>%
  ungroup() %>%
  arrange(study, type, -group_stability) %>%
  group_by(study, type) %>%
  mutate(group_name = row_number() %>% as.factor()) %>%
  ungroup() %>%
  unite(study, type, col = "study_type", remove = FALSE) %>%
  ggplot(aes(x = group_name, y = group_stability)) +
  geom_hline(aes(yintercept = mean, color = study, linetype = type),
             size = 0.2,
             data = horizontal_line_summary_data, show.legend = FALSE) +
  geom_line(aes(group = study_type, color = study, linetype = type),
            size = 0.75, alpha = 0.8, show.legend = FALSE) +
  geom_point(aes(color = study), shape = 21, fill = "white", size = 3,
             show.legend = FALSE) +
  theme_classic() +
  theme(plot.subtitle = element_text(face = "bold")) +
  lims(y = c(0, 1)) +
  labs(x = NULL,
       y = "Jaccard stability",
       linetype = "resampling strategy",
       subtitle = "Group stability") +
  scale_color_manual(values = this_palette[c(6, 11, 12)])

p_stability <-
  p_stability_plot_individual / p_stability_plot_group

ggsave(filename = "img/study_2_jaccard_stability.png",
       plot = p_stability,
       height = 3, width = 5, units = "in", scale = 3, dpi = 300)

# Looking at the groups ---------------------------------------------------

# OK I am going to attempt to write a function to find "most similar" by Jaccard
# index groups

# This function takes two data frames that must each contain a $group column and
# a $members list-column that is a character vector.  The input_data$members
# column is matched via Jaccard coefficient against the target_data$members
# columns.  The function returns a new data frame containing all data from the
# input_data, as well as $target_group indexing the matched data, and a
# $jaccard_index column listing the match coefficient.

jaccard_alignment <- function(input_data, target_data){

  output_data <-
    input_data %>%
    transmute(group,
              members,
              target_group = 0,
              jaccard_index = 0)

  for(i in 1:nrow(input_data)){

    members_i <-
      input_data[i, 2]$members %>%
      unlist()

    for(j in 1:nrow(target_data)){

      members_j <-
        target_data[j, 2]$members %>%
        unlist()

      group_j <-
        target_data[j, 1]$group

      jaccard_match <-
        length(intersect(members_i, members_j)) / length(union(members_i, members_j))

      if(jaccard_match > output_data[i, ]$jaccard_index){
        output_data[i, 3] <- group_j
        output_data[i, 4] <- jaccard_match
      }
    }
  }

  return(output_data)

}

# Now we set up the proper data frames for alignment

recursive_partition_vectors_complete_sorting <-
  recursive_partition_groups_complete_sorting %>%
  group_by(group) %>%
  summarize(members = str_c(sample, collapse = ",")) %>%
  mutate(members = map(members, ~str_split(.x, ",", simplify = TRUE) %>% as.vector()))

recursive_partition_vectors_incomplete_sorting <-
  recursive_partition_groups_incomplete_sorting %>%
  group_by(group) %>%
  summarize(members = str_c(sample, collapse = ",")) %>%
  mutate(members = map(members, ~str_split(.x, ",", simplify = TRUE) %>% as.vector()))

recursive_partition_vectors_incomplete_linking <-
  recursive_partition_groups_incomplete_linking %>%
  group_by(group) %>%
  summarize(members = str_c(sample, collapse = ",")) %>%
  mutate(members = map(members, ~str_split(.x, ",", simplify = TRUE) %>% as.vector()))

partition_vectors_chocolate_wheel <-
  chocolate_wheel %>%
  transmute(group = as.factor(`Cocoa of excellence group`) %>% as.numeric(),
            sample = `Sample name`) %>%
  group_by(group) %>%
  summarize(members = str_c(sample, collapse = ",")) %>%
  mutate(members = map(members, ~str_split(.x, ",", simplify = TRUE) %>% as.vector()))


# We can align all the derived groups with the original wheel groups.  NB: the
# magnitude of the Jaccard Index is biased downward when there are more groups
# in one set.

aligned_groups_wheel <-
  bind_rows(
    jaccard_alignment(recursive_partition_vectors_complete_sorting, partition_vectors_chocolate_wheel) %>%
      mutate(members = map_chr(members, ~str_c(.x, collapse = ", ")),
             study = "complete sorting"),
    jaccard_alignment(recursive_partition_vectors_incomplete_linking, partition_vectors_chocolate_wheel) %>%
      mutate(members = map_chr(members, ~str_c(.x, collapse = ", ")),
             study = "incomplete linking"),
    jaccard_alignment(recursive_partition_vectors_incomplete_sorting, partition_vectors_chocolate_wheel) %>%
      mutate(members = map_chr(members, ~str_c(.x, collapse = ", ")),
             study = "incomplete sorting")
  ) %>%
  left_join(partition_vectors_chocolate_wheel, by = c("target_group" = "group")) %>%
  transmute(study,
            members = members.x,
            wheel_group = map_chr(members.y, ~str_c(.x, collapse = ", ")),
            jaccard_index)

aligned_groups_wheel %>%
  arrange(study, wheel_group)

# This tabular output probably needs to get written to CSV to format in
# something like Word.  It can be edited via gt but there need to be a number of
# cell merges, etc, which I am not sure whether gt can do.


# First, we output a readable version of the chocolate wheel for reproduction in
# our paper.

ft_aligned_groups <-
  chocolate_wheel %>%
  arrange(`Cocoa of excellence group`) %>%
  group_by(`Cocoa of excellence group`) %>%
  summarize(wheel_group = str_c(`Sample name`, collapse = ", ")) %>%
  left_join(
    aligned_groups_wheel %>%
      select(-jaccard_index) %>%
      pivot_wider(names_from = study, values_from = members)
  ) %>%
  unnest(everything()) %>%
  flextable()

ft_aligned_groups %>%
  colformat_char(na_str = "NA") %>%
  merge_v() %>%
  set_header_labels(`Cocoa of excellence group` = "Group Name",
                    wheel_group = "Group Members") %>%
  save_as_docx(path = "tables/study_2_group_alignment.docx")


# Quality check with GPA --------------------------------------------------

res_gpa <-
  diss_matrix_complete_sorting %>%
  cmdscale() %>%
  as_tibble(rownames = "sample") %>%
  left_join(
    diss_matrix_incomplete_linking %>%
      cmdscale() %>%
      as_tibble(rownames = "sample"),
    by = "sample") %>%
  left_join(
    diss_matrix_incomplete_sorting %>%
      cmdscale() %>%
      as_tibble(rownames = "sample"),
    by = "sample"
  ) %>%
  column_to_rownames("sample") %>%
  FactoMineR::GPA(group = c(2, 2, 2), name.group = c("complete sort", "incomplete link", "incomplete sort"))

res_gpa$RV
res_gpa$simi
res_gpa$dep
res_gpa$correlations
res_gpa$PANOVA


































