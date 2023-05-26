# This is a set of common helper functions for incomplete linking analysis

library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(DistatisR)
library(ape)
library(readxl)
library(AddDistTreeSplit)
library(patchwork)
library(tictoc)

# Helper functions --------------------------------------------------------

# From a graph geodesic distance between two samples in [0, Inf), return a
# dissimilarity in [0, 1]
graph_dissimilarity <- function(geodesic_distance){
  1 - ifelse(is.infinite(geodesic_distance),
             0,
             ifelse(geodesic_distance == 0,
                    1,
                    1 / geodesic_distance)
  )
}

# From a table of dissimilarities (derived from graph/geodesic distances on a
# set of linking results), return a K x K matrix of dissimilarities
get_dissimilarity_matrix <- function(dissimilarity_table){

  dissimilarity_table %>%
    select(from, to, dissimilarity) %>%
    group_by(from, to) %>%
    summarize(total_dissimilarity = sum(dissimilarity)) %>%
    ungroup() %>%
    pivot_wider(names_from = to, values_from = total_dissimilarity) %>%
    column_to_rownames("from") %>%
    as.matrix()

}

# From a K x K dissimilarity matrix (like from get_dissimilarity_matrix()
# above), return a partition based on the Koenig et al 2021 recursive
# partitioning algorithm for trees
get_recursive_partition_groups <- function(dissimilarity_matrix){

  partitions <-
    dissimilarity_matrix %>%
    dist(method = "maximum") %>%
    recursive_partitioning()

  partitions[[length(partitions)]] %>%
    enframe() %>%
    unnest(value) %>%
    rename(group = 1, sample = 2) %>%
    arrange(sample)

}

# From a K x 2 partition table with .$group as the column of partitions and
# .$sample as the column of samples/objects that were partitioned, return a
# binary K x K matrix of similarities based on group/partition membership.
# (NB: This function is rewritten from DistatisR)
recursive_partition_pairwise <- function(partition){
  groups <- partition$group
  labels <- partition$sample
  similarity <- matrix(as.matrix(groups), nrow = length(groups), ncol = length(groups))
  similarity <- (similarity == t(similarity)) - 0
  dimnames(similarity) <- list(labels, labels)

  similarity

}

get_dissimilarity_table <- function(cooccurrence_table){

  cooccurrence_table %>%
    select(from, to, presence) %>%
    filter(from != to) %>%
    graph_from_data_frame(directed = FALSE) %>%
    as_tbl_graph() %>%
    activate(edges) %>%
    filter(presence != 0) %>%
    distances() %>%
    as_tibble(rownames = "from") %>%
    pivot_longer(cols = -from, names_to = "to", values_to = "geodesic_distance") %>%
    mutate(dissimilarity = graph_dissimilarity(geodesic_distance))

}

# From a simple, undirected, unweighted graph, get the graph dissimilarities in
# a table with .$from, .$to, .$geodesic_distance, and .$dissimilarity
get_dissimilarity_from_graph <- function(graph){

  graph %>%
    distances(weights = NA) %>%
    as_tibble(rownames = "from") %>%
    pivot_longer(cols = -from, names_to = "to", values_to = "geodesic_distance") %>%
    mutate(dissimilarity = graph_dissimilarity(geodesic_distance))

}

# Given a full data set from a sorting or linking study, with .$from, .$to, and
# .$presence columns, draw a specified subset of
# samples from that dataset, and sample the observed links between the samples
# in order to draw an accurately simulated resample without propagating NAs.
# Optionally, if transitivity is required for the samples, "as_sorting = TRUE"
# will force all groups to be cliques/completely transitive.
simulate_graph <- function(original_results, selected_samples, as_sorting = FALSE){

  # Here we simulate the graph based on the original results
  simulated_graph <-
    original_results %>%
    select(from, to, presence) %>%
    filter(from != to) %>%
    graph_from_data_frame(directed = FALSE) %>%
    induced_subgraph(vids = selected_samples) %>%
    as_data_frame() %>%
    group_by(from, to) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    graph_from_data_frame(directed = FALSE) %>%
    as_tbl_graph() %>%
    activate(edges) %>%
    filter(presence == 1)

  # If we are assuming that this is "linking like" we can allow non-transitive similarity
  if(!as_sorting){
    return(simulated_graph)
  }
  # Otherwise we will wrangle the components of the graph into fully connected cliques
  if(as_sorting){
    simulated_graph %>%
      components() %>%
      .$membership %>%
      as_tibble(rownames = "sample") %>%
      rename(sample = 1, group = 2) %>%
      recursive_partition_pairwise() %>%
      graph_from_adjacency_matrix(mode = "undirected", diag = FALSE) %>%
      as_tbl_graph()
  }
}


# A little helper to get the number of components in a graph, so that we can see
# how connected or disconnected an individual sort is.
component_count <- function(graph) graph %>% components() %>% .$no


# This function is a copy of simulate_graph() that also returns the presence ==
# 0 links purely for use in visualization, and so isn't to be called in general.
simulate_graph2 <- function(original_results, selected_samples, as_sorting = FALSE){

  # Here we simulate the graph based on the original results
  simulated_graph <-
    original_results %>%
    select(from, to, presence) %>%
    filter(from != to) %>%
    graph_from_data_frame(directed = FALSE) %>%
    induced_subgraph(vids = selected_samples) %>%
    as_data_frame() %>%
    group_by(from, to) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    graph_from_data_frame(directed = FALSE) %>%
    as_tbl_graph() %>%
    activate(edges)

  # If we are assuming that this is "linking like" we can allow non-transitive similarity
  if(!as_sorting){
    return(simulated_graph)
  }
  # Otherwise we will wrangle the components of the graph into fully connected cliques
  if(as_sorting){
    simulated_graph %>%
      filter(presence == 1) %>%
      components() %>%
      .$membership %>%
      as_tibble(rownames = "sample") %>%
      rename(sample = 1, group = 2) %>%
      recursive_partition_pairwise() %>%
      graph_from_adjacency_matrix(mode = "undirected", diag = FALSE) %>%
      as_tbl_graph() %>%
      activate(edges) %>%
      mutate(presence = 1)
  }
}

# Frequently, we want to take a table of links with .$from, .$to, and .$presence
# %in% c(0, 1) and return a graph that contains all observed links where
# .$presence == 1 AND any isolates.  Simply filtering for .$presence == 1 and
# then converting to a graph will drop those isolated nodes, so this function
# does the work for us.

# A function to generate a graph from an edgelist with a .$presence (in c(0, 1))
# variable, which indicates whether a link is actually drawn or not.

graph_from_links <-
  function(links){
    links %>%
      select(from, to, presence) %>%
      graph_from_data_frame(directed = FALSE) %>%
      simplify(edge.attr.comb = "first") %>%
      as_tbl_graph() %>%
      activate(edges) %>%
      filter(presence == 1)
  }

# This function takes two (indicator) matrices of the same dimensions that
# represent a (pairwise) partition of items, represented by the rows/columns,
# and calculates the jaccard index for each row (sample) in the the two
# matrices, which will range from 0 (the sample has no group members in common
# between the two partitions) and 1 (perfect agreement).

jaccard_index <- function(m1, m2){

  res <- numeric(nrow(m1))

  for(i in 1:nrow(m1)){

    j1 <- names(m1[i, ])[m1[i, ] == 1]
    j2 <- names(m2[i, ])[m2[i, ] == 1]

    res[i] <- length(intersect(j1, j2)) / length(union(j1, j2))

  }

  tibble(sample = row.names(m1),
         jaccard_index = res)

}
