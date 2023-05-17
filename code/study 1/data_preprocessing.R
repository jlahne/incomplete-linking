# This script imports and cleans the individual data sheets from the incomplete
# linking task and gives a single table that lists samples in each set along
# with whether or not a link was drawn.

library(tidyverse)
library(igraph)
library(tidygraph)
library(readxl)
library(here)



# Make a list of lower-triangular adjacency matrices 

adjacency_matrices <- tibble(id = character(),
                             adjacency_matrix = list())

# Make a table of links

link_table <- tibble(id = character(),
                     from = character(),
                     to = character(),
                     presence = double())


for(filename in dir("data/study 1/chocolate sheets/")){
  
  id <- str_extract(filename, "\\(.*\\)") %>%
    str_remove_all("\\(|\\)")
  
  raw_data <- read_excel(file.path("data/study 1", "chocolate sheets", filename))
  
  # Reduce the design to the actual incomplete block
  
  excluded_samples <- 
    raw_data %>%
    column_to_rownames("Sample") %>%
    is.na() %>%
    colSums() %>%
    as_tibble(rownames = "sample") %>%
    left_join(
      raw_data %>%
        column_to_rownames("Sample") %>%
        is.na() %>%
        rowSums() %>%
        as_tibble(rownames = "sample"),
      by = "sample"
    ) %>%
    # If the sample isn't included in the block, will have k - 1 = 9 NAs in both
    # row and column
    filter(value.x == 9 & value.y == 9) %>%
    pull(sample)
  
  adjacency_matrix <-   
    raw_data %>%
    filter(!(Sample %in% excluded_samples)) %>%
    select(!one_of(excluded_samples)) %>%
    column_to_rownames("Sample") %>%
    as.matrix()
  
  edge_list <- 
    adjacency_matrix %>%
    as_tibble(rownames = "sample") %>%
    pivot_longer(-sample) %>%
    rename(from = 1, to = 2, presence = 3) %>%
    filter(from != to) %>%
    drop_na(presence) %>%
    mutate(id = id) %>%
    relocate(id)
  
  adjacency_matrices <- 
    bind_rows(adjacency_matrices, 
              tibble(id = id, adjacency_matrix = list(adjacency_matrix)))
  
  link_table <- 
    bind_rows(link_table,
              edge_list)
  
}

write_csv(link_table, "data/study 1/study_1_edgelist_incomplete_linking.csv")

# Get the previous data wrangled into edgelists

load("data/study 1/previous_study_data.rds")

# complete sorting

complete_sorting_links <- tibble(sample_set = character(), 
                                 from = character(), 
                                 to = character(), 
                                 presence = numeric())

for(i in 1:dim(sorting_distance_array)[3]){
  complete_sorting_links <- 
    sorting_distance_array[, , i] %>%
    apply(c(1, 2), function(x) 1 - x) %>%
    as_tibble(rownames = "from") %>%
    pivot_longer(-from, names_to = "to", values_to = "presence") %>%
    mutate(sample_set = paste0("J", i)) %>%
    relocate(sample_set) %>%
    bind_rows(complete_sorting_links, .)
}

write_csv(complete_sorting_links, "data/study 1/study_1_edgelist_complete_sorting.csv")

# complete linking

complete_linking_links <- tibble(sample_set = character(), 
                                 from = character(), 
                                 to = character(), 
                                 presence = numeric())


for(i in 1:dim(linking_distance_array)[3]){
  complete_linking_links <- 
    linking_distance_array[, , i] %>%
    apply(c(1, 2), function(x) 1 - x) %>%
    as_tibble(rownames = "from") %>%
    pivot_longer(-from, names_to = "to", values_to = "presence") %>%
    mutate(sample_set = paste0("J", i)) %>%
    relocate(sample_set) %>%
    bind_rows(complete_linking_links, .)
  
}

write_csv(complete_linking_links, "data/study 1/study_1_edgelist_complete_linking.csv")

# sample table

chocolate_ids <- 
  read_excel("data/study 1/previous study Chocolate Study Codes.xlsx") %>%
  drop_na()

chocolate_ids %>%
  rename(sample = 1, brand = 2, type = 3, cocoa = 4, name = 5) %>%
  write_csv("data/study 1/study_1_samples.csv")



