library(tidyverse)
library(crossdes)
library(numbers)
library(here)

set.seed(1234)

# Design generation -------------------------------------------------------

### DO NOT RUN -- IT WILL TAKE A WHILE ###
# 
# # This should generate a crossover balanced design, with trt = 62 treatments (words),
# # each seen 24 times (adding a multiple of 3 helps with factorization)
# 
# # Number of possible blocks (without regard to order):
# 62 * 24 # = 1488
# 
# # Now get the prime factors of 1488
# primeFactors(1488) # = 2^4 * 3 * 31, so we are going to choose 2^4 as the block size and 3*31 as the number of blocks
# 
# generated_design <- find.BIB(trt = 62, b = 93, k = 16) # This is very slow if b > 10
#
# save(generated_design, file = here("outputs", "study_2_generated_design.RData"))

# Design check ------------------------------------------------------------

load(here("outputs", "study_2_generated_design.RData"))

# Check the generated design

isGYD(generated_design, tables = T, type = T) # This provides us with a design that is partially balanced in that every treatment occurs 24 times
