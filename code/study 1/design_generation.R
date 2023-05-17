library(tidyverse)
library(crossdes)
library(numbers)
library(here)


# Design generation -------------------------------------------------------

set.seed(1234)

# This should generate a crossover balanced design, with trt = 10 treatments (types of chocolate),
# each seen 20 times (adding a multiple of 3 helps with factorization)

# Number of possible blocks (without regard to order):
10 * 20 # = 200

# Now get the prime factors of 200
primeFactors(200) # = 2^3 * 5^2

generated_design <- find.BIB(trt = 10, b = 40, k = 6) # This is very slow if b > 10

# Check the generated design

isGYD(generated_design, tables = T, type = T) # This provides us with a design that is partially balanced in that every treatment occurs 24 times

save(generated_design, file = here("outputs", "study_1_generated_design.RData"))
