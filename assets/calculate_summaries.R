library(tidyverse)
source("assets/get_1D_summaries.R")
simulations <- readRDS("files_for_attendees/simulations.rds")
parameters <- readRDS("files_for_attendees/parameters.rds")
sim_summaries <- apply(simulations, FUN = all_measures, MARGIN = 1)

sim_summaries <- data.frame(matrix(unlist(sim_summaries), nrow = nrow(simulations), byrow = T))
sim_summaries <- sim_summaries %>% 
  as_tibble() %>% 
  magrittr::set_colnames(c("R", "A", "TPF", "D", "S", "RNG", "gzip")) %>% 
  bind_cols(parameters, .)

saveRDS(sim_summaries, "files_for_attendees/sim_summaries.rds")
