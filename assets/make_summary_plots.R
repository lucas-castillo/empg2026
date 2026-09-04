library(tidyverse)
library(jsonlite)
select <- dplyr::select
parse_file <- function(data_path, file_name){
  id_number <- file_name %>% 
    str_remove("experiment_results_") %>% 
    str_remove(".csv") %>% 
    as.numeric()
  
  
  df <- read_csv(paste(data_path, file_name, sep=""))
  pseudonym <- jsonlite::fromJSON(df[[2,"response"]])$Q0
  df %>% filter(trial_index > 3) %>% 
    mutate(response = sapply(response, \(x) jsonlite::fromJSON(x)$Q0)) %>% 
    mutate(response = as.numeric(response)) %>% 
    mutate(index = 1:n()) %>% 
    mutate(id_number, pseudonym) %>% 
    select(id_number, pseudonym, index, response, rt)
}


DF <- tibble()
for (f in list.files(data_path)){
  DF <- bind_rows(DF, parse_file(data_path, f) %>% mutate(rt = as.numeric(rt)))
}

# Clean data --------------------------------------------------------------
DF <- DF %>% 
  group_by(id_number, pseudonym) %>% 
  mutate(ub = mean(response, na.rm = T) + 4 * sd(response, na.rm=T)) %>%
  mutate(lb = mean(response, na.rm = T) - 4 * sd(response, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(implausible = response > ub | response < lb) %>% 
  mutate(response2 = ifelse(implausible, NA, response)) %>% 
  select(-c(ub, lb))



# Flag non-randoms --------------------------------------------------------
source("assets/RQA.R")
DET_df <- DF %>% 
  group_by(id_number, pseudonym) %>% 
  do(
    tibble(
      DET1 = RQA(.$response2)$DET,
      DET2 = RQA(diff(.$response2))$DET,
    )
  )
  
DET_guilty <- DET_df %>% 
  filter(DET1 > .85 | DET2 > .85) %>% 
  pull(id_number, pseudonym)

source("assets/get_1D_summaries.R")
obs_summaries <- DF %>% 
  group_by(id_number, pseudonym) %>% 
  do(
    all_measures(.$response2)
  ) %>% 
  ungroup

r_expectations <- DF %>% 
  group_by(id_number, pseudonym) %>% 
  do(
    random_expectations(.$response2)
  ) %>% 
  ungroup

r_expectations %>% 
  summarise(across(R:gzip, mean))

obs_summaries.p <- obs_summaries
r_expectations.p <- r_expectations

obs_summaries.p$S <- obs_summaries.p$S * 1000
r_expectations.p$S <- r_expectations.p$S * 1000

fancy_names <- c(
  "R" = "Repetitions",
  "A" = "Adjacencies",
  "TPF" = "Turning Points",
  "D" = "Distances",
  "S" = "Shape (x1000)",
  "RNG" = "RNG",
  "gzip" = "gzip"
)

obs_summaries.p <- obs_summaries.p %>% 
  pivot_longer(R:gzip) %>% 
  mutate(name = factor(fancy_names[name], levels=fancy_names)) %>% 
  mutate(x=rnorm(n(), sd=.1))

most_random <- bind_rows(
  obs_summaries %>% mutate(expected=F),
  r_expectations %>% mutate(expected=T)
) %>% 
  pivot_longer(R:gzip) %>% 
  group_by(name) %>% 
  mutate(M = mean(value, na.rm=T), SD = sd(value, na.rm=T)) %>% 
  mutate(value = (value - M) / SD) %>% 
  pivot_wider(
    names_from = expected, 
    names_prefix = "e_", 
    values_from = value, id_cols = c(pseudonym, name)
  ) %>% 
  mutate(ss = (e_FALSE - e_TRUE)**2) %>% 
  group_by(pseudonym) %>% 
  summarise(sss = sum(ss)) %>% 
  filter(sss == min(sss)) %>% 
  pull(pseudonym)


p0 <- obs_summaries.p %>% 
  mutate(mr = pseudonym == most_random) %>% 
  ggplot(aes(x, value)) + 
  geom_point(
    aes(alpha = name, fill=mr),
    shape=21, colour="black", size=6
  ) + 
  facet_wrap(vars(name), ncol=7, scales="free", strip.position="bottom") +
  geom_hline(
    data = r_expectations.p %>% 
      summarise(across(R:gzip, mean)) %>% 
      pivot_longer(R:gzip) %>% 
      mutate(name = factor(fancy_names[name], levels=fancy_names)),
    mapping=aes(yintercept = value)
  ) + 
  scale_x_continuous(breaks=NULL, minor_breaks = NULL, limits = c(-1, 1), expand=c(0,0)) +
  theme_minimal(22) + 
  theme(
    strip.clip = "off", 
    panel.background = element_rect(fill = "#eeeeee", color="transparent"), 
    plot.background = element_rect(fill = "#eeeeee"), 
    panel.grid = element_blank(), 
    legend.position = "none"
  ) +
  xlab("Measure") + ylab("Value") + 
  scale_alpha_manual(
    values = c(
      0, 0, 0, 0, 0, 0, 0
    )
  ) + guides(alpha="none") + 
  scale_fill_manual(values = c("#375C8C", "#375C8C"))

p1 <- p0 + scale_alpha_manual(
  values = c(
    1, 0, 0, 0, 0, 0, 0
  )
)

p2 <- p0 + scale_alpha_manual(
  values = c(
    1, 1, 0, 0, 0, 0, 0
  )
)
p3 <- p0 + scale_alpha_manual(
  values = c(
    1, 1, 1, 0, 0, 0, 0
  )
)
p4 <- p0 + scale_alpha_manual(
  values = c(
    1, 1, 1, 1, 0, 0, 0
  )
)
p5 <- p0 + scale_alpha_manual(
  values = c(
    1, 1, 1, 1, 1, 0, 0
  )
)
p6 <- p0 + scale_alpha_manual(
  values = c(
    1, 1, 1, 1, 1, 1, 0
  )
)
p7 <- p0 + scale_alpha_manual(
  values = c(
    1, 1, 1, 1, 1, 1, 1
  )
)
p8 <- p7 + scale_fill_manual(values = c("#375C8C", "#F2D541"))

saveRDS(object = obs_summaries, file = paste("files_for_attendees/", str_remove(data_path, "data/"), "obs_summaries.rds", sep=""))
saveRDS(object = DF, file = paste("files_for_attendees/", str_remove(data_path, "data/"), "obs_data.rds", sep=""))



