###############################################################################
#'
#' validation_simulation/vis_sims.R
#' Author: Ben R. Goldstein
#' Description: Visualize results from the validation simulation
#' 
###############################################################################

library(tidyverse)

result_files <- list.files('intermediate/validation_sim', full.names = TRUE)

summary_list <- list()
for (i in 1:length(result_files)) {
  thisres <- readRDS(result_files[[i]])
  
  summary_list[[i]] <- thisres$result
}

theming <- list(
  theme_bw(),
  theme(axis.ticks = element_blank()),
  ylab("% of simulated effects"),
  xlab(""),
  scale_fill_manual("Correct\ninference?",
                    values = c("TRUE" = "#444444", "FALSE" = "#dd4040"))
)


summary_all <- bind_rows(summary_list) %>% 
  filter(!grepl("logDens", param)) %>% 
  mutate(
    `50%` = (mean < 0.5 & truth == 0) | (mean > 0.5 & truth == 1),
    `90%` = (mean < 0.9 & truth == 0) | (mean > 0.9 & truth == 1),
    truth = ifelse(truth, "SVC", "Non-SVC")
  ) %>% 
  select(`50%`, `90%`, species, scenario, truth, parnum, type, iter) %>% 
  pivot_longer(cols = c("50%", "90%"))

### OVERALL FP and FN RATES
summary_all %>% 
  group_by(name, truth) %>% 
  summarize(false_rate = 1 - mean(value))


# Which variables have multiple SVCs?
multiple_SVC_vars <- summary_all %>% 
  filter(name == "50%") %>% 
  group_by(species, scenario, parnum, iter) %>% 
  summarize(num_SVC = sum(truth == "SVC"), result = mean(value)) %>% 
  filter(num_SVC > 1) %>% 
  mutate(id = paste0(species, scenario, parnum, iter))

# Get overall FP and FN for these
summary_all %>% 
  filter(paste0(species, scenario, parnum, iter) %in% multiple_SVC_vars$id) %>% 
  group_by(name, truth) %>% 
  summarize(false_rate = 1 - mean(value))
# ...and broken out by svc type
summary_all %>% 
  filter(paste0(species, scenario, parnum, iter) %in% multiple_SVC_vars$id) %>% 
  group_by(name, truth, type) %>% 
  summarize(false_rate = 1 - mean(value))

# Overall counts, by cutoff level
ggplot(summary_all) +
  geom_bar(aes(name, fill = value), position = "fill") +
  theme_bw() + theme(axis.ticks = element_blank()) + ylab("% of simulated effects") + xlab("") +
  scale_fill_manual("Correct inference?",
                    values = c("TRUE" = "#444444", "FALSE" = "#dd4040"))

# Use only 50% cutoff, display by truth
summary_all %>% 
  filter(name == "50%") %>% 
  ggplot() +
  geom_bar(aes(truth, fill = value), position = "fill") +
  theme_bw() + theme(axis.ticks = element_blank()) + ylab("% of simulated effects") + xlab("") +
  scale_fill_manual("Correct inference?",
                    values = c("TRUE" = "#444444", "FALSE" = "#dd4040"))


simulation_results_fig <- summary_all %>% 
  mutate(outcome = ifelse(value, "Correct", 
                          ifelse(truth == "SVC", "False negative", "False positive"))) %>%
  mutate(type = recode(type, "eco" = "Ecoregion", "lin" = "Lineage")) %>% 
  filter(name == "50%") %>% 
  ggplot() +
  geom_bar(aes(truth, fill = outcome), position = "stack") +
  facet_grid(type~paste0("Scenario ", scenario)) +
  scale_fill_manual("",
                    values = c("Correct" = "#e0e0e0", 
                               "False negative" = "#d1ac17",
                               "False positive" = "#dd4040"
                               )) +
  xlab("What was the true effect?") + ylab("Number of effects") +
  theme_bw() +
  theme(axis.ticks = element_blank(), panel.grid = element_blank(),
        strip.background =element_rect(fill="white"))

ggsave(filename = "paper_figs/simulation_results.jpg", 
       simulation_results_fig, width = 5, height = 3.5, dpi = 600)
ggsave(filename = "paper_figs/simulation_results.pdf", 
       simulation_results_fig, width = 5, height = 3.5, dpi = 600)


# Look at scenario 3 more closely
summary_all %>% 
  filter(name == "50%", scenario == 3) %>% 
  View()

