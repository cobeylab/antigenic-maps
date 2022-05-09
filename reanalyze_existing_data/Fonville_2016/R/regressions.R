source('analysis_funs.R')


ferret_titers = import_ferret_titers()
ferret_distances = import_map_distances()
child_titers = import_human_titer_data()

individual_distances <- import_human_titer_data() %>%
  group_by(subject, season, age_mos) %>%
  extract(season, into = 'season_y1', regex = '(\\d\\d\\d\\d)/\\d+', convert = T, remove = F) %>% # Extract first year of season in which sample was collected
  mutate(earliest_possible_YOB = season_y1 - ceiling(age_mos/12),
         latest_possible_YO_sample = season_y1 + 2,
         circulated_in_infection_window = strain_year >= earliest_possible_YOB & strain_year <= latest_possible_YO_sample,
         max_titer_strain = ifelse(titer == max(titer), strain, NA),
         putative_primary_strain = ifelse(strain %in% get_putative_primary_strains(strain, titer, circulated_in_infection_window), as.character(strain), NA),
         subject_label = get_subject_factor_by_YOB(subject, earliest_possible_YOB, season),
         xval = jitter(strain_year),
         s_max = max(logtiter[circulated_in_infection_window]),
         individual_distance = s_max - logtiter)


merged_distances = individual_distances %>%
  group_by(subject, putative_primary_strain) %>%
  filter(!is.na(putative_primary_strain)) %>%
  summarise() %>%
  distinct() %>%
  ungroup() %>% ## Get a df of subjects and primary starins
  merge(  ## Merge with the ferret distances from each primary strain to all other strains
    ferret_distances %>%
      rename(putative_primary_strain = strain1,
             strain = strain2,
             ferret_map_distance = distance),
    all.y = T
  ) %>%
  merge(individual_distances %>% select(-putative_primary_strain)) %>%
  mutate(strain = factor(strain, levels = get_strains_chronologically(strain, strain_year))) %>%
  ungroup() %>% group_by(subject) %>%
  mutate(primary_strain_id = as.numeric(factor(putative_primary_strain)),
         x_dodge = as.numeric(strain) + (primary_strain_id-mean(primary_strain_id))*.2) %>%
  ungroup() %>%
  group_by(putative_primary_strain) %>%
  mutate(primary_strain_year = extract_strain_year(putative_primary_strain)) %>%
  mutate(temporal_distance = abs(primary_strain_year - strain_year)) %>%
  mutate(ferret_overestimate = ferret_map_distance - individual_distance)

ferret_distance_fit = lm(ferret_overestimate~0+ferret_map_distance, data = merged_distances)
ferret_distance_fit_grouped = lm(ferret_overestimate~0+ferret_map_distance*is_high_responder, data = merged_distances)
temporal_distance_fit = lm(ferret_overestimate~0+temporal_distance, data = merged_distances)
temporal_distance_fit_grouped = lm(ferret_overestimate~0+temporal_distance*is_high_responder, data = merged_distances)
anova(ferret_distance_fit_grouped, ferret_distance_fit)
anova(temporal_distance_fit_grouped, temporal_distance_fit)


tibble(model = c('temporal_distance', 'temporal_distance_grouped', 'ferret_distance', 'ferret_distance_grouped'),
       AIC = c(AIC(temporal_distance_fit), AIC(temporal_distance_fit_grouped), AIC(ferret_distance_fit), AIC(ferret_distance_fit_grouped))) %>%
  mutate(del_AIC = AIC-min(AIC)) %>%
  arrange(del_AIC)

## The ferret distance model grouped by fit is the best, but the high/low responder groups are not generalizable.
## Instead, save the ferret_distance model, which is the best of the non-grouped models

write_rds(x = ferret_distance_fit, '../processed_data/primary_ferret_human_distance_lm.rds')

merged_distances %>%
  ggplot(aes(x = ferret_map_distance, y = ferret_overestimate)) +
  geom_point(alpha = .5) +
  geom_smooth(method = lm, formula = 'y~0+x', se = TRUE)
ggsave('../plots/regression_fit.png')