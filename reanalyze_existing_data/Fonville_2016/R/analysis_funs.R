titers_to_numeric = function(titer_col){
  ## Input a column vector containing titer observations
  ## Assign undetectable ("<10") and endpoint titers (">=1280") a numeric value
  titer_col = ifelse(titer_col == '*', NA, titer_col)
  titer_col = ifelse(titer_col == '<10', 5, titer_col)
  titer_col = ifelse(titer_col == '>=1280', 1280, titer_col)
  as.numeric(titer_col)
}


to_log = function(titer){
  ## Convert to log2 titer
  ## In Fonville et al they define log titer as log2(titer/10)
  stopifnot(is.numeric(titer))
  log2(titer/10)
}


clean_years = function(years){
  ## The years reported at the end of strain names follow all different formats
  ## This converts strings to numeric and converts 2-digit abbreviated years to 4-digit years
  years = as.numeric(years)
  ifelse(years<20, years+2000, 
         ifelse(years<100, years+1900,
                years))
}

extract_strain_year = function(strain){
  extract(tibble(strain), strain, into = 'strain_year', regex = '.+/(\\d\\d\\d?$)', convert = T) %>%
    mutate(strain_year = clean_years(strain_year)) %>%
    pull(strain_year)
}



import_ferret_titers = function(){
  ferret_titers = read_csv('../raw_data/ferret_titers_Dataset_S3.csv', skip = 1)
  names(ferret_titers)[1] = 'strain' ## Name the strain column
  ferret_titers %>% pivot_longer(-strain, names_to = 'serum', values_to = 'titer') %>% # Convert to long format
    mutate(strain_year = extract_strain_year(strain),
           titer = titers_to_numeric(titer), ## Extract numeric log and absolute titers
           logtiter = to_log(titer)) %>%
           mutate(strain = factor(strain, levels = get_strains_chronologically(strain, strain_year)))
}


import_human_titer_data <- function(){
  raw_dat = read_csv('../raw_data/human_titers_Dataset_S2.csv', skip = 2) %>%
    t()
  colnames(raw_dat) = raw_dat[1,]
  raw_dat = raw_dat[-1,]
  subject_ids = rownames(raw_dat)
  raw_dat = as_tibble(raw_dat) %>%
    mutate(subject = subject_ids) %>%
    pivot_longer(contains('/'), names_to = 'strain', values_to = 'titer') %>%
    mutate(strain_year = extract_strain_year(strain),
           titer = titers_to_numeric(titer),
           logtiter = to_log(titer),
           `age (months)` = as.numeric(`age (months)`)) %>%
    rename(age_mos = `age (months)`) %>%
    mutate(strain = factor(strain, levels = get_strains_chronologically(strain, strain_year)),
           is_high_responder = subject %in% c('S95-1',	'S00-1',	'S00-2',	'S04-1',	'S04-2',	'S04-3'))
  raw_dat
}


get_euclidean_distance = function(v1, v2){
  sqrt(sum((v1-v2)^2))
}

get_one_distance <- function(strain1, strain2, ag_coords){
  v1 = ag_coords[ag_coords$strain == strain1, 2:3]
  v2 = ag_coords[ag_coords$strain == strain2, 2:3]
  get_euclidean_distance(v1, v2)
}

import_map_distances <- function(){
  ag_coords = read_csv('../processed_data/kg_scraped_ferret_coordinates_Fig3.csv')
  ## Get all possible combinations of strains with coordinates
  with(ag_coords %>% filter(!is.na(x)),
       expand_grid(strain1 = strain,
                   strain2 = strain)) %>%
    rowwise() %>%
    mutate(distance = get_one_distance(strain1, strain2, ag_coords))
}

get_strains_chronologically <- function(strain, strain_year){
  tibble(strain = strain, 
         strain_year = strain_year) %>%
    distinct() %>%
    arrange(strain_year) %>%
    pull(strain)
}


get_putative_primary_strains = function(strain, titer, circulated_in_sampling_window){
  tibble(strain, titer, circulated_in_sampling_window) %>%
    filter(circulated_in_sampling_window) %>%
    filter(titer == max(titer)) %>%
    pull(strain)
}


get_subject_factor_by_YOB <- function(subject, YOB, season){
  tibble(subject = subject,
         YOB = YOB, 
         season = season) %>%
    distinct() %>%
    arrange(YOB) %>%
    mutate(label = get_subject_label(subject, YOB, season)) %>%
    pull(label)
}

get_subject_label <- function(subject, YOB, season){
  tibble(subject = subject,
         YOB = YOB, 
         season = season) %>%
    mutate(label = sprintf('%s, YOB:%s, S:%s', subject, YOB, season)) %>%
    pull(label)
}


adjust_ferret_distances <- function(ferret_map_distance, 
                                    lm_fit_file = '../processed_data/primary_ferret_human_distance_lm.rds'){
  
  ## Using a previously fitted regression: 'ferret_overestimate = 0 + beta*ferret_map_distance', predict the individual distance as:
  ##   ferret_map_distance - ferret_overestimate
  if(!file.exists(lm_fit_file)){
    source('regressions.R') 
    cat('Sourcing regressions.R\n')
    cat('Fitting the regression: ferret_overestimate = 0 + beta*ferret_map_distance \n')
    cat('Saving outputs as ../processed_data/primary_ferret_human_distance_lm.rds\n')
    lm_fit_file = '../processed_data/primary_ferret_human_distance_lm.rds'
  }
  ferret_map_distance - predict.lm(read_rds(lm_fit_file), newdata = tibble(ferret_map_distance = ferret_map_distance))
}
