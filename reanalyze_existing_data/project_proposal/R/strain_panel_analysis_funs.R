## Import data and format into the following columns:

# * dataset - dataset name
# * id - participant id
# * homologous_strain - name of strain responsible for recent infection or vaccination
# * homologous_strain_year - year in which homologous strain circulated
# * test_strain - name of strain to which titer is measured
# * test_strain_year - year in which the test strain circulated
# * dt - temporal distance (years) between the homologous strain and vaccine strain
# * past_future - 'past' if the test strain circulated before the homologous strain. 'future' otherwise.
# * logtiter - observed log titer
# * host_type - ferret or human
# * exposure_no - primary or secondary
# * exposure_type - vaccine or infection

to_log = function(titer){
  log2(titer/10)
}


select_standard_cols = function(df){
  df %>%
    select(id, 
           homologous_strain_year, 
           test_strain, 
           test_strain_year, 
           logtiter, 
           is_cross_reactive,
           dt, 
           past_future, 
           host_type, 
           exposure_no, 
           exposure_type, 
           dataset) %>%
    mutate(id = as.character(id), 
           homologous_strain_year = as.integer(homologous_strain_year), 
           test_strain = as.character(test_strain), 
           test_strain_year = as.integer(test_strain_year), 
           logtiter = as.numeric(logtiter), 
           is_cross_reactive = as.logical(is_cross_reactive),
           dt = as.numeric(dt), 
           past_future = factor(past_future, levels = c('past', 'future', 'concurrent')), 
           host_type = as.character(host_type), 
           exposure_no = as.character(exposure_no), 
           exposure_type = as.character(exposure_type), 
           dataset = as.character(dataset))
}

get_past_future = function(homologous_strain_year, test_strain_year){
  ifelse(homologous_strain_year>test_strain_year, 'past', 
         ifelse(homologous_strain_year<test_strain_year, 'future', 'concurrent'))
}

import_HaNam_PCR_infections <- function(cross_reactive_threshold = 40){
  source('../../Ha_Nam/R/analysis_functions.R')
  raw_df = read_csv('../../Ha_Nam/raw_data/Fonville_2014_TableS3.csv', show_col_types = F) %>% ## HaNam cohort titers
    mutate_at(-c(1,2), .funs = ~HaNam_titers_to_numeric(.)) %>% ## Convert titer values to numeric. Report undetectable = 5, and endpoint = 1280
    pivot_longer(-c(1,2), names_to = 'test_strain', values_to = 'titer') %>%
    dplyr::filter(!test_strain == 'A/PHILIPPINES/472/2002_EGG') %>% ## Filter out the Egg adapted version of A/PHILIPPINES/472/2002
    mutate(test_strain = ifelse(test_strain == 'A/PHILIPPINES/427/2002_MDCK', 'A/PHILIPPINES/472/2002', test_strain)) %>%
    mutate(logtiter = to_log(titer)) %>%
    tidyr::extract(col = test_strain, into = 'test_strain_year', regex = '/(\\d+)$', remove = F) %>%
    mutate(test_strain_year = clean_years(test_strain_year)) %>%
    merge(read_csv('../../Ha_Nam/processed_data/YOB_inferred.csv', show_col_types = F), by = 'Subject number') %>% ## Merge with Kangchon's annotation of PCR-confirmed infections %>%
    merge(read_csv('../../Ha_Nam/processed_data/HaNam_strains.csv', show_col_types = F), by = 'test_strain') %>%
    mutate(test_strain = full_name) %>%
    dplyr::filter(`PCR+` == 'PCR+') %>%
    dplyr::filter(`Sample year` == `Infection Time`) %>%
    dplyr::filter(!is.na(titer)) %>%
    mutate(is_infecting_year = test_strain_year == `Infection Time`) %>%
    group_by(`Subject number`) %>%
    mutate(is_homologous_strain = if(any(is_vietnam_strain & is_infecting_year)){ is_vietnam_strain & is_infecting_year } else {is_infecting_year},
           homologous_titer = mean(logtiter[is_homologous_strain])) %>%
    rename(id = `Subject number`,
           sample_year = `Sample year`,
           homologous_strain_year = `Infection Time`) %>%
      mutate(dataset = 'HaNam_PCR_infections',
             dt = abs(homologous_strain_year-test_strain_year),
             past_future = get_past_future(homologous_strain_year, test_strain_year),
             host_type = 'human',
             exposure_no = 'secondary',
             exposure_type = 'infection',
             is_cross_reactive = ifelse(logtiter >= to_log(cross_reactive_threshold), TRUE, FALSE)) 
    select_standard_cols(raw_df)
}

import_HaNam_seroconversions <- function(cross_reactive_threshold = 40){
  source('../../Ha_Nam/R/analysis_functions.R')
  raw_df = read_csv('../../Ha_Nam/raw_data/Fonville_2014_TableS3.csv', show_col_types = F) %>% ## HaNam cohort titers
    mutate_at(-c(1,2), .funs = ~HaNam_titers_to_numeric(.)) %>% ## Convert titer values to numeric. Report undetectable = 5, and endpoint = 1280
    pivot_longer(-c(1,2), names_to = 'test_strain', values_to = 'titer') %>%
    dplyr::filter(!test_strain == 'A/PHILIPPINES/472/2002_EGG') %>% ## Filter out the Egg adapted version of A/PHILIPPINES/472/2002
    mutate(test_strain = ifelse(test_strain == 'A/PHILIPPINES/427/2002_MDCK', 'A/PHILIPPINES/472/2002', test_strain)) %>%
    mutate(logtiter = to_log(titer)) %>%
    tidyr::extract(col = test_strain, into = 'test_strain_year', regex = '/(\\d+)$', remove = F) %>%
    mutate(test_strain_year = clean_years(test_strain_year)) %>%
    merge(read_csv('../../Ha_Nam/processed_data/YOB_inferred.csv', show_col_types = F), by = 'Subject number') %>% ## Merge with Kangchon's annotation of PCR-confirmed infections %>%
    merge(read_csv('../../Ha_Nam/processed_data/HaNam_strains.csv', show_col_types = F), by = 'test_strain') %>%
    dplyr::filter(Seroconversion) %>%
    dplyr::filter(`Sample year` == `Infection Time`) %>%
    dplyr::filter(!is.na(titer)) %>%
    rename(id = `Subject number`,
           sample_year = `Sample year`,
           homologous_strain_year = `Infection Time`) %>%
    mutate(dataset = 'HaNam_seroconversions',
           dt = abs(homologous_strain_year-test_strain_year),
           past_future = get_past_future(homologous_strain_year, test_strain_year),
           host_type = 'human',
           exposure_no = 'secondary',
           exposure_type = 'infection',
           is_cross_reactive = ifelse(logtiter >= to_log(cross_reactive_threshold), TRUE, FALSE),
           homologous_strain = NA) 
  select_standard_cols(raw_df)
}

import_HaNam_vaccine_study <- function(cross_reactive_threshold = 40){
  source('../../Ha_Nam/R/analysis_functions.R')
 raw_df =  bind_rows(
    ## 1995 vaccine study
    read_csv('../../Ha_Nam/processed_data/vaccinations_1995.csv', col_types = 'ccdciiciiiiiliicd', show_col_types = F)[,1:12] %>%
      mutate(homologous_strain_year = 1995),
    read_csv('../../Ha_Nam/processed_data/vaccinations_1997.csv', col_types = 'cciciiciiiiili', show_col_types = F)[,1:12] %>%
      mutate(homologous_strain_year = 1997)
  ) %>%
    mutate(dataset = 'HaNam_vaccine',
           is_cross_reactive = ifelse(logtiter >= to_log(cross_reactive_threshold), TRUE, FALSE),
           dt = abs(homologous_strain_year - test_strain_year),
           past_future = get_past_future(homologous_strain_year, test_strain_year),
           exposure_no = 'secondary',
           exposure_type = 'vaccine',
           host_type = 'human',
           id = as.numeric(sprintf('%s-%i', subject, homologous_strain_year)))
 select_standard_cols(raw_df)
}

import_Fonville_2016_human_data <- function(cross_reactive_threshold = 40){
  source('../../Fonville_2016/R/analysis_funs.R')
  setwd('../../Fonville_2016/R/')
  raw_data = import_human_titer_data()
  setwd('../../project_proposal/R/')
  raw_data %>%
    tidyr::extract(season, into = 'homologous_strain_year', regex = '\\d{4}/(\\d{4})', convert = T) %>%
    rename(test_strain = strain,
           test_strain_year = strain_year) %>%
    mutate(id = subject,
           is_cross_reactive = ifelse(logtiter >= to_log(cross_reactive_threshold), TRUE, FALSE),
           dt = abs(homologous_strain_year - test_strain_year),
           past_future = get_past_future(homologous_strain_year, test_strain_year),
           host_type = 'human',
           exposure_no = 'primary',
           exposure_type = 'infection',
           dataset = 'Fonville_2016_children'
    ) %>%
      select_standard_cols()
  
}


import_Fonville_2016_ferret_data <- function(cross_reactive_threshold = 40){
  source('../../Fonville_2016/R/analysis_funs.R')
  setwd('../../Fonville_2016/R/')
  raw_data = import_ferret_titers()
  setwd('../../project_proposal/R/')
  raw_data %>%
    mutate(id = serum,
           homologous_strain_year = extract_strain_year(serum),
           test_strain = strain,
           test_strain_year = strain_year,
           is_cross_reactive = ifelse(logtiter >= to_log(cross_reactive_threshold), TRUE, FALSE),
           dt = abs(homologous_strain_year-test_strain_year),
           past_future = get_past_future(homologous_strain_year, test_strain_year),
           host_type = 'ferret',
           exposure_no = 'primary',
           exposure_type = 'infection',
           dataset = 'Fonville_2016_ferrets') %>%
    select_standard_cols()
}



import_Bedford_2014_data <- function(cross_reactive_threshold=40){
  parse_undetectable = function(und_titer){
    min_dilution = gsub('<(\\d+)', '\\1', und_titer) %>% as.numeric() ## Get the minimum measured dilution
    ifelse(min_dilution <= 10, 
           to_log((min_dilution)/2), ## Report half the measured min dilution
           ifelse(min_dilution <= 40, ## if between 10 and 40, report half the min dilution minus 0.5
                  to_log(min_dilution/2)-.5,
           NA)) ## Else don't report
  }
  titers_to_numeric_Bedford = function(titers){
    is_numeric = !grepl("<", titers)
    titers[is_numeric] = to_log(as.numeric(titers[is_numeric]))
    titers[!is_numeric] = parse_undetectable(titers[!is_numeric])
    as.numeric(titers)
  }
  read_delim(file = '../../map_refitting/raw-data/Bedford_2014_data.txt', show_col_types = F) %>%
    rename(test_strain = virusStrain,
           test_strain_year = virusYear,
           homologous_strain = serumStrain,
           homologous_strain_year = serumYear) %>%
    mutate(id = as.numeric(as.factor(sprintf('%s-%s', serumIsolate, source))),
           logtiter = titers_to_numeric_Bedford(titer),
           is_cross_reactive = ifelse(logtiter>=to_log(cross_reactive_threshold), TRUE, FALSE),
           dt = abs(homologous_strain_year-test_strain_year),
           past_future = get_past_future(homologous_strain_year, test_strain_year),
           host_type = 'ferret',
           exposure_no = 'primary',
           exposure_type = 'infection',
           dataset = 'Bedford_2014') %>%
    select_standard_cols() %>%
    filter(!is.na(logtiter))
}

plot_cr_bars <- function(df){
  df %>%
    filter(!is.na(logtiter)) %>%
    mutate(dt = homologous_strain_year-test_strain_year) %>%
    group_by(dt, is_cross_reactive) %>%
    summarise(count = n()) %>%
    group_by(dt) %>%
    mutate(fraction = count/sum(count),
           lower = fraction - 1.96*sqrt(fraction*(1-fraction)/sum(count)),
           upper = fraction + 1.96*sqrt(fraction*(1-fraction)/sum(count)),
           lower = ifelse(is_cross_reactive, lower, NA),
           upper = ifelse(is_cross_reactive, upper, NA))%>%
    ggplot() +
    geom_bar(aes(x = dt, y = fraction, fill = is_cross_reactive), stat = 'identity', alpha = .5)+
    geom_segment(aes(x = dt, xend = dt, y = lower, yend = upper))+
    geom_vline(aes(xintercept = 0), lty = 2) +
    ggtitle(sprintf('%s, %s response', unique(df$host_type), unique(df$exposure_no)), subtitle = sprintf('%s', unique(df$dataset)))+
    xlim(c(-40, 40))
}




