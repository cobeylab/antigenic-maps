HaNam_titers_to_numeric = function(titer_col){
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
  extract(tibble(strain), strain, into = 'strain_year', regex = '.+/(\\d?\\d?\\d\\d$)', convert = T) %>%
    mutate(strain_year = clean_years(strain_year)) %>%
    pull(strain_year)
}


get_hanam_homologous_titer_for_individual = function(titers,  # column vector of titers to strains from a single year (group by year before passing titer inputs)
                                                     is_vietnam_titer # column vector indicating if the titer is to a vietnam isolate
){
  if(any(is_vietnam_titer)){ # If vietnam-specific strains from the year of interest exist in the panel, label them homologous
    homologous_titer = mean(titers[is_vietnam_titer])
  }else{ # If there are no Vietnam-specific strains from the year of interst, then any strain from the year of interest is considered homologous
    homologous_titer = mean(titers)
  }
  homologous_titer
}


get_HaNam_strains_chronologically <- function(){
  print('get_HaNam_strains_chronologically() is deprecated. Use get_strains_chronologicaly() instead.')
  read_csv('../processed_data/HaNam_strains.csv')  %>%
    tidyr::extract(col = test_strain, into = 'test_strain_year', regex = '/(\\d+)$', remove = F) %>%
    mutate(test_strain_year = clean_years(test_strain_year)) %>%
    arrange(test_strain_year) %>%
    pull(full_name)
}


get_vaccine_study_strains_chronologically <- function(study_year){
  print('get_vaccine_study_strains_chronologically() is deprecated. Use get_strains_chronologicaly() instead.')
  read_csv(paste0('../processed_data/vaccine_study_', study_year, '_strains.csv'), show_col_types = FALSE)  %>%
    arrange(test_strain_year) %>%
    pull(full_name)
}


get_strains_chronologically <- function(test_strain, # Character vector of strain names
                                        test_strain_year # Character vector of strain years
){
  tibble(test_strain = test_strain,
         test_strain_year = test_strain_year) %>%
    arrange(test_strain_year) %>%
    select(test_strain) %>%
    distinct() %>%
    pull(test_strain)
}


extract_pairwise_ditance <- function(strain1_name, 
                                     strain2_name, 
                                     ag_coords ## Data frame containing antigen names and corrdinates from Fonville 2014 Table S3
){
  ## Calculate and return the Euclidean distance between strains 1 and 2
  # Extract the coordinates for strain 1 as a vector
  v1 = ag_coords %>%
    filter(viruses == strain1_name) %>%
    select(-viruses) %>%
    unlist()
  # Extract the coordinates for strain 2 as a vector
  v2 = ag_coords %>%
    filter(viruses == strain2_name) %>%
    select(-viruses) %>%
    unlist()
  stopifnot(length(v1) == length(v2))
  # Get the Eudlidean distance and return
  sqrt(sum((v1-v2)^2))
}


calculate_and_save_map_distances <- function(){
  ag_coords = read_csv('../raw_data/Fonville_antigen_coords.csv')
  all_test_strains = bind_rows(read_csv('../processed_data/vaccine_study_1995_strains.csv'),
                               read_csv('../processed_data/vaccine_study_1997_strains.csv'),
                               read_csv('../processed_data/vaccine_study_2009_strains.csv'),
                               read_csv('../processed_data/vaccine_study_2010_strains.csv'),
                               read_csv('../processed_data/HaNam_strains.csv')) %>%
    mutate(full_name = toupper(full_name)) %>%
    distinct() %>% 
    arrange(full_name)
  
  ## Verify that all strains have estimated coordinates
  all_test_strains %>% 
    mutate(has_coord = full_name %in% ag_coords$viruses) %>%
    filter(!has_coord)
  
  ## Get pairwise distances
  ## This takes a minute to run and could probably be sped up
  dist_table = expand_grid(strain1 = ag_coords$viruses,
                           strain2 = ag_coords$viruses) %>%
    rowwise() %>%
    mutate(ferret_distance = extract_pairwise_ditance(strain1, strain2, ag_coords)) %>%
    ungroup()
  
  write_csv(dist_table, '../processed_data/calculated_strain_distances.csv')
  cat('saved ../processed_data/calculated_strain_distances.csv')
}



get_recent_ferret_distances <- function(sample_year,
                                        last_n_years,
                                        dist_table ## Table of antigenic distances between each pair of strains in the panel
){
  # Get the mean and median ferret (Euclidean) distances between each test strain in the panel, and recent strains circulating in the last n years
  # If last_n_years is 0, recent strains only include the previous season
  # If last_n_years is 1, recent strains include those from the 2 previous seasons...etc
  dist_table %>%
    # Extract circulation year from strain 1 names
    tidyr::extract(col = strain1, into = 'strain1_year', regex = '/(\\d+)$', remove = F) %>%
    mutate(strain1_year = clean_years(strain1_year)) %>%
    ## Keep strain 1s from the past n years
    filter(strain1_year %in% (sample_year - (0:(last_n_years)))) %>%
    distinct() %>%
    ## Calculate the mean, and median distance beteen valid strain1s and each strain2
    group_by(strain2) %>%
    summarise(mean_ferret_distance = mean(ferret_distance),
              median_ferret_distance = median(ferret_distance)) %>%
    mutate(sample_year = sample_year,
           last_n_years = last_n_years)
}



get_recent_individual_distances <- function(sample_year,
                                            last_n_years,
                                            titer_df ## Subset of the HaNam titers dataset excluding individuals with PCR-confirmed infections
){
  # Get the mean and median individual (logtiter) distances between each test strain in the panel, and recent strains circulating in the last n years
  foreach(subject_id = unique(titer_df$`Subject number`),
          .export = 'titer_df') %do% {
            xx = filter(titer_df,
                        `Subject number` == subject_id) %>%
              filter(`Sample year` == sample_year) %>%
              rename(sample_year = `Sample year`,
                     id = `Subject number`)
            
            merge(
              ## Recent strains
              xx %>%
                filter(test_strain_year %in% (sample_year - (0:(last_n_years)))) %>%
                rename(recent_strain = test_strain,
                       recent_strain_titer = logtiter) ,
              ## All panel strains
              xx %>% 
                rename(test_strain_titer = logtiter),
              by = c('id', 'YOB'), all = TRUE
            ) %>%
              ## For each recent-panel strain, get the titer difference
              mutate(pairwise_distance = recent_strain_titer - test_strain_titer) %>%
              group_by(test_strain) %>%
              summarise(median_dist_to_recent_strains = median(pairwise_distance), YOB = unique(YOB)) %>%
              ungroup() %>%
              mutate(`Subject number` = subject_id,
                     `Sample year` = sample_year,
                     last_n_years = last_n_years) 
          } %>%
    bind_rows()
}


plot_non_infections = function(n_years_back, ## How many years back from the sample year should be included in recent strain comparison?
                               titer_df, ## Data frame contaiing columns: `Sample year`, test_strain, test_strain_year, logtiter
                               plotname_flag,
                               dist_table
){
  stopifnot(all(c('Sample year', 'test_strain', 'test_strain_year', 'logtiter') %in% names(titer_df)))
  ## For each strain in the panel, plot the ferret distance to recent strains, and a violin plot of individual distances to recent strains
  sample_years = unique(titer_df$`Sample year`)
  xx = merge(
    ## recent_ferret_distances
    lapply(sample_years, FUN = function(ll) get_recent_ferret_distances(ll, last_n_years = n_years_back, dist_table = dist_table)) %>%
      bind_rows() %>%
      rename(test_strain = strain2),
    ## recent individual distances
    lapply(sample_years, FUN = function(ll) get_recent_individual_distances(ll, last_n_years = n_years_back, titer_df = titer_df)) %>%
      bind_rows() %>%
      tidyr::extract(col = test_strain, into = 'test_strain_year', regex = '/(\\d+)$', remove = F) %>%
      mutate(test_strain_year = clean_years(test_strain_year)) 
  ) 
  
  
  
  for(yy in sample_years){
    recent_start_year = yy - n_years_back
    pdat = xx %>%
      filter(sample_year == yy) %>%
      mutate(test_strain = factor(test_strain, levels = xx %>% select(starts_with('test_strain')) %>%
                                    arrange(test_strain_year) %>%
                                    distinct() %>%
                                    pull(test_strain)),
             xval = as.numeric(test_strain),
             last_n_years = factor(last_n_years),
             is_sample_year = sample_year == test_strain_year,
             YOB = as.factor(YOB),
             ylim_cutoff = ifelse(median_ferret_distance>15, 15, NA),
             adjusted_distance = adjust_ferret_distances(median_ferret_distance, `past/future`)) 
    rawplot = pdat %>%
      ggplot() +
      geom_violin(aes(x = test_strain, y = median_dist_to_recent_strains), draw_quantiles = .5) +
      geom_hline(aes(yintercept = 0), lty = 3)+
      geom_point(aes(x = jitter(xval), y = median_dist_to_recent_strains, color = YOB), alpha = .6) +
      geom_segment(aes(x = xval-.3, xend = xval+0.3, y = median_ferret_distance, yend = median_ferret_distance), color = 'indianred2')+
      geom_segment(aes(x = xval-.3, xend = xval+0.3, y = adjusted_distance, yend = adjusted_distance), color = 'gray')+
      ylim(c(-10, 15))+
      theme(axis.text.x = element_text(angle = 70, hjust = 1))+
      # facet_wrap(.~sample_year, labeller = label_both) +
      theme(legend.position = 'bottom')+
      ylab(sprintf('median distance to recent strains\n(%i seasons prior to sample)', n_years_back))+
      xlab('test_strain') +
      ggtitle(sprintf('%s from sample year: %i \nRecent strains include %i-%i', plotname_flag, yy, recent_start_year, yy)) +
      guides(color = guide_legend(nrow = 3, byrow = TRUE))+
      scale_color_viridis_d()
    if(any(!is.na(pdat$ylim_cutoff))) {
      rawplot +
        geom_segment(aes(x = xval, xend = xval, y = ylim_cutoff-2, yend = ylim_cutoff), color = 'indianred2', arrow = arrow(length = unit(0.05, "inches")))
    }else{
      rawplot
    }
    filename = sprintf('../plots/HaNam_%s_%iyearsback_sampleyear%i.png', plotname_flag, n_years_back, yy)
    cat(sprintf('saved plot to %s \n', filename))
    ggsave(filename = filename, width = 12, height = 9, units = 'in', dpi = 200)
  }
} 



prelim_vaccine_study_import <- function(csv_filename){
  ## This function imports vaccine study data frames (Tables S5, S6, S11, S12), fixes their names, which are spread across 2 rows in the original file, and returns a clean version of the data frame with unique column names
  
  repair_names = function(name, ## Original (duplicated name)
                          df){ # Raw data frame containing extra name info in row 1
    paste0(
      gsub(pattern = '(.+\\w)...\\d?\\d?\\d$', replacement = '\\1', x = name),
      '-',
      df[1,name]
    )
  }
  
  ## Import the raw data frame
  this_data_frame =read_csv(csv_filename)
  is_strain_column = grepl('/', x = names(this_data_frame))
  ## Fix the names
  names(this_data_frame) = c(names(this_data_frame)[!is_strain_column],
                             sapply(names(this_data_frame)[is_strain_column], FUN = function(xx){
                               repair_names(xx, this_data_frame)}
                             ))
  names(this_data_frame)[1] = 'subject'
  if(names(this_data_frame)[2] == 'Age'){
    names(this_data_frame)[2] = 'age'
  }else{
    this_data_frame = this_data_frame %>% mutate(age = NA)
  }
  this_data_frame %>% as_tibble()# return
}


process_vaccination_study_data <- function(titer_filename,
                                           strain_fiename,
                                           homologous_strain_name){
  
  study_reftable <- tibble(vaccine_strain_year = c(1995, 1997, 2009, 2010),
                           sample_year = c(1997, 1998, 2009, 2010),
                           vaccine_strain = c('A/NANCHANG/933/95', 'A/SYDNEY/5/97', 'A/PERTH/16/2009', 'A/VIETNAM/53/2010'))
  stopifnot(homologous_strain_name %in% study_reftable$vaccine_strain)
  
  
  prelim_vaccine_study_import(titer_filename) %>%  ## Import raw vaccine study titer csv
    slice(-1) %>% ## Drop the extraneous row with pre/post column name labels
    mutate(across(.cols = contains('/'), .funs = ~HaNam_titers_to_numeric(.))) %>% ## Convert titer values to numeric. Report undetectable = 5 --> logtiter = 0
    pivot_longer(contains('/'), names_to = c('test_strain', 'timepoint'), values_to = 'titer', names_sep = '-') %>% ## Pivot to long
    filter(!test_strain == 'A/PHILIPPINES/472/2002_EGG') %>% ## Filter out the Egg adapted version of A/PHILIPPINES/472/2002
    mutate(test_strain = ifelse(test_strain == 'A/PHILIPPINES/427/2002_MDCK', 'A/PHILIPPINES/472/2002', test_strain)) %>% ## Rename the MDCK-passaged version using the standard format
    mutate(vaccine_strain = homologous_strain_name) %>%
    merge(study_reftable, by = 'vaccine_strain') %>% # Merge in notes on study year and smaple year 
    mutate(YOB = as.numeric(sample_year) - as.numeric(age)) %>%
    merge(read_csv(strain_fiename)) %>% ## Merge with a strain name dictionary in which I've manually standardized the names to follow the same format
    mutate(test_strain = full_name) %>% ## Replace names with human-edited names in the standard format
    select(-full_name) %>% ## Drop the redundant column from the dictionary
    mutate(titer = as.numeric(titer),
           logtiter = to_log(titer)) %>% ## Get the log titer
    tidyr::extract(col = test_strain, into = 'test_strain_year', regex = '/(\\d+)$', remove = F) %>%
    mutate(test_strain_year = clean_years(test_strain_year))  %>%  ## Extract and clean the year of test strain circulation from the strain name
    group_by(subject, YOB, timepoint) %>% ## Get the homologous titer and distance
    mutate(is_homologous_strain = test_strain == homologous_strain_name,
           homologous_strain = homologous_strain_name,
           homologous_titer = mean(logtiter[is_homologous_strain]),
           homologous_strain_year = extract_strain_year(homologous_strain),
           individual_distance = homologous_titer-logtiter,
           dy = test_strain_year-homologous_strain_year,
           `past/future` = ifelse(dy<0, 'past', 'future'))  %>%
    ## Merge with ferret distances
    merge(dist_table %>% 
            dplyr::filter(strain1 == homologous_strain_name) %>% 
            rename(test_strain = strain2, 
                   homologous_strain = strain1),
          by = c('homologous_strain', 'test_strain'), all_x = TRUE, all_y = FALSE
    )
  
}


plot_vaccine_study_by_strain <- function(titer_data_frame,
                                         study_year,
                                         return_plot = TRUE){
  plot_data = titer_data_frame %>%
    ungroup() %>%
    mutate(subject = as.factor(subject),
           test_strain = factor(test_strain, levels = get_strains_chronologically(test_strain, test_strain_year)),
           test_strain_year = as.factor(test_strain_year),
           x_val = jitter(as.numeric(test_strain)),
           xval_homologous_label = ifelse(is_homologous_strain, test_strain, NA),
           ylim_cutoff = ifelse(ferret_distance>15, 15, NA),
           adjusted_distance = adjust_ferret_distances(ferret_distance, `past/future`))
  
  count_table = plot_data %>%
    group_by(test_strain, timepoint) %>%
    dplyr::summarise(n_individuals = n()) %>%
    ungroup() %>%
    mutate(label = sprintf('n=%i', n_individuals))
  
  this_plot =  plot_data %>%
    ggplot() +
    geom_violin(aes(x = test_strain, y = individual_distance)) +
    geom_point(aes(x = x_val, y = jitter(individual_distance), color = YOB), alpha = .6) +
    geom_violin(aes(x = test_strain, y = individual_distance), fill = NA) +
    geom_point(aes(x = xval_homologous_label, y = -7.5), color = 'red', pch = 8)+
    geom_segment(aes(x = as.numeric(test_strain)-0.3, xend = as.numeric(test_strain)+0.3, y = ferret_distance, yend = ferret_distance), color = 'indianred3')+
    geom_segment(aes(x = as.numeric(test_strain)-0.3, xend = as.numeric(test_strain)+0.3, y = adjusted_distance, yend = adjusted_distance), color = 'gray')+
    geom_text(aes(x = test_strain, y = -8.5, label = label), data = count_table, angle = 45, size = 2.5) +
    geom_hline(aes(yintercept = 0), lty = 3) +
    facet_wrap(timepoint~., labeller = 'label_both',ncol = 1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(c(-9.5, max(plot_data$individual_distance))) +
    ylab('individual distance\n[homoogous titer]-[heterologous titer]')+
    xlab('test strain') +
    scale_color_viridis_c() +
    ylim(c(-10, 15))
  if(any(!is.na(plot_data$ylim_cutoff))) {
    this_plot = this_plot +
      geom_segment(aes(x = test_strain, xend = test_strain, y = ylim_cutoff-2, yend = ylim_cutoff), color = 'indianred2', arrow = arrow(length = unit(0.05, "inches"))) +
      geom_hline(aes(yintercept = 15), color = 'indianred2', lty = 3)
  }
  filename = paste0('../plots/vaccine_study_', study_year, '_bystrain.png')
  ggsave(filename = filename, width = 17, height = 9, units = 'in', dpi = 200)
  cat(sprintf('saved plot to %s', filename))
  if(return_plot == TRUE) return(this_plot)
}



plot_vaccine_study_by_year <- function(titer_data_frame,
                                       study_year,
                                       return_plot = TRUE){
  plot_data = titer_data_frame %>%
    ungroup() %>%
    mutate(subject = as.factor(subject),
           test_strain_year = factor(test_strain_year),
           x_val = jitter(as.numeric(test_strain_year)),
           xval_homologous_label = ifelse(test_strain_year==study_year, test_strain_year, NA),
           ylim_cutoff = ifelse(ferret_distance>15, 15, NA),
           adjusted_distance = adjust_ferret_distances(ferret_distance))
  
  count_table = plot_data %>%
    group_by(test_strain_year, timepoint) %>%
    dplyr::summarise(n_individuals = n(),
                     n_test_strains = length(unique(test_strain))) %>%
    ungroup() %>%
    mutate(label_indiv = sprintf('n=%i', n_individuals),
           label_strains = sprintf('s=%i', n_test_strains))
  
  this_plot = plot_data %>%
    ggplot() +
    geom_violin(aes(x = as.factor(test_strain_year), y = individual_distance), draw_quantiles = c(0.5)) +
    #geom_boxplot(aes(x = test_strain_year, y = ferret_distance), color = 'indianred3', fill = NA)+
    geom_segment(aes(x = as.numeric(test_strain_year)-.2, xend = as.numeric(test_strain_year)+.2, y = ferret_distance, yend = ferret_distance), color = 'indianred3')+
    geom_segment(aes(x = as.numeric(test_strain_year)-.2, xend = as.numeric(test_strain_year)+.2, y = adjusted_distance, yend = adjusted_distance), color = 'gray')+
    geom_point(aes(x = x_val, y = jitter(individual_distance), color = YOB), alpha = .3) +
    geom_point(aes(x = xval_homologous_label, y = -7.5), color = 'red', pch = 8)+
    #geom_segment(aes(x = as.numeric(test_strain_year)-0.3, xend = as.numeric(test_strain_year)+0.3, y = ferret_distance, yend = ferret_distance), color = 'indianred3')+
    geom_text(aes(x = test_strain_year, y = -8.5, label = label_indiv), data = count_table, angle = 45, size = 2.5) +
    geom_text(aes(x = test_strain_year, y = max(plot_data$individual_distance)-5, label = label_indiv), data = count_table, angle = 45, size = 2.5, color = 'indianred3') +
    geom_hline(aes(yintercept = 0), lty = 3) +
    facet_wrap(timepoint~., labeller = 'label_both',ncol = 1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(c(-9, max(plot_data$individual_distance))) +
    ylab('individual distance\n[homoogous titer]-[heterologous titer]')+
    xlab('test strain') +
    scale_color_viridis_c()
  if(any(!is.na(plot_data$ylim_cutoff))) {
    this_plot = this_plot +
      geom_segment(aes(x = test_strain_year, xend = test_strain_year, y = ylim_cutoff-2, yend = ylim_cutoff), color = 'indianred2', arrow = arrow(length = unit(0.05, "inches"))) +
      geom_hline(aes(yintercept = 15), color = 'indianred2', lty = 3)+
      ylim(c(-10, 15))
  }
  filename = paste0('../plots/vaccine_study_', study_year, '_byyear.png')
  ggsave(filename = filename, width = 17, height = 9, units = 'in', dpi = 200)
  cat(sprintf('saved plot to %s', filename))
  if(return_plot == TRUE) return(this_plot)
}


get_xval_for_year_of_birth = function(test_strain, # vector of test strain names
                                      subject_id, # unique subject id
                                      YOB){ # unique year of birth
  ## This function inputs columns from a plotting data frame grouped by subject, in which test_strain and test_strain_year are factors ordered chronologically
  ## It outputs the x-value at which to plot a vertical line indicating the subject's YOB
  ## We aim to plot the vertical line immediately to the left of the frist strain in the sero panel that circulated after the subject was born
  stopifnot(length(unique(YOB))==1) ## Make sure there's only one YOB 
  stopifnot(is.factor(test_strain))
  
  tibble(test_strain = factor(levels(test_strain), levels = levels(test_strain)), ## using levels instead of values ensures all strains are included, even if this subject doesn't have a titer measurement
         subject_id = subject_id,
         x_val = as.numeric(test_strain),
         YOB = unique(YOB)) %>% ## Extract test_strain_years
    tidyr::extract(col = test_strain, into = 'test_strain_year', regex = '/(\\d+)$', remove = F) %>%
    mutate(test_strain_year = clean_years(test_strain_year)) %>%
    arrange(test_strain_year) %>% ## Make sure the data are arranged chronologically
    distinct() %>%
    mutate(dx = c(1, diff(x_val))) %>%
    filter(test_strain_year >= YOB) %>% ## Only keep strains that circulated in the individual's lifetime
    slice(1) %>% ## Keep the first strain
    mutate(this_xval = x_val - (dx/2)) %>%  ## Return the first post-birth strain's xvalue, shifted left 
    pull(this_xval)
}





adjust_ferret_distances <- function(ferret_map_distance, 
                                    past_future,
                                    lm_fit_file = '../../Fonville_2016/processed_data/primary_ferret_human_distance_lm.rds'){
  
  ## Using a previously fitted regression: 'ferret_overestimate = 0 + beta*ferret_map_distance', predict the individual distance as:
  ##   ferret_map_distance - ferret_overestimate
  if(!file.exists(lm_fit_file)){
    stop(sprintf('regression file %s does not exist. This WD is: %s', lm_fit_file, getwd()))
    # setwd('../../Fonville_2016/R/')
    # source('regressions.R')
    # cat('Sourcing ../../Fonville_2016/R/regressions.R\n')
    # cat('Fitting the regression: ferret_overestimate = 0 + beta*ferret_map_distance \n')
    # cat('Saving outputs as ../../Fonville_2016/processed_data/primary_ferret_human_distance_lm.rds\n')
    # lm_fit_file = '../../Fonville_2016/processed_data/primary_ferret_human_distance_lm.rds'
  }
  ferret_overestimate = predict.lm(read_rds(lm_fit_file), newdata = tibble(ferret_map_distance = ferret_map_distance,
                                                                           `past/future` = past_future))
  ferret_map_distance - ferret_overestimate ## Return adjusted ferret distance
}












plot_vaccine_study_by_age <- function(titer_data_frame,
                                         study_year,
                                         return_plot = TRUE){
  plot_data = titer_data_frame %>%
    filter(timepoint == 'POST') %>%
    ungroup() %>%
    mutate(subject = as.factor(subject),
           `age<=20` = YOB>=1977,
           age_offset = ifelse(`age<=20`, -.2, .2),
           test_strain = factor(test_strain, levels = get_strains_chronologically(test_strain, test_strain_year)),
           test_strain_year = as.factor(test_strain_year),
           x_val = jitter(as.numeric(test_strain))+age_offset,
           xval_homologous_label = ifelse(is_homologous_strain, test_strain, NA),
           ylim_cutoff = ifelse(ferret_distance>15, 15, NA),
           adjusted_distance = adjust_ferret_distances(ferret_distance, `past/future`))
  
  count_table = plot_data %>%
    group_by(test_strain,  `age<=20`) %>%
    dplyr::summarise(n_individuals = n()) %>%
    ungroup() %>%
    mutate(label = sprintf('n=%i', n_individuals))
  
  this_plot =  plot_data %>%
    ggplot() +
    geom_violin(aes(x = test_strain, y = individual_distance)) +
    geom_point(aes(x = x_val, y = jitter(individual_distance), color = YOB), alpha = .6) +
    geom_violin(aes(x = test_strain, y = individual_distance), fill = NA) +
    geom_point(aes(x = xval_homologous_label, y = -7.5), color = 'red', pch = 8)+
    geom_segment(aes(x = as.numeric(test_strain)-0.3, xend = as.numeric(test_strain)+0.3, y = ferret_distance, yend = ferret_distance), color = 'indianred3')+
    geom_segment(aes(x = as.numeric(test_strain)-0.3, xend = as.numeric(test_strain)+0.3, y = adjusted_distance, yend = adjusted_distance), color = 'gray')+
    geom_text(aes(x = test_strain, y = -8.5, label = label), data = count_table, angle = 45, size = 2.5) +
    geom_hline(aes(yintercept = 0), lty = 3) +
    facet_wrap(`age<=20`~., labeller = 'label_both',ncol = 1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(c(-9.5, max(plot_data$individual_distance))) +
    ylab('individual distance\n[homoogous titer]-[heterologous titer]')+
    xlab('test strain') +
    scale_color_viridis_c() +
    ylim(c(-10, 15))
  if(any(!is.na(plot_data$ylim_cutoff))) {
    this_plot = this_plot +
      geom_segment(aes(x = test_strain, xend = test_strain, y = ylim_cutoff-2, yend = ylim_cutoff), color = 'indianred2', arrow = arrow(length = unit(0.05, "inches"))) +
      geom_hline(aes(yintercept = 15), color = 'indianred2', lty = 3)
  }
  filename = paste0('../plots/vaccine_study_', study_year, '_byage.png')
  ggsave(filename = filename, width = 17, height = 9, units = 'in', dpi = 200)
  cat(sprintf('saved plot to %s', filename))
  if(return_plot == TRUE) return(this_plot)
}


