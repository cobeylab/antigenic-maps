parse_numeric_titers <- function(titer){
  ## For titers reported in the supplement of Bedford et al. 2014:
  ##  If a titer is reported as "<x", report a numeric point estimate equal to x/2 
  ## Check that all titers are either numeric or contain "<"
  stopifnot(all(is.numeric(as.numeric(titer)) | grepl('<', titer)))
  
  parsed_titers = as.numeric(titer)
  for(ii in 1:length(parsed_titers)){
    if(is.na(parsed_titers[ii])){
      ## For "<" titers, extract the upper bound, and report half its value
      parsed_titers[ii] = gsub(pattern = '<(.+)', replacement = '\\1',x = titer[ii]) %>%
        as.numeric()/2
    }
  }
  stopifnot(all(is.numeric(parsed_titers)))
  parsed_titers
}


