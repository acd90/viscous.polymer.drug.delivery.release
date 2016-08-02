# Compare script to COPASUtils


# ---------------------------------------------------------------------
#                           File import
# ---------------------------------------------------------------------
import_96wellplates <- function(filepath, pattern, remove = 0, skip = 16,  ...){
  require(readxl)
  require(dplyr)
  require(tidyr)
  filelist <- list.files(path       = filepath,
                         pattern    = pattern, 
                         full.names = T)
 # Sometimes theres files that make your life miserable;
 # remove takes care of them
   if(remove != 0){
    list.import <- lapply(filelist[-remove], read_excel, skip = skip,...)
    filelist <- filelist[-remove]
  }else{
    list.import <-  lapply(filelist, read_excel , skip= skip,...)
  }
  
well <- paste0(LETTERS[1:8], rep(1:12,each=8))

 # Binding all imported list and selecting the squareish data frame  
wide_df <-  list.import %>% 
  lapply(., '[',1:8 , 2:13) %>% 
  lapply(.,unlist) %>% 
  do.call(rbind, .) %>% 
  data.frame(filename = basename(filelist), .) 

names(wide_df) <- c("filename", well)
tidy_df <- gather(wide_df, key = "well", value = "absorbance", 2:97, -filename)
tidy_df[] <-lapply(tidy_df, function(x) if(is.factor(x)==T) as.character(x) else x) 
return(tidy_df)
}
# ---------------------------------------------------------------------
#                         Annotation
# ---------------------------------------------------------------------
wells <- function(initialrow, finalrow, initialcolumn, finalcolumn){
  wells <- paste0(rep(LETTERS[initialrow:finalrow], each = finalcolumn-initialcolumn+1),
                  rep(initialcolumn:finalcolumn, times = finalrow-initialrow+1))
  return(wells)
}

multi_join <- function(..., by){

  jointhis <- list(...) 
  joint_df <- Reduce(function(x,y) merge(x,y, by = by),jointhis)

  return(joint_df)}

get_time <- function(mjoin_df, pattern = "(?<=D)[0-9]+"){
  library(tidyr)
  library(dplyr)
  library(stringr)
  library(gtools)
  
  # get_time takes a mjoin_df output from import_96wellplates
  # extract the filenames and turns the filename into numeric
  # with a dayplate column it determines to which day each well 
  # corresponds in a file
  # dayplate is a dataframe with a column of values and a column of wells
  mjoin_df <- mjoin_df[mixedorder(mjoin_df$well), ]
  
  timelist <- mjoin_df$filename %>%
    str_extract_all(., pattern = pattern) %>% 
    lapply(., as.numeric)
  
  mjoin_df <- mutate(mjoin_df, TID = 1:nrow(mjoin_df))
  
  timetest <- vector("list", length = nrow(mjoin_df))
  
  for(i in 1:nrow(mjoin_df)){
    timetest[[i]] <-timelist[[i]][mjoin_df$Day.Plate[i]]
  }
  
  df_time <- mjoin_df %>% 
    data.frame(., Time = as.numeric(unlist(timetest)), stringsAsFactors = F) 
  
  return(df_time)
}

annotate96 <- function(wellgrid, values, varname){
  library(tidyr)
  library(stringr)
  library(dplyr)
  
  # This function allows you to annotate 96 well plates in a more readable fashion
  # wellgrid is a character rectangular grid of the kind A1:H12
  # more than one well grid can be supplied in a given character
  # ie. wellgrid = "A1:B3, A3:B6", value = 1
  # value is the iterated value for that grid
  # the vector length of value and wellgrid must be the same
  # ie. wellgrid = c("A3:B6, A4:B10", "E7:H7"), value = c(1,2)
  
  # values/variables need better names
  
  stringvector <- wellgrid %>% 
    str_extract_all(., pattern = "[A-Ha-h][0-9]+\\:[A-Ha-h][0-9]+")
  
  # ULQ stands for Upper Left Quadrant  (Top left)
  ULQ <- stringvector %>% 
    str_extract_all(., pattern = "(?<!\\:)[a-hA-H][0-9]+") # Negative lookbehind
  
  
  # LRQ stands for Lower Right Quadrant (Bottom Right)
  LRQ <- stringvector %>% 
    str_extract_all(., pattern = "(?<=[0-9]\\:)[a-hA-H][0-9]+") # ?<= lookbehind assertion
  
  initialrow <- ULQ %>% 
    str_extract_all(., pattern = "[A-H]" ) %>% 
    lapply(., match, table = LETTERS) 
  
  
  finalrow <- LRQ %>% 
    str_extract_all(., pattern = "[A-H]") %>% 
    lapply(., match, table = LETTERS)  
  
  initialcolumn <- ULQ %>%
    str_extract_all(., pattern = "[0-9]+") %>% 
    lapply(., as.numeric) 
  
  finalcolumn <- LRQ %>% 
    str_extract_all(., pattern = "[0-9]+") %>% 
    lapply(., as.numeric) 
  
  # Initializing values
  listtest <- vector("list", length(wellgrid))
  value <- vector("list", length(wellgrid))
 
   plate_96wells <- wells(initialrow = 1,
                         finalrow = 8, 
                         initialcolumn = 1,
                         finalcolumn = 12) %>%
    data.frame(well = ., stringsAsFactors = F)
  
  for(i in 1:length(wellgrid)){
    listtest[[i]] <- Map(wells, initialrow[[i]], finalrow[[i]], initialcolumn[[i]], finalcolumn[[i]])
    value[[i]] <- rep(values[i], length(unlist(listtest[[i]])))
  }
  
  df <- data.frame(well = unlist(listtest),
                   values = unlist(value),
                   stringsAsFactors = F) 
  
  annotated_df<- full_join(df, plate_96wells, by = "well")
  names(annotated_df)[2] <- varname
  return(annotated_df)
}

# ---------------------------------------------------------------------
#                         Calculations
# ---------------------------------------------------------------------

standard_curve <- function(df){ 
  library(dplyr)
  library(broom)
  
  b.m  <- df %>%
    do(model = lm(absorbance~concentration, data = .)) %>% 
    broom::tidy(model) %>%
    select(1:3) %>% 
    spread(key = term, estimate)

  r.sq <- calibration_df %>% do(model = lm(absorbance~concentration, data = .)) %>%
    glance(model) %>% 
    select(1:2)
  sc <- inner_join(b.m, r.sq, by = "caltime")
  names(sc) <- c("caltime", "intercept", "slope", "r.squared")
  return(sc)
}

dydx <- function(y,x){
  dydx <- (y-lag(y))/(x-lag(x))
  return(dydx)
}

mass_conservation <- function(cumulative.release, timepoint.release) {
  index <- timepoint.release<0
  sumindex <- sum(index)
  
  while(isTRUE(sumindex > 0)){
    index <- timepoint.release<0
    leadindex <- lead(index,default = F)
    cumulative.release[index] <- cumulative.release[leadindex]
    cumulative.release <- cumulative.release-lag(cumulative.release,default = 0) 
    sumindex <- sum(index)
  }
  return(cumulative.release)
}



drug_release <- function(df, intercept, slope, vial.volume, sample.volume){
  library(dplyr)
  # df is a data.frame with filename, wellname/SID, and absorbance/96wellouput
  # slope, vial.volume and sample.volume must be dimensionally correct (have the same units)

  release_df <- df %>%
    mutate(concentration = (absorbance - intercept)/slope) %>%
    mutate(drug.in.vial = concentration*vial.volume) %>%
    mutate(drug.removed = concentration*sample.volume) %>%
    mutate(total.removed = cumsum(drug.removed)) %>%         
    mutate(cumulative.release = drug.in.vial + total.removed - drug.removed) %>%
    mutate(timepoint.release = cumulative.release - lag(cumulative.release,default = F)) 
  return(release_df)}
  


        
# ---------------------------------------------------------------------
#                         Statistics
# ---------------------------------------------------------------------
delivery_stats <- function(grouped_df, column){
  
  stats  <- summarise(grouped_df, N = length(column),
                      Average = mean(column),
                      StDev = sd(column))
  return(stats)                   
}

# ---------------------------------------------------------------------
#                        Plotting
# ---------------------------------------------------------------------

