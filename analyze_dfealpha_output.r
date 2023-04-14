rm(list=ls())
library(RColorBrewer)
library(tidyverse)
source("/nas/longleaf/home/adaigle/DFESelfing/calculate_pi.r")
library(reshape2)

dfealpha_dir <- "/nas/longleaf/home/adaigle/rerun_dfealpha/"
dfealpha_output_dirs <- paste(paste(dfealpha_dir, dir(dfealpha_dir, pattern = "DFE_alpha_output"), sep = ""),
    "/", dir(
    paste(dfealpha_dir, dir(dfealpha_dir, pattern = "DFE_alpha_output"), sep = ""), pattern = "selected")
    , sep = "")


#table of results with columns signifying selfing%, DFE, and output
dfealpha_colnames <- as.data.frame(c("N1", "N2", "t2", "Nw", "b", "Es", "f0","L"))
dfealpha_raw_results <- tibble(
    name = list.files(path = dfealpha_output_dirs, pattern = "est_dfe.out"),
    fullpath = paste(dfealpha_output_dirs, "/", name,sep=""),
    selfing = str_extract(fullpath, "(?<=DFE_alpha_output_)\\d+") %>% if_else(is.na(.), "0", .), #0's were empty
    matchname = str_extract(fullpath, "(DFE)\\d+(output)\\d"), #will be useful for comparing results later, or combining with other programs
    DFE = str_extract(matchname, "(DFE)\\d+"),
    data = lapply(fullpath, function(x) 
        as.data.frame(t(data.frame(matrix(unlist(strsplit(readLines(x), " ")), ncol = 2, byrow = TRUE)))) %>%
            rownames_to_column() %>%
            `colnames<-`(.[1,]) %>%
            .[-1,] %>%
            `rownames<-`(NULL))
)

df <- dfealpha_raw_results$data %>% bind_rows()
df <- df %>% mutate(row_number = row_number())
dfealpha_raw_results <- dfealpha_raw_results  %>% 
    mutate(row_number = row_number())

# Join the original table with the new data frame using the row numbers as the key
dfealpha_raw_results <- dfealpha_raw_results %>% 
    left_join(df, by = "row_number") %>% 
    select(-data, -row_number) %>%
    as.data.frame %>% 
    mutate_at(vars(7:14), as.numeric) %>%
    mutate(gamma = -2*Nw*Es) %>%
    rename(f0_fromoutput = f0) #f0 is a dfealpha output, but I also need to name a class f0 later

DFE_proportions_dfe_alpha <- function(meanGamma,beta, Nw) { 
    # modified function to calc vector of dfe classes given gamma and beta
    # assumes Nw is 100 bc grapes doesn't have two step
    # need to confirm it has no modified Ne
    meanS <- meanGamma/(2*Nw)

    # this code adapted from get_dfe_class_proportions
    s_shape <- beta
    s_rate <- s_shape/abs(meanS)
    x1 <- 1/(2*Nw)
    x10 <- 10/(2*Nw)
    x100 <- 100/(2*Nw)
    f0 <- pgamma(x1, shape=s_shape, rate=s_rate)
    f1 <- pgamma(x10, shape=s_shape, rate=s_rate) - f0
    f2 <- pgamma(x100, shape=s_shape, rate=s_rate) - f1 - f0
    f3 <- 1.0 - f2 - f1 - f0
    #print(paste(f0, f1, f2, f3, sep=", "))
    #return(c(list(f0), list(f1), list(f2), list(f3)))
    return(c(f0 = f0, f1 = f1, f2 = f2, f3 = f3))
} # for non DFEalpha 
# to do this ina  way that makes more sense
#piemp, get Nemperical, and use that ias your Nw for now one
true_gammas <- list(
  DFE1 = 5,
  DFE2 = 50,
  DFE3 = 1000
)
true_betas <- list(
  DFE1 = 0.9,
  DFE2 = 0.5,
  DFE3 = 0.3
)

dfealpha_raw_results_wtruth <- dfealpha_raw_results %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE]))

#now i run the class generator function on all outputs
dfealpha_raw_results_wclasses <- dfealpha_raw_results_wtruth %>%
  mutate(output = pmap(list(gamma, b, Nw), DFE_proportions_dfe_alpha)) %>%
  unnest_wider(output) %>%
  rename(z0 = f0) %>%  #had to rename columns bc i rerun the function
  rename(z1 = f1) %>% 
  rename(z2 = f2) %>% 
  rename(z3 = f3) %>% 
  #rename( c(t0,t1,t2,t3) = c(f0,f1,f2,f3)) %>%
  mutate(output = pmap(list(gamma, b, Nw), DFE_proportions_dfe_alpha)) %>%
  unnest_wider(output) 
  #mutate(across(.fns = ~ if(all(!is.na(as.numeric(.x)))) as.numeric(.x) else .x))

dfealpha_summary <- dfealpha_raw_results_wclasses %>% 
    group_by(selfing, DFE) %>%
    summarize(across(where(is.numeric), list(avg = mean, sd = sd))) %>%
    rename(z0 = z0_avg) %>%  #had to rename columns bc i rerun the function
    rename(z1 = z1_avg) %>% 
    rename(z2 = z2_avg) %>% 
    rename(z3 = z3_avg) %>%
    rename(f0 = f0_avg) %>%  #had to rename columns bc i rerun the function
    rename(f1 = f1_avg) %>% 
    rename(f2 = f2_avg) %>% 
    rename(f3 = f3_avg) 

  # convert the data frame to a tidy format
dfealpha_tidy <- gather(dfealpha_summary, 
    key = "generation", value = "value", c(z0,z1,z2,z3)) %>%
    mutate(generation = recode(generation,
    z0 = 'f0', 
    z1 = 'f1', 
    z2 = 'f2', 
    z3 = 'f3')) %>%
    mutate(selfing = as.character(selfing))

#_________________________________________________________________________________

#read in all selfing %, dfes, and experiments
grapes_dir <- "/nas/longleaf/home/adaigle/DFESelfing/grapes/"
output_dirs <- paste(grapes_dir, dir(grapes_dir, pattern = "output_\\d"), "/all", sep = "")

#table of results with columns signifying selfing%, DFE, and output
grapes_raw_results <- tibble(
    name = list.files(path = output_dirs, pattern = ".csv$"),
    matchname = sub(".txt.csv","",name), #will be useful for comparing results later, or combining with other programs
    DFE = str_extract(matchname, "(DFE)\\d+"),
    fullpath = paste(output_dirs, "/", name,sep=""),
    selfing = str_extract(fullpath, "(?<=grapes_output_)\\d+"),
    data = lapply(fullpath,read.csv)
)

#now I want to clean this up because it has too much going on
# I will keep just the gammazero row in each DF
grapes_gammazero_raw_results <- grapes_raw_results %>% 
    mutate(data = map(data, ~ { #start of a mapping funciton, which is applied to each df in column
    df_filtered <- .x %>% filter(model == 'GammaZero') %>% #select only GammaZero row
    select(1:52) %>%# removing prediction columns for now
    select_if(~ all(!is.na(.))) # remove NA columns
    return(df_filtered)
}))

#now i summarize the outputs by averaging them and finding the standard deviation
#first I turn my table of dataframes into a dataframe, because my df's only had one row 
# Extract the column of data frames and bind them into a single data frame
df <- grapes_gammazero_raw_results$data %>% bind_rows()

# Add a column to the data frame that contains the row numbers
df <- df %>% mutate(row_number = row_number())
grapes_gammazero_raw_results <- grapes_gammazero_raw_results  %>% 
    mutate(row_number = row_number())

# Join the original table with the new data frame using the row numbers as the key
grapes_gammazero_raw_results <- grapes_gammazero_raw_results %>% 
    left_join(df, by = "row_number") %>% 
    select(-data, -row_number) %>%
    as.data.frame

# summarize using group_by and summarize
# I originally tried using table of dataframes, but got weird
# for now I'll take one experiment at a time

true_gammas <- list(
  DFE1 = 5,
  DFE2 = 50,
  DFE3 = 1000
)
true_betas <- list(
  DFE1 = 0.9,
  DFE2 = 0.5,
  DFE3 = 0.3
)

grapes_raw_results_wtruth <- grapes_gammazero_raw_results %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE]))

#now make dfe class prediction for each gmean(gamma) and shape(b)
DFE_proportions_grapes <- function(meanGamma,beta) { 
    # modified function to calc vector of dfe classes given gamma and beta
    # assumes Nw is 100 bc grapes doesn't have two step
    # need to confirm it has no modified Ne
    Nw <- 100 
    meanS <- meanGamma/(2*Nw)

    # this code adapted from get_dfe_class_proportions
    s_shape <- beta
    s_rate <- s_shape/abs(meanS)
    x1 <- 1/(2*Nw)
    x10 <- 10/(2*Nw)
    x100 <- 100/(2*Nw)
    f0 <- pgamma(x1, shape=s_shape, rate=s_rate)
    f1 <- pgamma(x10, shape=s_shape, rate=s_rate) - f0
    f2 <- pgamma(x100, shape=s_shape, rate=s_rate) - f1 - f0
    f3 <- 1.0 - f2 - f1 - f0
    #print(paste(f0, f1, f2, f3, sep=", "))
    #return(c(list(f0), list(f1), list(f2), list(f3)))
    return(c(f0 = f0, f1 = f1, f2 = f2, f3 = f3))
}
grapes_raw_results_wclasses <- grapes_gammazero_raw_results %>%
  mutate(output = map2(GammaZero.negGmean, GammaZero.negGshape, DFE_proportions_grapes)) %>%
  unnest_wider(output) %>%
  rename(t0 = f0) %>%  #had to rename columns bc i rerun the function
  rename(t1 = f1) %>% 
  rename(t2 = f2) %>% 
  rename(t3 = f3) 
  #rename( c(t0,t1,t2,t3) = c(f0,f1,f2,f3)) %>%
  #mutate(output = map2(true_mean, true_shape, DFE_proportions_grapes)) %>%
  #unnest_wider(output) 

grapes_gammazero_summary <- grapes_raw_results_wclasses %>% 
    group_by(selfing, DFE) %>%
    summarize(across(where(is.numeric), list(avg = mean, sd = sd)))

grapes_gammazero_simple_summary <- 
    grapes_gammazero_summary[
        c('DFE','selfing','t0_avg','t1_avg', 't2_avg', 't3_avg',
            't0_sd','t1_sd', 't2_sd', 't3_sd')] %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE])) %>% 
    mutate(output = pmap(list(true_mean, true_shape, 100), DFE_proportions_dfe_alpha)) %>%
    unnest_wider(output) 

#### TRUE DATAFRAME
# I am making an extra df for true values to make things less confusing
# I will Have a mean and shape column, then run program to get classes, then make into a tidy df for easy plotting

dataframe_of_truth <- #must run calculate_pi.r
summary_table %>% 
    mutate(row_names = row.names(summary_table),
    DFE = str_extract(row_names, "(DFE)\\d+"),
    selfing = str_extract(row_names, "(?<=selfing)\\d+"),) %>% 
    tibble() %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE]))
dataframe_of_truth2 <-
  tidy_summary_table %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE]))
###making gamma and beta summary table for the non adaptive situation
dataframe_of_truth3 <- dataframe_of_truth2 %>%
  mutate(adjusted_gamma = B * true_mean) %>%
  mutate(matchname = paste0(DFE,"output",output)) %>%
  filter(!str_detect(fullpath, "pos")) %>%
  group_by(DFE,selfing,matchname) %>%
  select(c(matchname,B, empirical_Ne, DFE, selfing, true_mean, true_shape, adjusted_gamma)) 

dfealpha_summary2 <- dfealpha_raw_results %>%
  mutate(selfing = ifelse(selfing == "100", "99", selfing)) %>%
  group_by(DFE,selfing,matchname) %>%
  select(c(matchname,selfing,DFE,gamma,b)) 

grapes_gammazero_summmary2 <- grapes_gammazero_raw_results %>%
  group_by(DFE,selfing) %>%
  select(c(matchname,selfing,DFE,GammaZero.negGmean, GammaZero.negGshape)) 

prediction_accuracy_table <- dataframe_of_truth3 %>%
  left_join(dfealpha_summary2) %>%
  left_join(grapes_gammazero_summmary2) %>%
  arrange(DFE) %>%
  mutate (F = (as.numeric(selfing)/100)/(2-(as.numeric(selfing)/100)),
    selfing_Ne = 5000/(1+F),
    selfing_B = selfing_Ne / 5000,
    F_adjusted_gamma = selfing_B * true_mean,
    deltaNe = selfing_Ne - empirical_Ne,
    newNE = 5000 - deltaNe,
    newNegamma = (newNE/5000) * true_mean) %>% group_by(DFE,selfing)



## adding adjusted truths to plots
# will just append to end now so I don't break things
F_adjusted_classes <- prediction_accuracy_table %>%
  mutate(output = map2(F_adjusted_gamma, true_shape, DFE_proportions_grapes)) %>%
  unnest_wider(output) %>% 
  select(DFE, selfing, f0, f1, f2, f3) %>% 
  group_by(selfing, DFE) %>%
    summarize(across(where(is.numeric), list(avg = mean, sd = sd))) %>%
    rename(f0 = f0_avg) %>%  #had to rename columns bc i rerun the function
    rename(f1 = f1_avg) %>% 
    rename(f2 = f2_avg) %>% 
    rename(f3 = f3_avg) 
  
F_adjusted_classes_tidy <- gather(F_adjusted_classes, 
    key = "generation", value = "value", c(f0, f1, f2, f3)) %>%
    mutate(selfing = case_when(
        selfing == '0'  ~ 'F_adjusted_0', 
        selfing == '50' ~ 'F_adjusted_50', 
        selfing == '80' ~ 'F_adjusted_80', 
        selfing == '90' ~ 'F_adjusted_90',
        selfing == '95' ~ 'F_adjusted_95',
        selfing == '99' ~ 'F_adjusted_99',
        TRUE ~ selfing
    ))

DFE_proportions_truth <- function(B, meanGamma,beta) { 
    # modified function to calc vector of dfe classes given gamma and beta
    # assumes Nw is 100 bc grapes doesn't have two step
    # need to confirm it has no modified Ne
    Nw <- 100 # assuming this is what we should do? 
    meanS <- B * (meanGamma/(2*Nw))

    # this code adapted from get_dfe_class_proportions
    s_shape <- beta
    s_rate <- s_shape/abs(meanS)
    x1 <- 1/(2*Nw)
    x10 <- 10/(2*Nw)
    x100 <- 100/(2*Nw)
    f0 <- pgamma(x1, shape=s_shape, rate=s_rate)
    f1 <- pgamma(x10, shape=s_shape, rate=s_rate) - f0
    f2 <- pgamma(x100, shape=s_shape, rate=s_rate) - f1 - f0
    f3 <- 1.0 - f2 - f1 - f0
    #print(paste(f0, f1, f2, f3, sep=", "))
    #return(c(list(f0), list(f1), list(f2), list(f3)))
    return(c(f0 = f0, f1 = f1, f2 = f2, f3 = f3))
}

dataframe_of_truth <- dataframe_of_truth %>% filter(!str_detect(row_names, "pos")) %>%
mutate(output = pmap(list(B, true_mean, true_shape), DFE_proportions_truth)) %>%
    unnest_wider(output) 

dataframe_of_truth2_nopos <- dataframe_of_truth2 %>% filter(!str_detect(fullpath, "pos")) %>%
mutate(output = pmap(list(B, true_mean, true_shape), DFE_proportions_truth)) %>%
    unnest_wider(output) %>%
    group_by(selfing, DFE) %>%
    summarize(across(where(is.numeric), list(avg = mean, sd = sd))) %>%
    rename(f0 = f0_avg) %>%  #had to rename columns bc i rerun the function
    rename(f1 = f1_avg) %>% 
    rename(f2 = f2_avg) %>% 
    rename(f3 = f3_avg) 

truth_tidy <- gather(dataframe_of_truth2_nopos, 
    key = "generation", value = "value", c(f0, f1, f2, f3)) %>%
    mutate(selfing = case_when(
        selfing == '0'  ~ 'true0_recalc', 
        selfing == '50' ~ 'true50', 
        selfing == '80' ~ 'true80', 
        selfing == '90' ~ 'true90',
        selfing == '95' ~ 'true95',
        selfing == '99' ~ 'true99',
        TRUE ~ selfing
    ))

grapes_gammazero_simple_summary <- grapes_gammazero_simple_summary %>%
    rename(
        'z0' = f0,
        'z1'= f1 ,
        'z2'= f2 ,
        'z3'= f3 ,
       'f0' = t0_avg , 
       'f1' = t1_avg , 
       'f2' = t2_avg , 
        'f3' = t3_avg ,
        'f0_sd'= t0_sd , 
        'f1_sd'= t1_sd , 
        'f2_sd'= t2_sd , 
        'f3_sd' = t3_sd ) 
# convert the data frame to a tidy format
df_tidy <- gather(grapes_gammazero_simple_summary, 
    key = "generation", value = "value", c(f0, f1, f2, f3)) 
    #mutate(generation = recode(generation,
    #t0_avg = 'f0', 
    #t1_avg = 'f1', 
    #t2_avg = 'f2', 
    #t3_avg = 'f3')) 

df_true_tidy <- gather(grapes_gammazero_simple_summary, 
    key = "generation", value = "value", c(z0:z3)) %>%
    mutate(generation = recode(generation,
        z0 = 'f0', 
        z1 = 'f1', 
        z2 = 'f2', 
        z3 = 'f3'
    ))  %>% 
    select(1,13:14) %>% distinct() %>%
    mutate(selfing = "True0")

#grapes_for_plotting <- bind_rows(df_tidy, df_true_tidy) %>% 
#  mutate(generation = recode(generation,
#    t0_avg = 'f0', 
#    t1_avg = 'f1', 
#    t2_avg = 'f2', 
#    t3_avg = 'f3',
#    t0_sd = 'f0_sd', 
#    t1_sd = 'f1_sd', 
#    t2_sd = 'f2_sd', 
#    t3_sd = 'f3_sd')) %>% 
#    melt() %>%
#    filter(variable == "value" | paste0(generation, "_sd") == variable)

#grapes_for_plotting$variable <- ifelse(grepl("_sd", 
#    grapes_for_plotting$generation), "sd", "value")  

grapes_for_plotting <- bind_rows(df_true_tidy,df_tidy) %>%
    select(1:8) %>%
    melt() %>% 
    mutate(value = ifelse(is.na(value), 0, value)) %>% 
    #mutate(generation = recode(generation,
    #  t0_avg = 'f0', 
    #  t1_avg = 'f1', 
    #  t2_avg = 'f2', 
    #  t3_avg = 'f3',
    #  t0_sd = 'f0_sd', 
    #  t1_sd = 'f1_sd', 
    #  t2_sd = 'f2_sd', 
    #  t3_sd = 'f3_sd')) %>% 
    filter(variable == "value" | paste0(generation, "_sd") == variable)

grapes_for_plotting$variable <- ifelse(grepl("_sd", grapes_for_plotting$variable), "sd", "value")  


voodoo_grapes <- pivot_wider(grapes_for_plotting, id_cols = c("generation","DFE","selfing"), names_from = "variable", values_from = "value")  
#dplyr::group_by(generation, DFE, selfing, variable) %>%
#    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
#    dplyr::filter(n > 1L)
# create a bar plot with separate panels for each DFE
#ggplot(bind_rows(df_true_tidy,df_tidy) , aes(x = generation, y = value, fill = selfing)) +
#  geom_bar(stat = "identity", position = "dodge", colour = "black") +
#  facet_wrap(~ DFE, nrow = 2) +
#  labs(x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") 
#

# Define the order of the selfing levels
#selfing_order <- c("True0", 0, 50, 80, 90, 95, 99)
#
## Create the grouped bar chart with custom selfing order
#ggplot(bind_rows(df_true_tidy,df_tidy), aes(x = generation, y = value, fill = factor(selfing, levels = c("True0", "0", "50", "80", "90", "95", "99")))) +
#  geom_bar(stat = "identity", position = "dodge", colour = "black") +
#  facet_wrap(~ DFE, nrow = 2) +
#  labs(x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
#  scale_fill_manual(values = c("black", "red","blue","purple", "goldenrod4", "hotpink", "green"))

  
dfealpha_for_plotting <- bind_rows(df_true_tidy,dfealpha_tidy) %>%
    select(1:4,32,34,36,38) %>%
    melt() %>% 
    mutate(value = ifelse(is.na(value), 0, value))

# filter sd's
dfealpha_for_plotting <- dfealpha_for_plotting %>% 
  filter(variable == "value" | paste0(generation, "_sd") == variable)

dfealpha_for_plotting$variable <- ifelse(grepl("_sd", dfealpha_for_plotting$variable), "sd", "value")  

voodoo <- pivot_wider(dfealpha_for_plotting, id_cols = c("generation","DFE","selfing"), names_from = "variable", values_from = "value")
# Define the order of the selfing levels
#selfing_order <- c("True0", 0, 50, 80, 90, 95, 100)
selfing_order <- c("True0", 0, 50, 80, 90, 95, 100)

# Create the grouped bar chart with custom selfing order
ggplot(voodoo, aes(x = generation, y = value, fill = factor(selfing, levels = c("True0", "0", "50", "80", "90", "95", "100")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 2) +
  labs(title = "DFEalpha", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  #scale_fill_manual(values = c("black", "red","blue","purple", "goldenrod4", "hotpink", "green")) +
  scale_fill_manual(values = c("black", "red","blue","purple", "goldenrod4", "hotpink", "green")) +
  expand_limits(y=c(0,1))
#now I take the average, getting avg and sd for my classes
#I could also probably just ggplot but I think we will always 
#be using the averages for now


#grapes plot
selfing_order <- c("True0", 0, 50, 80, 90, 95, 99)

# Create the grouped bar chart with custom selfing order
ggplot(voodoo_grapes, aes(x = generation, y = value, fill = factor(selfing, levels = c("True0", "0", "50", "80", "90", "95", "99")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 2) +
  labs(title = "Grapes", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("black", "red","blue","purple", "goldenrod4", "hotpink", "green")) +
  expand_limits(y=c(0,1))
  

voodoo$name <- "DFEalpha"
voodoo_grapes$name <- "Grapes"

######stuff below here is using new truth
#### DFE ALPHA ONLY
dfealpha_for_plotting <- bind_rows(df_true_tidy,dfealpha_tidy,F_adjusted_classes_tidy,truth_tidy) %>%
    select(1:4,32,34,36,38) %>%
    melt() %>% 
    mutate(value = ifelse(is.na(value), 0, value)) %>% 
    filter(variable == "value" | paste0(generation, "_sd") == variable)
    
dfealpha_for_plotting$variable <- ifelse(grepl("_sd", dfealpha_for_plotting$variable), "sd", "value")  

voodoo <- pivot_wider(dfealpha_for_plotting, 
  id_cols = c("generation","DFE","selfing"), names_from = "variable", values_from = "value") %>%
mutate(selfing = case_when(
        selfing == 'True0'  ~ 'truth', 
        selfing == 'true0_recalc' ~ 'true0',
        selfing == '100' ~ '99',
        TRUE ~ selfing 
))
# Define the order of the selfing levels
#selfing_order <- c("True0", 0, 50, 80, 90, 95, 100)
#selfing_order <- c("True0", "true0_recalc", 0, "true50", 50, "true80", 80, "true90", 90, "true95", 95, "true99", 100)

#removing adjusted truths for rainbow plots for now
dfealpha_rainbow_plot <- voodoo %>% 
  filter(!grepl("F_adjusted|true", selfing)) 
  #filter(grepl("truth", selfing))
selfing_order <- c("truth", 0, 50, 80, 90, 95, 99)

ggplot(dfealpha_rainbow_plot, aes(x = generation, y = value, fill = factor(selfing, 
    levels = selfing_order))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 2) +
  labs(title = "DFEalpha", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))

#voodoo2 <- voodoo %>% mutate(selfing = case_when(
#        selfing == 'True0'  ~ 'truth', 
#        selfing == 'true0_recalc' ~ 'true0',
#        selfing == '100' ~ '99',
#        TRUE ~ selfing 
#))
#

# Split data frame into "truth" and "non-truth" data frames
voodoo2_truth <- voodoo %>% filter(grepl("truth", selfing))
voodoo2_not_truth <- voodoo %>% filter(!grepl("truth", selfing))

# Replicate "truth" data frame
df_truth_rep <- voodoo2_truth %>% slice(rep(1:n(), each = 6))

# Assign selfing class to replicated "truth" data frame
df_truth_rep_self <- df_truth_rep %>%
  mutate(
    replicate = rep(1:6, length.out = n()),
    selfing_class = case_when(
      replicate == 1 ~ "0% Selfing",
      replicate == 2 ~ "50% Selfing",
      replicate == 3 ~ "80% Selfing",
      replicate == 4 ~ "90% Selfing",
      replicate == 5 ~ "95% Selfing",
      replicate == 6 ~ "99% Selfing",
      TRUE ~ "error"
    )
  )

selfing_nums <- c(0, 50, 80, 90, 95, 99)
voodoo3 <- voodoo2_not_truth %>%
mutate(
    selfing_class = case_when(
      selfing == paste0("true", selfing_nums[1]) ~ "0% Selfing",
      selfing == "0" ~ "0% Selfing",
      selfing == "F_adjusted_0" ~ "0% Selfing",
      selfing == paste0("true", selfing_nums[2]) ~ "50% Selfing",
      selfing == "50" ~ "50% Selfing",
      selfing == "F_adjusted_50" ~ "50% Selfing",
      selfing == paste0("true", selfing_nums[3]) ~ "80% Selfing",
      selfing == "80" ~ "80% Selfing",
      selfing == "F_adjusted_80" ~ "80% Selfing",
      selfing == paste0("true", selfing_nums[4]) ~ "90% Selfing",
      selfing == "90" ~ "90% Selfing",
      selfing == "F_adjusted_90" ~ "90% Selfing",
      selfing == paste0("true", selfing_nums[5]) ~ "95% Selfing",
      selfing == "95" ~ "95% Selfing",
      selfing == "F_adjusted_95" ~ "95% Selfing",
      selfing == paste0("true", selfing_nums[6]) ~ "99% Selfing",
      selfing == "99" ~ "99% Selfing",
      selfing == "F_adjusted_99" ~ "99% Selfing",
      TRUE ~ "truth"
    )
  )

voodoo3 <- bind_rows(df_truth_rep_self, voodoo3)
simulated_truth <- voodoo3 %>%
  filter(grepl("true|F_adjusted", selfing))

selfing_order <- c("truth", "F_adjusted_0", "true0", "0", "F_adjusted_50","true50", "50", "F_adjusted_80","true80", "80", "F_adjusted_90","true90", "90", "F_adjusted_95","true95", "95", "F_adjusted_99","true99", "99")

ggplot(voodoo3, aes(x = generation, y = value, fill = factor(selfing, 
    levels = selfing_order))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "DFEalpha", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D"),6)))+ 
  theme(legend.position="none")

### grapes
  grapes_for_plotting <- bind_rows(df_true_tidy,df_tidy,F_adjusted_classes_tidy, truth_tidy) %>%
    select(1:8) %>%
    melt() %>% 
    mutate(value = ifelse(is.na(value), 0, value)) %>% 
    filter(variable == "value" | paste0(generation, "_sd") == variable)

grapes_for_plotting$variable <- ifelse(grepl("_sd", grapes_for_plotting$variable), "sd", "value")  


voodoo_grapes <- pivot_wider(grapes_for_plotting, id_cols = c("generation","DFE","selfing"), names_from = "variable", values_from = "value")  %>%
    mutate(selfing = case_when(
        selfing == "True0" ~ "truth",
        selfing == "true0_recalc" ~ "true0",
        TRUE ~ selfing
    ))
#test plot with all the truth
#selfing_order <- c("truth", "F_adjusted_0","true0", 0, "F_adjusted_50","true50", 50, "F_adjusted_80","true80", 80, "F_adjusted_90","true90", 90, "F_adjusted_95","true95", 95, "F_adjusted_99","true99", 99)

#now we are removing all the rainbows from the plots, so this code removes the adjusted truths 
#also cleaning up names for better legend
grapes_rainbow_plot <- voodoo_grapes %>%
  filter(!grepl("F_adjusted|true", selfing))
selfing_order <- c("truth", 0, 50, 80, 90, 95, 99)

# Create the grouped bar chart with custom selfing order
ggplot(grapes_rainbow_plot, aes(x = generation, y = value, fill = factor(selfing, 
    levels = selfing_order))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 1) +
  labs(title = "Grapes", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))


voodoo_grapes2 <- voodoo_grapes %>%
    mutate(selfing = case_when(
        selfing == 0~ "0_grapes",
        selfing == 50 ~ "50_grapes",
        selfing == 80 ~ "80_grapes",
        selfing == 90 ~ "90_grapes",
        selfing == 95 ~ "95_grapes",
        selfing == 99 ~ "99_grapes"
    )) %>%
    na.omit()


selfing_nums <- c(0, 50, 80, 90, 95, 99)
voodoo3_grapes <- voodoo_grapes2 %>%
mutate(
    selfing_class = case_when(
      selfing == "0_grapes" ~ "0% Selfing",
      selfing == "50_grapes" ~ "50% Selfing",
      selfing == "80_grapes" ~ "80% Selfing",
      selfing == "90_grapes" ~ "90% Selfing",
      selfing == "95_grapes" ~ "95% Selfing",
      selfing == "99_grapes" ~ "99% Selfing",
    )
  )

just_grapes_plot <- bind_rows(df_truth_rep_self, voodoo3_grapes, simulated_truth) %>%
filter(!grepl("F_adjusted|true", selfing))
#grapes only plot:
ggplot(just_grapes_plot, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("truth", "F_adjusted_0", "true0", "0_grapes", "F_adjusted_50", "true50", "50_grapes", "F_adjusted_80", "true80", "80_grapes", "F_adjusted_90", "true90", "90_grapes", "F_adjusted_95", "true95", "95_grapes", "F_adjusted_99", "true99", "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "Grapes", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values = c("#404040", rep(c("purple"),6)))+ 
  theme(legend.position="none")

F
# combo dfealpha and grapes plot:

combo_plot <- bind_rows(voodoo3,voodoo3_grapes) 
  #filter(grepl("truth|true", selfing))
  #filter(!grepl("F_adjusted", selfing))
ggplot(combo_plot, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("truth", "F_adjusted_0", "true0", 0, "0_grapes",
        "F_adjusted_50", "true50", 50, "50_grapes", "F_adjusted_80", "true80", 80, "80_grapes",
        "F_adjusted_90", "true90", 90, "90_grapes", "F_adjusted_95", "true95", 95, "95_grapes",
        "F_adjusted_99", "true99", 99, "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "DFEalpha and Grapes ", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  theme(legend.position="none", axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), plot.title= element_text(size=25))
##################################################
#Analyzing grapes GammaExpo results


grapes_dir <- "/nas/longleaf/home/adaigle/DFESelfing/grapes/"
pos_output_dirs <- paste(grapes_dir, dir(grapes_dir, pattern = "output_pos"), "/all", sep = "")

#table of results with columns signifying selfing%, DFE, and output
pos_grapes_raw_results <- tibble(
    name = list.files(path = pos_output_dirs, pattern = ".csv$"),
    matchname = sub(".txt.csv","",name), #will be useful for comparing results later, or combining with other programs
    DFE = str_extract(matchname, "(DFE)\\d+"),
    fullpath = paste(pos_output_dirs, "/", name,sep=""),
    selfing = str_extract(fullpath, "(?<=grapes_output_pos_)\\d+"),
    data = lapply(fullpath,read.csv)
)

#now I want to clean this up because it has too much going on
# I will keep just the gammazero row in each DF
pos_grapes_gammaexpo_raw_results <- pos_grapes_raw_results %>% 
    mutate(data = map(data, ~ { #start of a mapping funciton, which is applied to each df in column
    df_filtered <- .x %>% filter(model == 'GammaExpo') %>% #select only GammaZero row
    select(1:52) %>%# removing prediction columns for now
    select_if(~ all(!is.na(.))) # remove NA columns
    return(df_filtered)
}))

#now i summarize the outputs by averaging them and finding the standard deviation
#first I turn my table of dataframes into a dataframe, because my df's only had one row 
# Extract the column of data frames and bind them into a single data frame
df <- pos_grapes_gammaexpo_raw_results$data %>% bind_rows() %>% mutate(row_number = row_number())

pos_grapes_gammaexpo_raw_results <- pos_grapes_gammaexpo_raw_results  %>% 
    mutate(row_number = row_number()) %>% 
    left_join(df, by = "row_number") %>% 
    select(-data, -row_number) %>%
    as.data.frame

pos_grapes_gammaexpo_raw_results_wtruth <- pos_grapes_gammaexpo_raw_results %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE]))

pos_grapes_gammaexpo_raw_results_wclasses <- pos_grapes_gammaexpo_raw_results %>%
  mutate(output = map2(GammaExpo.negGmean, GammaExpo.negGshape, DFE_proportions_grapes)) %>%
  unnest_wider(output) %>%
  rename(t0 = f0) %>%  #had to rename columns bc i rerun the function
  rename(t1 = f1) %>% 
  rename(t2 = f2) %>% 
  rename(t3 = f3) %>% group_by(selfing,DFE)
  #rename( c(t0,t1,t2,t3) = c(f0,f1,f2,f3)) %>%
  #mutate(output = map2(true_mean, true_shape, DFE_proportions_grapes)) %>%
  #unnest_wider(output) 

alphas <- read.table( file = "/nas/longleaf/home/adaigle/DFESelfing/alpha_list.txt", header = TRUE)
#rename selfing levels to be compatable with graph
alphas <- alphas %>% mutate(selfing=recode(selfing, "0" = 'true_alpha_0',
       '50' ='true_alpha_50', 
        '80' = 'true_alpha_80', 
        '90' = 'true_alpha_90',
        '95' = 'true_alpha_95',
        '99' = 'true_alpha_99'
    ))%>% group_by(selfing, DFE)

pos_grapes_gammaexpo_raw_results_wclasses2<-rbind(pos_grapes_gammaexpo_raw_results_wclasses, alphas)


grapes_gammaexpo_summary <- pos_grapes_gammaexpo_raw_results_wclasses2 %>% 
group_by(selfing, DFE) %>%
summarize(across(where(is.numeric), list(avg = mean, sd = sd)))



grapes_gammaexpo_simple_summary <- 
    grapes_gammaexpo_summary[
        c('DFE','selfing','t0_avg','t1_avg', 't2_avg', 't3_avg',
            't0_sd','t1_sd', 't2_sd', 't3_sd')] %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE])) %>% 
    mutate(output = pmap(list(true_mean, true_shape, 100), DFE_proportions_dfe_alpha)) %>% #modified so I assume gamma = 2Nes
    unnest_wider(output) %>%
    rename(
        'z0' = f0,
        'z1'= f1 ,
        'z2'= f2 ,
        'z3'= f3 ,
        'f0' = t0_avg , 
        'f1' = t1_avg , 
        'f2' = t2_avg , 
        'f3' = t3_avg ,
        'f0_sd'= t0_sd , 
        'f1_sd'= t1_sd , 
        'f2_sd'= t2_sd , 
        'f3_sd' = t3_sd)
# convert the data frame to a tidy format
posgrapes_df_tidy <- gather(grapes_gammaexpo_simple_summary, 
    key = "generation", value = "value", c(f0, f1, f2, f3)) 
    #mutate(generation = recode(generation,
    #t0_avg = 'f0', 
    #t1_avg = 'f1', 
    #t2_avg = 'f2', 
    #t3_avg = 'f3')) 

posgrapes_df_true_tidy <- gather(grapes_gammaexpo_simple_summary, 
    key = "generation", value = "value", c(z0:z3)) %>%
    mutate(generation = recode(generation,
        z0 = 'f0', 
        z1 = 'f1', 
        z2 = 'f2', 
        z3 = 'f3'
    ))  %>% 
    select(1,13:14) %>% distinct() %>%
    mutate(selfing = "True0")

pos_dataframe_of_truth <- #must run calculate_pi.r
pos_summary_table %>% 
    mutate(row_names = row.names(pos_summary_table),
    DFE = str_extract(row_names, "(DFE)\\d+"),
    selfing = str_extract(row_names, "(?<=selfing)\\d+"),) %>% 
    tibble() %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE]))

pos_dataframe_of_truth <- pos_dataframe_of_truth %>% 
mutate(output = pmap(list(B, true_mean, true_shape), DFE_proportions_truth)) %>%
    unnest_wider(output) 

##new tidy version with multiple replicates. Can remove other two versions above if this is working out
dataframe_of_truth2_pos <- dataframe_of_truth2 %>% filter(str_detect(fullpath, "pos")) %>%
mutate(output = pmap(list(B, true_mean, true_shape), DFE_proportions_truth)) %>%
    unnest_wider(output) %>%
    group_by(selfing, DFE) %>%
    summarize(across(where(is.numeric), list(avg = mean, sd = sd))) %>%
    rename(f0 = f0_avg) %>%  #had to rename columns bc i rerun the function
    rename(f1 = f1_avg) %>% 
    rename(f2 = f2_avg) %>% 
    rename(f3 = f3_avg) 

pos_truth_tidy <- gather(dataframe_of_truth2_pos, 
    key = "generation", value = "value", c(f0, f1, f2, f3)) %>%
    mutate(selfing = case_when(
        selfing == '0'  ~ 'true0_recalc', 
        selfing == '50' ~ 'true50', 
        selfing == '80' ~ 'true80', 
        selfing == '90' ~ 'true90',
        selfing == '95' ~ 'true95',
        selfing == '99' ~ 'true99',
        TRUE ~ selfing
    ))

pos_grapes_for_plotting <- bind_rows(posgrapes_df_tidy,posgrapes_df_true_tidy,pos_truth_tidy) %>%
    select(1:6,13:14) %>%
    melt() %>% 
    mutate(value = ifelse(is.na(value), 0, value)) %>% 
    filter(variable == "value" | paste0(generation, "_sd") == variable)

pos_grapes_for_plotting$variable <- ifelse(grepl("_sd", pos_grapes_for_plotting$variable), "sd", "value")  

pos_voodoo_grapes <- pivot_wider(pos_grapes_for_plotting, id_cols = c("generation","DFE","selfing"), 
    names_from = "variable", values_from = "value") 
pos_voodoo_grapes <- pos_voodoo_grapes %>% 
    mutate(selfing = case_when(
        selfing == "True0" ~ "truth",
        selfing == "true0_recalc" ~ "true0",
        selfing == 0~ "0_grapes",
        selfing == 50 ~ "50_grapes",
        selfing == 80 ~ "80_grapes",
        selfing == 90 ~ "90_grapes",
        selfing == 95 ~ "95_grapes",
        selfing == 99 ~ "99_grapes",
        TRUE ~ selfing
    ))
pos_voodoo_grapes <- filter(pos_voodoo_grapes, !(selfing %in% c('true90', 'true80', 'true95'))) #leaving filtering for last so that I can toggle what I want in the chart
#pos_voodoo_grapes <- filter(pos_voodoo_grapes, !(DFE %in% c('DFE3'))) # was used before DFE2 run was complete
selfing_order <- c("truth", "true0", "0_grapes", "true50", "50_grapes", "true99", "99_grapes")

# Create the grouped bar chart with custom selfing order
ggplot(pos_voodoo_grapes, aes(x = generation, y = value, fill = factor(selfing, levels = c("truth", "true0", "0_grapes", "true50", "50_grapes", "true99", "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 2) +
  labs(title = "Grapes positive selection deleterious DFEs", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1))

############# summary tables
#write.csv(dfealpha_summary, file = "/nas/longleaf/home/adaigle/DFESelfing/dfealpha_summary_unfolded.csv")
#dfealpha_summary_simple <- select(dfealpha_summary, DFE, b_avg, Es_avg, N1_avg, N2_avg, t2_avg) %>%
#    rename_at(vars(ends_with("_avg")), ~ str_remove(.x, "_avg")) %>%
#    arrange(DFE) %>% rename(s = Es, beta = b, Selfing_rate = selfing)
#write.csv(dfealpha_summary_simple, row.names = FALSE,
#    file = "/nas/longleaf/home/adaigle/DFESelfing/dfealpha_summary_simple.csv")
#
#write.csv(grapes_gammaexpo_summary, file = "/nas/longleaf/home/adaigle/DFESelfing/grapes_gammaexpo_summary.csv")
#write.csv(grapes_gammazero_summary, file = "/nas/longleaf/home/adaigle/DFESelfing/grapes_gammazero_summary.csv")
#combination plot
#ggplot(pos_voodoo_grapes, aes(x = generation, y = value, fill = factor(selfing, 
#    levels = c("truth", "true0", 0, "0_grapes",
#        "true50", 50, "50_grapes", "true80", 80, "80_grapes",
#        "true90", 90, "90_grapes", "true95", 95, "95_grapes",
#        "true99", 99, "99_grapes")))) +
#  geom_bar(stat = "identity", position = "dodge", colour = "black") +
#  labs(title = "Grapes positive selection deleterious DFEs", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
#  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
#  expand_limits(y=c(0,1)) +
#  facet_grid(rows = vars(DFE), cols = vars(case_when(
#    as.numeric(gsub("[^0-9]", "", as.character(selfing))) <= 0 ~ "0% Selfing",
#    as.numeric(gsub("[^0-9]", "", as.character(selfing))) <= 50 ~ "50% Selfing",
#    as.numeric(gsub("[^0-9]", "", as.character(selfing))) <= 80 ~ "80% Selfing",
#    as.numeric(gsub("[^0-9]", "", as.character(selfing))) <= 90 ~ "90% Selfing",
#    as.numeric(gsub("[^0-9]", "", as.character(selfing))) <= 95 ~ "95% Selfing",
#    as.numeric(gsub("[^0-9]", "", as.character(selfing))) <= 99 ~ "99% Selfing",
#    TRUE ~ "uncorrected truth"
#  )))

pos_voodoo_grapes2_truth <- pos_voodoo_grapes %>% filter(grepl("truth", selfing))
pos_voodoo_grapes2_not_truth <- pos_voodoo_grapes %>% filter(!grepl("truth", selfing))

# Replicate "truth" data frame
df_truth_rep <- pos_voodoo_grapes2_truth %>% slice(rep(1:n(), each = 6))

# Assign selfing class to replicated "truth" data frame
df_truth_rep_self <- df_truth_rep %>%
  mutate(
    replicate = rep(1:6, length.out = n()),
    selfing_class = case_when(
      replicate == 1 ~ "0% Selfing",
      replicate == 2 ~ "50% Selfing",
      replicate == 3 ~ "80% Selfing",
      replicate == 4 ~ "90% Selfing",
      replicate == 5 ~ "95% Selfing",
      replicate == 6 ~ "99% Selfing",
      TRUE ~ "error"
    )
  ) %>%
  filter(!(replicate %in% c(3, 4, 5))) %>% #leaving this here in case we add more selfing levels
  select(-replicate)
  


pos_voodoo_grapes2 <- pos_voodoo_grapes2_not_truth %>%
mutate(
    selfing_class = case_when(
      selfing == "0_grapes" ~ "0% Selfing",
      selfing == "50_grapes" ~ "50% Selfing",
      #selfing == "80_grapes" ~ "80% Selfing",
      #selfing == "90_grapes" ~ "90% Selfing",
      #selfing == "95_grapes" ~ "95% Selfing",
      selfing == "99_grapes" ~ "99% Selfing",
      selfing == paste0("true", selfing_nums[1]) ~ "0% Selfing",
      selfing == paste0("true", selfing_nums[2]) ~ "50% Selfing",
      #selfing == paste0("true", selfing_nums[3]) ~ "80% Selfing",
      #selfing == paste0("true", selfing_nums[4]) ~ "90% Selfing",
      #selfing == paste0("true", selfing_nums[5]) ~ "95% Selfing",
      selfing == paste0("true", selfing_nums[6]) ~ "99% Selfing",
    )
  ) 


#plot gammaexpo results
ggplot(rbind(pos_voodoo_grapes2, df_truth_rep_self), aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("truth", "true0", "0_grapes", "true50", "50_grapes", "true80", "80_grapes", "true90", "90_grapes", "true95", "95_grapes", "true99", "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "Grapes (with beneficial mutations)", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values = c("#404040", rep(c("#619CFF", "purple"),6))) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))


#### plotting adv muts using grapes_gammaexpo_summary
# dividing mean by 2? I am assuming it is parameterized by 4Nes and I want it to match our 2NeS
ggplot(grapes_gammaexpo_summary, aes(x = DFE, y = GammaExpo.posGmean_avg/2, fill = selfing)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = GammaExpo.posGmean_avg/2 - GammaExpo.posGmean_sd/2,
                    ymax = GammaExpo.posGmean_avg/2 + GammaExpo.posGmean_sd/2),
                position = position_dodge(width = 0.9), width = 0.2) +
  labs(x = "DFE", y = "GammaExpo.posGmean_avg")

ggplot(grapes_gammaexpo_summary, aes(x = DFE, y = GammaExpo.pos_prop_avg, fill = selfing)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = GammaExpo.pos_prop_avg - GammaExpo.pos_prop_sd,
                    ymax = GammaExpo.pos_prop_avg + GammaExpo.pos_prop_sd),
                position = position_dodge(width = 0.9), width = 0.2) +
  labs(x = "DFE", y = "GammaExpo.pos_prop_avg")

#alpha estimates gamma
alphaplot <- grapes_gammaexpo_summary %>% 
    mutate(selfing = case_when(
        selfing == 'true_alpha_0'  ~ 'True alpha, 0% Selfing', 
        selfing == 'true_alpha_50' ~ 'True alpha, 50% Selfing', 
        selfing == 'true_alpha_99' ~ 'True alpha, 99% Selfing',
        selfing == '0' ~ 'Estimated alpha, 0% Selfing', 
        selfing == '50' ~ 'Estimated alpha, 50% Selfing', 
        selfing == '99' ~ 'Estimated alpha, 99% Selfing', 
        TRUE ~ selfing
    ))

alpha_order <- c('True alpha, 0% Selfing', 'Estimated alpha, 0% Selfing', 
  'True alpha, 50% Selfing', 'Estimated alpha, 50% Selfing', 
  'True alpha, 99% Selfing', 'Estimated alpha, 99% Selfing')
ggplot(alphaplot, aes(x = DFE, y = alpha_avg, fill = factor(selfing, alpha_order))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = alpha_avg - alpha_sd,
                    ymax = alpha_avg + alpha_sd),
                position = position_dodge(width = 0.9), width = 0.2) +
  labs(title = "Grapes GammaExpo model, alpha prediction", x = "DFE", y = "alpha", fill = "") + 
  theme(axis.text.x = element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15),legend.text = element_text(size=15))

  #############
  #plots for Galtier
#noselfingplot <- rbind(
#        filter(bind_rows(voodoo2, voodoo_grapes2), selfing == 'true0'),
#        filter(bind_rows(voodoo2, voodoo_grapes2), selfing == '0_grapes')
#        ) %>%
#    mutate(selfing = recode(selfing,
#        'true0' = 'True DFE',
#        '0_grapes' = 'Grapes estimated DFE'
#    ))
#
#ggplot(noselfingplot, aes(x = generation, y = value, fill = factor(selfing, 
#    levels = c("True DFE", 'Grapes estimated DFE')))) +
#  geom_bar(stat = "identity", position = "dodge", colour = "black") +
#  labs(title = "Grapes estimation vs truth value (Using gamma = NeS)", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
#  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
#  expand_limits(y=c(0,1)) +
#  facet_grid(rows = vars(case_when(
#    as.numeric(gsub("[^0-9]", "", as.character(selfing))) <= 0 ~ "",
#  )), cols = vars(DFE))


###making gamma and beta summary table for the non adaptive situation
dataframe_of_truth3 <- dataframe_of_truth2 %>%
  mutate(adjusted_gamma = B * true_mean) %>%
  mutate(matchname = paste0(DFE,"output",output)) %>%
  filter(!str_detect(fullpath, "pos")) %>%
  group_by(DFE,selfing,matchname) %>%
  select(c(matchname,B, empirical_Ne, DFE, selfing, true_mean, true_shape, adjusted_gamma)) 

dfealpha_summary2 <- dfealpha_raw_results %>%
  #mutate(selfing = ifelse(selfing == "100", "99", selfing)) %>%
  group_by(DFE,selfing,matchname) %>%
  select(c(matchname,selfing,DFE,gamma,b)) 

grapes_gammazero_summmary2 <- grapes_gammazero_raw_results %>%
  group_by(DFE,selfing) %>%
  select(c(matchname,selfing,DFE,GammaZero.negGmean, GammaZero.negGshape)) 

prediction_accuracy_table <- dataframe_of_truth3 %>%
  left_join(dfealpha_summary2) %>%
  left_join(grapes_gammazero_summmary2) %>%
  arrange(DFE) %>%
  mutate (F = (as.numeric(selfing)/100)/(2-(as.numeric(selfing)/100)),
    selfing_Ne = 5000/(1+F),
    selfing_B = selfing_Ne / 5000,
    F_adjusted_gamma = selfing_B * true_mean,
    deltaNe = selfing_Ne - empirical_Ne,
    newNE = 5000 - deltaNe,
    newNegamma = (newNE/5000) * true_mean) %>% group_by(DFE,selfing)



## adding adjusted truths to plots
# will just append to end now so I don't break things
F_adjusted_classes <- prediction_accuracy_table %>%
  mutate(output = map2(F_adjusted_gamma, true_shape, DFE_proportions_grapes)) %>%
  unnest_wider(output) %>% 
  select(DFE, selfing, f0, f1, f2, f3) %>% 
  group_by(selfing, DFE) %>%
    summarize(across(where(is.numeric), list(avg = mean, sd = sd))) %>%
    rename(f0 = f0_avg) %>%  #had to rename columns bc i rerun the function
    rename(f1 = f1_avg) %>% 
    rename(f2 = f2_avg) %>% 
    rename(f3 = f3_avg) 
  
F_adjusted_classes_tidy <- gather(F_adjusted_classes, 
    key = "generation", value = "value", c(f0, f1, f2, f3)) %>%
    mutate(selfing = case_when(
        selfing == '0'  ~ 'F_adjusted_0', 
        selfing == '50' ~ 'F_adjusted_50', 
        selfing == '80' ~ 'F_adjusted_80', 
        selfing == '90' ~ 'F_adjusted_90',
        selfing == '95' ~ 'F_adjusted_95',
        selfing == '99' ~ 'F_adjusted_99',
        TRUE ~ selfing
    ))


prediction_accuracy_table_gamma_melt <- prediction_accuracy_table %>% 
  summarize(across(c(true_mean,adjusted_gamma,F_adjusted_gamma, gamma, GammaZero.negGmean), list(avg = mean, sd = sd))) %>%
  rename("true gamma_avg" = true_mean_avg,
    "BGS-adjusted gamma_avg" = adjusted_gamma_avg,
    "F-adjusted gamma_avg" = F_adjusted_gamma_avg,
    "DFEalpha gamma_avg" = gamma_avg,
    "Grapes gamma_avg" = GammaZero.negGmean_avg,
    "true gamma_sd" = true_mean_sd,
    "BGS-adjusted gamma_sd" = adjusted_gamma_sd,
    "F-adjusted gamma_sd" = F_adjusted_gamma_sd,
    "DFEalpha gamma_sd" = gamma_sd,
    "Grapes gamma_sd" = GammaZero.negGmean_sd) %>% 
  pivot_longer(
    cols = c(3:12), 
    names_to = c(".value", "variable"), 
    names_sep = "_", 
    values_drop_na = TRUE
  ) %>%

  melt(variable.name = "gamma_method", id.vars = c("DFE", "selfing", "variable")) %>%
    pivot_wider(id_cols = c("DFE", "selfing", "gamma_method"), 
                        names_from = "variable", 
                        values_from = "value")

#making selfing class column to make graph labels prettier
gammaorder <- c("true gamma", "F-adjusted gamma", "BGS-adjusted gamma", "DFEalpha gamma", "Grapes gamma")

prediction_accuracy_table_gamma_melt <- prediction_accuracy_table_gamma_melt %>%
mutate(
    selfing_class = case_when(
      selfing == "0" ~ "0% Selfing",
      selfing == "50" ~ "50% Selfing",
      selfing == "80" ~ "80% Selfing",
      selfing == "90" ~ "90% Selfing",
      selfing == "95" ~ "95% Selfing",
      selfing == "99" ~ "99% Selfing"
    )
  ) %>%
  arrange(factor(gamma_method, levels = gammaorder))

#colorblind palete?
cbPalette <- c("grey", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
ggplot(prediction_accuracy_table_gamma_melt, aes(x = factor(gamma_method, levels = gammaorder), y = avg, fill = factor(gamma_method, levels = gammaorder))) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class), scales = "free_y") +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  #stat_summary(fun.data = "mean_sdl", geom = "bar", position = "dodge") +
  geom_errorbar(aes(ymin = avg-sd, ymax = avg+sd), position = position_dodge(width = 0.9), width = 0.2) +
  labs(x = "Gamma estimation method", y = "Gamma", fill = "Gamma estimation method", title = "Gamma estimation accuracy") + 
  scale_fill_manual(values=c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  theme(axis.text.x = element_blank(), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15),legend.text = element_text(size=15))

ggplot(bind_rows(voodoo3,voodoo3_grapes), aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("truth", "F_adjusted_0", "true0", 0, "0_grapes",
        "F_adjusted_50", "true50", 50, "50_grapes", "F_adjusted_80", "true80", 80, "80_grapes",
        "F_adjusted_90", "true90", 90, "90_grapes", "F_adjusted_95", "true95", 95, "95_grapes",
        "F_adjusted_99", "true99", 99, "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "DFEalpha and Grapes ", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6)))+ 
  theme(legend.position="none")
#beta stuff
prediction_accuracy_table_beta_melt <- prediction_accuracy_table %>% 
  summarize(across(c(true_shape, b, GammaZero.negGshape), list(avg = mean, sd = sd))) %>%
  rename("true shape_avg" = true_shape_avg,
    "DFEalpha beta_avg" = b_avg,
    "Grapes beta_avg" = GammaZero.negGshape_avg,
    "true shape_sd" = true_shape_sd,
    "DFEalpha beta_sd" = b_sd,
    "Grapes beta_sd" = GammaZero.negGshape_sd) %>%
  pivot_longer(
    cols = c(3:8), 
    names_to = c(".value", "variable"), 
    names_sep = "_", 
    values_drop_na = TRUE
  ) %>%
   melt(variable.name = "gamma_method", id.vars = c("DFE", "selfing", "variable")) %>%
    pivot_wider(id_cols = c("DFE", "selfing", "gamma_method"), 
                        names_from = "variable", 
                        values_from = "value") %>%
mutate(
    selfing_class = case_when(
      selfing == "0" ~ "0% Selfing",
      selfing == "50" ~ "50% Selfing",
      selfing == "80" ~ "80% Selfing",
      selfing == "90" ~ "90% Selfing",
      selfing == "95" ~ "95% Selfing",
      selfing == "99" ~ "99% Selfing"
    )
  ) 

ggplot(prediction_accuracy_table_beta_melt, aes(x = gamma_method, y = avg, fill = gamma_method)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class), scales = "free_y") +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  #stat_summary(fun.data = "mean_sdl", geom = "bar", position = "dodge") +
  geom_errorbar(aes(ymin = avg-sd, ymax = avg+sd), position = position_dodge(width = 0.9), width = 0.2) +
  labs(title = "Beta estimation accuracy", x = "Method", y = "Beta", fill = "") +
  scale_fill_manual(values=c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(axis.text.x = element_blank(), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15),legend.text = element_text(size=15))


#t_tests <- prediction_accuracy_table %>% group_by(DFE,selfing) %>%
#  summarize(dfealpha_v_true = list(t.test(gamma, true_mean, paired = TRUE)),
#    grapes_v_true = list(t.test(GammaZero.negGmean, true_mean, paired = TRUE)),
#    dfealpha_v_adjustedgamma = list(t.test(gamma, adjusted_gamma, paired = TRUE)),
#    grapes_v_adjustedgamma = list(t.test(GammaZero.negGmean, adjusted_gamma, paired = TRUE)),
#    dfealpha_v_Fgamma = list(t.test(gamma, F_adjusted_gamma, paired = TRUE)),
#    grapes_v_Fgamma = list(t.test(GammaZero.negGmean, F_adjusted_gamma, paired = TRUE)),
#    dfealpha_newNegamma = list(t.test(gamma, newNegamma, paired = TRUE)),
#    grapes_v_newNegamma = list(t.test(GammaZero.negGmean, newNegamma, paired = TRUE))
#  )
#t_tests_pval <- prediction_accuracy_table %>% 
#  group_by(DFE,selfing) %>%
#  summarize(dfealpha_v_true = t.test(gamma, true_mean, paired = TRUE)$p.value,
#            grapes_v_true = t.test(GammaZero.negGmean, true_mean, paired = TRUE)$p.value,
#            dfealpha_v_adjustedgamma = t.test(gamma, adjusted_gamma, paired = TRUE)$p.value,
#            grapes_v_adjustedgamma = t.test(GammaZero.negGmean, adjusted_gamma, paired = TRUE)$p.value,
#            dfealpha_v_Fgamma = t.test(gamma, F_adjusted_gamma, paired = TRUE)$p.value,
#            grapes_v_Fgamma = t.test(GammaZero.negGmean, F_adjusted_gamma, paired = TRUE)$p.value,
#            dfealpha_newNegamma = t.test(gamma, newNegamma, paired = TRUE)$p.value,
#            grapes_v_newNegamma = t.test(GammaZero.negGmean, newNegamma, paired = TRUE)$p.value,
#            across(where(is.numeric), list(avg = mean, sd = sd)))
#
## calc true alpha and plot with est_alpha
#count_sfs_names <- list.files(list.dirs(path = "/nas/longleaf/home/adaigle/SFS/pos"), pattern = "count.fixed")
#count_sfs_df_list <- lapply(
#    paste0(path_to_files,count_sfs_names), 
#    function(x) read.table(x, header=TRUE))
#
#fixed_table <- tibble(
#    name = list.files(path = path_to_files, pattern = "count.fixed"),
#    matchname = sub("_count.fixed","",name),
#    fullpath = paste0(path_to_files, "/", name), 
#    data = lapply(fullpath,read.table)
#)

## KS test
#set.seed(42)
#n <- 100
#test <- prediction_accuracy_table %>% mutate(
#  true_rgamma = list(rgamma(n,shape=true_shape, rate = 1/(true_mean/(2*100)))),
#  dfealpha_rgamma = list(rgamma(n,shape=as.numeric(b), rate = 1/(as.numeric(gamma)/(2*100)))),
#  grapes_rgamma = list(rgamma(n,shape=GammaZero.negGshape, rate = 1/(GammaZero.negGmean/(2*100)))),
#  adusted_rgamma = list(rgamma(n,shape=true_shape, rate = 1/(adjusted_gamma/(2*100)))),
#  F_rgamma = list(rgamma(n,shape=true_shape, rate = 1/(F_adjusted_gamma/(2*100)))),
#  newNe_rgamma = list(rgamma(n,shape=true_shape, rate = 1/(newNegamma/(2*100)))),
#) 
#func <- function(x,y) { ks.test(x,y)$p.value}
#mapply(func, test$true_rgamma,test$dfealpha_rgamma)
#mapply(func, test$true_rgamma,test$dfealpha_rgamma)
#mapply(func, test$true_rgamma,test$dfealpha_rgamma)
#mapply(func, test$true_rgamma,test$dfealpha_rgamma)
#
#ks_tests_pval <- test %>% 
#  group_by(DFE,selfing) %>%
#  summarize(dfealpha_v_true = list(mapply(func, dfealpha_rgamma, true_rgamma)),
#            grapes_v_true = list(mapply(func, grapes_rgamma, true_rgamma)),
#            dfealpha_v_adjustedgamma = list(mapply(func, dfealpha_rgamma, adusted_rgamma)),
#            grapes_v_adjustedgamma = list(mapply(func, grapes_rgamma, adusted_rgamma)),
#            dfealpha_v_Fgamma = list(mapply(func, dfealpha_rgamma, F_rgamma)),
#            grapes_v_Fgamma = list(mapply(func, grapes_rgamma, F_rgamma)),
#            dfealpha_newNegamma = list(mapply(func, dfealpha_rgamma, newNe_rgamma)),
#            grapes_v_newNegamma = list(mapply(func, grapes_rgamma, newNe_rgamma))
#  )
#            #across(where(is.numeric), list(avg = mean, sd = sd)))
#
#
#prediction_accuracy_table %>% filter(DFE == "DFE1", selfing == 99)


#gamma plot
#hist(rgamma(n=100000, shape = 100, rate = (100/.025)), type = 'line')

##################################################
#Analyzing grapes GammaExpo results


dfealpha_dom_dir <- "/nas/longleaf/home/adaigle/work/dominance_inputsandoutputs/"
dfealpha_dom_output_dirs <- paste(paste(dfealpha_dom_dir, dir(dfealpha_dom_dir, pattern = "DFE_alpha_dom_output"), sep = ""),
    "/", dir(
    paste(dfealpha_dom_dir, dir(dfealpha_dom_dir, pattern = "DFE_alpha_dom_output"), sep = ""), pattern = "selected")
    , sep = "")


#table of results with columns signifying selfing%, DFE, and output
dfealpha_dom_colnames <- as.data.frame(c("N1", "N2", "t2", "Nw", "b", "Es", "f0","L"))
dfealpha_dom_raw_results <- tibble(
    name = list.files(path = dfealpha_dom_output_dirs, pattern = "est_dfe.out"),
    fullpath = paste(dfealpha_dom_output_dirs, "/", name,sep=""),
    selfing = str_extract(fullpath, "(?<=DFE_alpha_dom_output_)\\d+") %>% if_else(is.na(.), "0", .), #0's were empty
    matchname = str_extract(fullpath, "(DFE)\\d+(output)\\d"), #will be useful for comparing results later, or combining with other programs
    DFE = str_extract(matchname, "(DFE)\\d+"),
    data = lapply(fullpath, function(x) 
        as.data.frame(t(data.frame(matrix(unlist(strsplit(readLines(x), " ")), ncol = 2, byrow = TRUE)))) %>%
            rownames_to_column() %>%
            `colnames<-`(.[1,]) %>%
            .[-1,] %>%
            `rownames<-`(NULL))
)

df <- dfealpha_dom_raw_results$data %>% bind_rows()
df <- df %>% mutate(row_number = row_number())
dfealpha_dom_raw_results <- dfealpha_dom_raw_results  %>% 
    mutate(row_number = row_number())

# Join the original table with the new data frame using the row numbers as the key
# GAMMA IS $NWeS because of h=0.25 and not .5
dfealpha_dom_raw_results <- dfealpha_dom_raw_results %>% 
    left_join(df, by = "row_number") %>% 
    select(-data, -row_number) %>%
    as.data.frame %>% 
    mutate_at(vars(7:14), as.numeric) %>%
    mutate(gamma = -4*Nw*Es) %>%
    rename(f0_fromoutput = f0) #f0 is a dfealpha_dom output, but I also need to name a class f0 later

DFE_proportions_dfe_alpha_dom <- function(meanGamma,beta, Nw) { 
    # modified function to calc vector of dfe classes given gamma and beta
    # assumes Nw is 100 bc grapes doesn't have two step
    # need to confirm it has no modified Ne
    meanS <- meanGamma/(2*Nw)

    # this code adapted from get_dfe_class_proportions
    s_shape <- beta
    s_rate <- s_shape/abs(meanS)
    x1 <- 1/(2*Nw)
    x10 <- 10/(2*Nw)
    x100 <- 100/(2*Nw)
    f0 <- pgamma(x1, shape=s_shape, rate=s_rate)
    f1 <- pgamma(x10, shape=s_shape, rate=s_rate) - f0
    f2 <- pgamma(x100, shape=s_shape, rate=s_rate) - f1 - f0
    f3 <- 1.0 - f2 - f1 - f0
    #print(paste(f0, f1, f2, f3, sep=", "))
    #return(c(list(f0), list(f1), list(f2), list(f3)))
    return(c(f0 = f0, f1 = f1, f2 = f2, f3 = f3))
} 
true_gammas <- list(
  DFE1 = 5,
  DFE2 = 50,
  DFE3 = 1000
)
true_betas <- list(
  DFE1 = 0.9,
  DFE2 = 0.5,
  DFE3 = 0.3
)
dfealpha_dom_raw_results_wtruth <- dfealpha_dom_raw_results %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE]))

#now i run the class generator function on all outputs
dfealpha_dom_raw_results_wclasses <- dfealpha_dom_raw_results_wtruth %>%
  mutate(output = pmap(list(gamma, b, Nw), DFE_proportions_dfe_alpha_dom)) %>%
  unnest_wider(output) %>%
  rename(z0 = f0) %>%  #had to rename columns bc i rerun the function
  rename(z1 = f1) %>% 
  rename(z2 = f2) %>% 
  rename(z3 = f3) %>% 
  mutate(output = pmap(list(gamma, b, Nw), DFE_proportions_dfe_alpha_dom)) %>%
  unnest_wider(output) 

dfealpha_dom_summary <- dfealpha_dom_raw_results_wclasses %>% 
    group_by(selfing, DFE) %>%
    summarize(across(where(is.numeric), list(avg = mean, sd = sd))) %>%
    rename(z0 = z0_avg) %>%  #had to rename columns bc i rerun the function
    rename(z1 = z1_avg) %>% 
    rename(z2 = z2_avg) %>% 
    rename(z3 = z3_avg) %>%
    rename(f0 = f0_avg) %>%  #had to rename columns bc i rerun the function
    rename(f1 = f1_avg) %>% 
    rename(f2 = f2_avg) %>% 
    rename(f3 = f3_avg) 

  # convert the data frame to a tidy format
dfealpha_dom_tidy <- gather(dfealpha_dom_summary, 
    key = "generation", value = "value", c(z0,z1,z2,z3)) %>%
    mutate(generation = recode(generation,
    z0 = 'f0', 
    z1 = 'f1', 
    z2 = 'f2', 
    z3 = 'f3')) %>%
    mutate(selfing = as.character(selfing))

grapes_dir <- "/nas/longleaf/home/adaigle/work/dominance_inputsandoutputs/"
dom_output_dirs <- paste0(grapes_dir, dir(grapes_dir, pattern = "grapes_dom_output"), "/all")

#table of results with columns signifying selfing%, DFE, and output
dom_grapes_raw_results <- tibble(
    name = list.files(path = dom_output_dirs, pattern = ".csv$"),
    matchname = sub(".txt.csv","",name), #will be useful for comparing results later, or combining with other programs
    DFE = str_extract(matchname, "(DFE)\\d+"),
    fullpath = paste(dom_output_dirs, "/", name,sep=""),
    selfing = str_extract(fullpath, "(?<=grapes_dom_output_)\\d+"),
    data = lapply(fullpath,read.csv)
)

#now I want to clean this up because it has too much going on
# I will keep just the gammazero row in each DF
dom_grapes_gammazero_raw_results <- dom_grapes_raw_results %>% 
    mutate(data = map(data, ~ { #start of a mapping funciton, which is applied to each df in column
    df_filtered <- .x %>% filter(model == 'GammaZero') %>% #select only GammaZero row
    select(1:52) %>%# removing prediction columns for now
    select_if(~ all(!is.na(.))) # remove NA columns
    return(df_filtered)
}))

#now i summarize the outputs by averaging them and finding the standard deviation
#first I turn my table of dataframes into a dataframe, because my df's only had one row 
# Extract the column of data frames and bind them into a single data frame
df <- dom_grapes_gammazero_raw_results$data %>% bind_rows() %>% mutate(row_number = row_number())

dom_grapes_gammazero_raw_results <- dom_grapes_gammazero_raw_results  %>% 
    mutate(row_number = row_number()) %>% 
    left_join(df, by = "row_number") %>% 
    select(-data, -row_number) %>%
    as.data.frame

dom_grapes_gammazero_raw_results_wtruth <- dom_grapes_gammazero_raw_results %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE]))

DFE_proportions_grapes_dom <- function(meanGamma,beta) { 
    # modified function to calc vector of dfe classes given gamma and beta
    # assumes Nw is 100 bc grapes doesn't have two step
    # need to confirm it has no modified Ne
    Nw <- 100 
    meanS <- meanGamma/(2*Nw)

    # this code adapted from get_dfe_class_proportions
    s_shape <- beta
    s_rate <- s_shape/abs(meanS)
    x1 <- 1/(2*Nw)
    x10 <- 10/(2*Nw)
    x100 <- 100/(2*Nw)
    f0 <- pgamma(x1, shape=s_shape, rate=s_rate)
    f1 <- pgamma(x10, shape=s_shape, rate=s_rate) - f0
    f2 <- pgamma(x100, shape=s_shape, rate=s_rate) - f1 - f0
    f3 <- 1.0 - f2 - f1 - f0
    #print(paste(f0, f1, f2, f3, sep=", "))
    #return(c(list(f0), list(f1), list(f2), list(f3)))
    return(c(f0 = f0, f1 = f1, f2 = f2, f3 = f3))
}
dom_grapes_gammazero_raw_results_wclasses <- dom_grapes_gammazero_raw_results %>%
  mutate(output = map2(GammaZero.negGmean, GammaZero.negGshape, DFE_proportions_grapes_dom)) %>%
  unnest_wider(output) %>%
  rename(t0 = f0) %>%  #had to rename columns bc i rerun the function
  rename(t1 = f1) %>% 
  rename(t2 = f2) %>% 
  rename(t3 = f3) %>% group_by(selfing,DFE)

dom_grapes_summary <- dom_grapes_gammazero_raw_results_wclasses %>% 
group_by(selfing, DFE) %>%
summarize(across(where(is.numeric), list(avg = mean, sd = sd)))



dom_grapes_simple_summary <- 
    dom_grapes_summary[
        c('DFE','selfing','t0_avg','t1_avg', 't2_avg', 't3_avg',
            't0_sd','t1_sd', 't2_sd', 't3_sd')] %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE])) %>% 
    mutate(output = pmap(list(true_mean, true_shape, 100), DFE_proportions_dfe_alpha)) %>% #modified so I assume gamma = 2Nes
    unnest_wider(output) %>%
    rename(
        'z0' = f0,
        'z1'= f1 ,
        'z2'= f2 ,
        'z3'= f3 ,
        'f0' = t0_avg , 
        'f1' = t1_avg , 
        'f2' = t2_avg , 
        'f3' = t3_avg ,
        'f0_sd'= t0_sd , 
        'f1_sd'= t1_sd , 
        'f2_sd'= t2_sd , 
        'f3_sd' = t3_sd)
# convert the data frame to a tidy format
domgrapes_df_tidy <- gather(dom_grapes_simple_summary, 
    key = "generation", value = "value", c(f0, f1, f2, f3)) 

domgrapes_df_true_tidy <- gather(dom_grapes_simple_summary, 
    key = "generation", value = "value", c(z0:z3)) %>%
    mutate(generation = recode(generation,
        z0 = 'f0', 
        z1 = 'f1', 
        z2 = 'f2', 
        z3 = 'f3'
    ))  %>% 
    select(1,13:14) %>% distinct() %>%
    mutate(selfing = "True0")



dfealpha_dom_summary2 <- dfealpha_dom_raw_results %>%
  group_by(selfing,DFE,matchname) %>%
  select(c(selfing,matchname,DFE,gamma,b)) 

## KEY DIFFERENCE HERE
# multiplying gammas by 2 because h = 0.25 not .5
domgrapes_gammazero_summmary2 <- dom_grapes_gammazero_raw_results %>%
  group_by(DFE,selfing) %>%
  select(c(matchname,DFE,selfing,GammaZero.negGmean, GammaZero.negGshape)) %>%
  mutate(
    GammaZero.negGmean = GammaZero.negGmean*2
  )

prediction_accuracy_table_dom <- dfealpha_dom_summary2 %>%
  left_join(domgrapes_gammazero_summmary2) %>% group_by(DFE,selfing) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE]))


prediction_accuracy_table_dom_gamma_melt <- prediction_accuracy_table_dom %>% 
  summarize(across(c(true_mean, gamma, GammaZero.negGmean), list(avg = mean, sd = sd))) %>%
  rename("true gamma_avg" = true_mean_avg,
    "DFEalpha gamma_avg" = gamma_avg,
    "Grapes gamma_avg" = GammaZero.negGmean_avg,
    "true gamma_sd" = true_mean_sd,
    "DFEalpha gamma_sd" = gamma_sd,
    "Grapes gamma_sd" = GammaZero.negGmean_sd) %>% 
  pivot_longer(
    cols = c(3:8), 
    names_to = c(".value", "variable"), 
    names_sep = "_", 
    values_drop_na = TRUE
  ) %>%

  melt(variable.name = "gamma_method", id.vars = c("DFE", "selfing", "variable")) %>%
    pivot_wider(id_cols = c("DFE", "selfing","gamma_method"), 
                        names_from = "variable", 
                        values_from = "value")
prediction_accuracy_table_dom_gamma_melt <- prediction_accuracy_table_dom_gamma_melt %>%
mutate(
    selfing_class = case_when(
      selfing == "0" ~ "0% Selfing",
      selfing == "50" ~ "50% Selfing",
      selfing == "80" ~ "80% Selfing",
      selfing == "90" ~ "90% Selfing",
      selfing == "95" ~ "95% Selfing",
      selfing == "99" ~ "99% Selfing"
    )
  ) 

gammaorder_dom <- c("true gamma", "DFEalpha gamma", "Grapes gamma")

ggplot(prediction_accuracy_table_dom_gamma_melt, aes(x = factor(gamma_method, levels = gammaorder_dom), y = avg, fill = factor(gamma_method, levels = gammaorder_dom))) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class), scales = "free_y") +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  #stat_summary(fun.data = "mean_sdl", geom = "bar", position = "dodge") +
  geom_errorbar(aes(ymin = avg-sd, ymax = avg+sd), position = position_dodge(width = 0.9), width = 0.2) +
  labs(x = "Gamma estimation method", y = "Gamma", fill = "Gamma estimation method", title = "Gamma estimation accuracy h=0.25") + 
  scale_fill_manual(values=c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(axis.text.x = element_blank(), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15),legend.text = element_text(size=15))

prediction_accuracy_table_dom_beta_melt <- prediction_accuracy_table_dom %>% 
  summarize(across(c(true_shape, b, GammaZero.negGshape), list(avg = mean, sd = sd))) %>%
  rename("true shape_avg" = true_shape_avg,
    "DFEalpha beta_avg" = b_avg,
    "Grapes beta_avg" = GammaZero.negGshape_avg,
    "true shape_sd" = true_shape_sd,
    "DFEalpha beta_sd" = b_sd,
    "Grapes beta_sd" = GammaZero.negGshape_sd) %>%
  pivot_longer(
    cols = c(3:8), 
    names_to = c(".value", "variable"), 
    names_sep = "_", 
    values_drop_na = TRUE
  ) %>%
   melt(variable.name = "gamma_method", id.vars = c("DFE", "selfing", "variable")) %>%
    pivot_wider(id_cols = c("DFE", "selfing", "gamma_method"), 
                        names_from = "variable", 
                        values_from = "value") %>%
mutate(
    selfing_class = case_when(
      selfing == "0" ~ "0% Selfing",
      selfing == "50" ~ "50% Selfing",
      selfing == "80" ~ "80% Selfing",
      selfing == "90" ~ "90% Selfing",
      selfing == "95" ~ "95% Selfing",
      selfing == "99" ~ "99% Selfing"
    )
  ) 
ggplot(prediction_accuracy_table_dom_beta_melt, aes(x = gamma_method, y = avg, fill = gamma_method)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class), scales = "free_y") +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  #stat_summary(fun.data = "mean_sdl", geom = "bar", position = "dodge") +
  geom_errorbar(aes(ymin = avg-sd, ymax = avg+sd), position = position_dodge(width = 0.9), width = 0.2) +
  labs(title = "Beta estimation accuracy h=0.25", x = "Method", y = "Beta", fill = "") +
  scale_fill_manual(values=c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(axis.text.x = element_blank(), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15),legend.text = element_text(size=15))





#work in progress
#dom_dataframe_of_truth <- #must run calculate_pi.r
#dom_summary_table %>% 
#    mutate(row_names = row.names(pos_summary_table),
#    DFE = str_extract(row_names, "(DFE)\\d+"),
#    selfing = str_extract(row_names, "(?<=selfing)\\d+"),) %>% 
#    tibble() %>%
#    group_by(DFE) %>%
#    mutate(true_mean = unlist(true_gammas[DFE])) %>%
#    mutate(true_shape = unlist(true_betas[DFE]))
#
#dom_dataframe_of_truth <- dom_dataframe_of_truth %>% 
#mutate(output = pmap(list(B, true_mean, true_shape), DFE_proportions_truth)) %>%
#    unnest_wider(output) 
#
###new tidy version with multiple replicates. Can remove other two versions above if this is working out
#dataframe_of_truth2_dom <- dataframe_of_truth2 %>% filter(str_detect(fullpath, "dom,")) %>%
#mutate(output = pmap(list(B, true_mean, true_shape), DFE_proportions_truth)) %>%
#    unnest_wider(output) %>%
#    group_by(selfing, DFE) %>%
#    summarize(across(where(is.numeric), list(avg = mean, sd = sd))) %>%
#    rename(f0 = f0_avg) %>%  #had to rename columns bc i rerun the function
#    rename(f1 = f1_avg) %>% 
#    rename(f2 = f2_avg) %>% 
#    rename(f3 = f3_avg) 

#pos_truth_tidy <- gather(dataframe_of_truth2_pos, 
#    key = "generation", value = "value", c(f0, f1, f2, f3)) %>%
#    mutate(selfing = case_when(
#        selfing == '0'  ~ 'true0_recalc', 
#        selfing == '50' ~ 'true50', 
#        selfing == '80' ~ 'true80', 
#        selfing == '90' ~ 'true90',
#        selfing == '95' ~ 'true95',
#        selfing == '99' ~ 'true99',
#        TRUE ~ selfing
#    ))
#
# currently lacks adjusted truth


dfealpha_dom_for_plotting <- bind_rows(df_true_tidy,dfealpha_dom_tidy) %>%
    select(1:4,32,34,36,38) %>%
    melt() %>% 
    mutate(value = ifelse(is.na(value), 0, value)) %>% 
    filter(variable == "value" | paste0(generation, "_sd") == variable)
    
dfealpha_dom_for_plotting$variable <- ifelse(grepl("_sd", dfealpha_dom_for_plotting$variable), "sd", "value")  

voodoo_dom <- pivot_wider(dfealpha_dom_for_plotting, 
  id_cols = c("generation","DFE","selfing"), names_from = "variable", values_from = "value") %>%
mutate(selfing = case_when(
        selfing == 'True0'  ~ 'truth', 
        selfing == 'true0_recalc' ~ 'true0',
        selfing == '100' ~ '99',
        TRUE ~ selfing 
)) %>% 
mutate(
    selfing_class = case_when(
      selfing == "0" ~ "0% Selfing",
      selfing == "50" ~ "50% Selfing",
      selfing == "80" ~ "80% Selfing",
      selfing == "90" ~ "90% Selfing",
      selfing == "95" ~ "95% Selfing",
      selfing == "99" ~ "99% Selfing",
    )
  )


dfealpha_dom_rainbow_plot <- voodoo_dom %>% 
  filter(!grepl("F_adjusted|true", selfing)) 
  #filter(grepl("truth", selfing))

selfing_order <- c("truth", 0, 50, 80, 90, 95, 99)
ggplot(dfealpha_dom_rainbow_plot, aes(x = generation, y = value, fill = factor(selfing, 
    levels = selfing_order))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 1) +
  labs(title = "DFEalpha (h=0.25)", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))


dom_grapes_for_plotting <- bind_rows(domgrapes_df_tidy,domgrapes_df_true_tidy) %>%
    select(1:6,13:14) %>%
    melt() %>% 
    mutate(value = ifelse(is.na(value), 0, value)) %>% 
    filter(variable == "value" | paste0(generation, "_sd") == variable)

dom_grapes_for_plotting$variable <- ifelse(grepl("_sd", dom_grapes_for_plotting$variable), "sd", "value")  

dom_voodoo_grapes <- pivot_wider(dom_grapes_for_plotting, id_cols = c("generation","DFE","selfing"), 
    names_from = "variable", values_from = "value") 
dom_voodoo_grapes <- dom_voodoo_grapes %>% 
    mutate(selfing = case_when(
        selfing == "True0" ~ "truth",
        selfing == "true0_recalc" ~ "true0",
        selfing == 0~ "0_grapes",
        selfing == 50 ~ "50_grapes",
        selfing == 80 ~ "80_grapes",
        selfing == 90 ~ "90_grapes",
        selfing == 95 ~ "95_grapes",
        selfing == 99 ~ "99_grapes",
        TRUE ~ selfing
    ))
dom_voodoo_grapes <- filter(dom_voodoo_grapes, !(selfing %in% c('true90', 'true80', 'true95'))) #leaving filtering for last so that I can toggle what I want in the chart
#pos_voodoo_grapes <- filter(pos_voodoo_grapes, !(DFE %in% c('DFE3'))) # was used before DFE2 run was complete
selfing_order <- c("truth", "true0", "0_grapes", "true50", "50_grapes", "true99", "99_grapes")

# Create the grouped bar chart with custom selfing order
ggplot(dom_voodoo_grapes, aes(x = generation, y = value, fill = factor(selfing, 
    levels = selfing_order))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 1) +
  labs(title = "DFEalpha (h=0.25)", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))


dom_voodoo_grapes2_truth <- dom_voodoo_grapes %>% filter(grepl("truth", selfing))
dom_voodoo_grapes2_not_truth <- dom_voodoo_grapes %>% filter(!grepl("truth", selfing))

# Replicate "truth" data frame
df_truth_rep <- dom_voodoo_grapes2_truth %>% slice(rep(1:n(), each = 6))

# Assign selfing class to replicated "truth" data frame
df_truth_rep_self <- df_truth_rep %>%
  mutate(
    replicate = rep(1:6, length.out = n()),
    selfing_class = case_when(
      replicate == 1 ~ "0% Selfing",
      replicate == 2 ~ "50% Selfing",
      replicate == 3 ~ "80% Selfing",
      replicate == 4 ~ "90% Selfing",
      replicate == 5 ~ "95% Selfing",
      replicate == 6 ~ "99% Selfing",
      TRUE ~ "error"
    )
  ) %>%
  filter(!(replicate %in% c(3, 4, 5))) %>% #leaving this here in case we add more selfing levels
  select(-replicate)

dom_voodoo_grapes2 <- dom_voodoo_grapes2_not_truth %>%
mutate(
    selfing_class = case_when(
      selfing == "0_grapes" ~ "0% Selfing",
      selfing == "50_grapes" ~ "50% Selfing",
      #selfing == "80_grapes" ~ "80% Selfing",
      #selfing == "90_grapes" ~ "90% Selfing",
      #selfing == "95_grapes" ~ "95% Selfing",
      selfing == "99_grapes" ~ "99% Selfing",
      selfing == paste0("true", selfing_nums[1]) ~ "0% Selfing",
      selfing == paste0("true", selfing_nums[2]) ~ "50% Selfing",
      #selfing == paste0("true", selfing_nums[3]) ~ "80% Selfing",
      #selfing == paste0("true", selfing_nums[4]) ~ "90% Selfing",
      #selfing == paste0("true", selfing_nums[5]) ~ "95% Selfing",
      selfing == paste0("true", selfing_nums[6]) ~ "99% Selfing",
    )
) 

#plot grapes results
ggplot(rbind(dom_voodoo_grapes2, df_truth_rep_self), aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("truth", "true0", "0_grapes", "true50", "50_grapes", "true80", "80_grapes", "true90", "90_grapes", "true95", "95_grapes", "true99", "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "Grapes (h=0.25)", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values = c("#404040", rep(c("purple"),6))) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))


combo_dom_plot <- bind_rows(voodoo_dom, dom_voodoo_grapes2) %>% na.omit()

ggplot(rbind(combo_dom_plot, df_truth_rep_self), aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("truth", "true0", "0", "0_grapes", "true50", "50", "50_grapes", "true80", "80_grapes", "true90", "90_grapes", "true95", "95_grapes", "true99",  "99", "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "DFE-alpha and Grapes (h=0.25)", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values = c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_blank(), legend.text = element_blank())





















#### Low_recomb results
dfealpha_dir <- "/nas/longleaf/home/adaigle/work/lowrecombination_inputsandoutputs/"
dfealpha_output_dirs <- paste(paste(dfealpha_dir, dir(dfealpha_dir, pattern = "DFE_alpha_lowrecom_output"), sep = ""),
    "/", dir(
    paste(dfealpha_dir, dir(dfealpha_dir, pattern = "DFE_alpha_lowrecom_output"), sep = ""), pattern = "selected")
    , sep = "")

#table of results with columns signifying selfing%, DFE, and output
dfealpha_colnames <- as.data.frame(c("N1", "N2", "t2", "Nw", "b", "Es", "f0","L"))
dfealpha_raw_results_lowrec <- tibble(
    name = list.files(path = dfealpha_output_dirs, pattern = "est_dfe.out"),
    fullpath = paste(dfealpha_output_dirs, "/", name,sep=""),
    selfing = str_extract(fullpath, "(?<=DFE_alpha_output_)\\d+") %>% if_else(is.na(.), "0", .), #0's were empty
    matchname = str_extract(fullpath, "(DFE)\\d+(output)\\d"), #will be useful for comparing results later, or combining with other programs
    DFE = str_extract(matchname, "(DFE)\\d+"),
    data = lapply(fullpath, function(x) 
        as.data.frame(t(data.frame(matrix(unlist(strsplit(readLines(x), " ")), ncol = 2, byrow = TRUE)))) %>%
            rownames_to_column() %>%
            `colnames<-`(.[1,]) %>%
            .[-1,] %>%
            `rownames<-`(NULL))
)

df <- dfealpha_raw_results_lowrec$data %>% bind_rows()
df <- df %>% mutate(row_number = row_number())
dfealpha_raw_results_lowrec <- dfealpha_raw_results_lowrec  %>% 
    mutate(row_number = row_number())

# Join the original table with the new data frame using the row numbers as the key
dfealpha_raw_results_lowrec <- dfealpha_raw_results_lowrec %>% 
    left_join(df, by = "row_number") %>% 
    select(-data, -row_number) %>%
    as.data.frame %>% 
    mutate_at(vars(7:14), as.numeric) %>%
    mutate(gamma = -2*Nw*Es) %>%
    rename(f0_fromoutput = f0)

dfealpha_raw_results_lowrec_wtruth <- dfealpha_raw_results_lowrec %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE]))

#now i run the class generator function on all outputs
dfealpha_raw_results_lowrec_wclasses <- dfealpha_raw_results_lowrec_wtruth %>%
  mutate(output = pmap(list(gamma, b, Nw), DFE_proportions_dfe_alpha)) %>%
  unnest_wider(output) %>%
  rename(z0 = f0) %>%  #had to rename columns bc i rerun the function
  rename(z1 = f1) %>% 
  rename(z2 = f2) %>% 
  rename(z3 = f3) %>% 
  #rename( c(t0,t1,t2,t3) = c(f0,f1,f2,f3)) %>%
  mutate(output = pmap(list(gamma, b, Nw), DFE_proportions_dfe_alpha)) %>%
  unnest_wider(output) 
  #mutate(across(.fns = ~ if(all(!is.na(as.numeric(.x)))) as.numeric(.x) else .x))

dfealpha_lowrec_summary <- dfealpha_raw_results_lowrec_wclasses %>% 
    group_by(selfing, DFE) %>%
    summarize(across(where(is.numeric), list(avg = mean, sd = sd))) %>%
    rename(z0 = z0_avg) %>%  #had to rename columns bc i rerun the function
    rename(z1 = z1_avg) %>% 
    rename(z2 = z2_avg) %>% 
    rename(z3 = z3_avg) %>%
    rename(f0 = f0_avg) %>%  #had to rename columns bc i rerun the function
    rename(f1 = f1_avg) %>% 
    rename(f2 = f2_avg) %>% 
    rename(f3 = f3_avg) 

  # convert the data frame to a tidy format
dfealpha_lowrec_tidy <- gather(dfealpha_lowrec_summary, 
    key = "generation", value = "value", c(z0,z1,z2,z3)) %>%
    mutate(generation = recode(generation,
    z0 = 'f0', 
    z1 = 'f1', 
    z2 = 'f2', 
    z3 = 'f3')) %>%
    mutate(selfing = as.character(selfing))

### begin grapes lowrec
grapes_lowrec_dir <- "/nas/longleaf/home/adaigle/work/lowrecombination_inputsandoutputs/"
output_dirs <- paste(grapes_lowrec_dir, dir(grapes_lowrec_dir, pattern = "grapes_lowrecom_output"), "/all", sep = "")

#table of results with columns signifying selfing%, DFE, and output
grapes_lowrec_raw_results <- tibble(
    name = list.files(path = output_dirs, pattern = ".csv$"),
    matchname = sub(".txt.csv","",name), #will be useful for comparing results later, or combining with other programs
    DFE = str_extract(matchname, "(DFE)\\d+"),
    fullpath = paste(output_dirs, "/", name,sep=""),
    selfing = 0,
    data = lapply(fullpath,read.csv)
)

#now I want to clean this up because it has too much going on
# I will keep just the gammazero row in each DF
grapes_lowrec_gammazero_raw_results <- grapes_lowrec_raw_results %>% 
    mutate(data = map(data, ~ { #start of a mapping funciton, which is applied to each df in column
    df_filtered <- .x %>% filter(model == 'GammaZero') %>% #select only GammaZero row
    select(1:52) %>%# removing prediction columns for now
    select_if(~ all(!is.na(.))) # remove NA columns
    return(df_filtered)
}))

#now i summarize the outputs by averaging them and finding the standard deviation
#first I turn my table of dataframes into a dataframe, because my df's only had one row 
# Extract the column of data frames and bind them into a single data frame
df <- grapes_lowrec_gammazero_raw_results$data %>% bind_rows()

# Add a column to the data frame that contains the row numbers
df <- df %>% mutate(row_number = row_number())
grapes_lowrec_gammazero_raw_results <- grapes_lowrec_gammazero_raw_results  %>% 
    mutate(row_number = row_number())

# Join the original table with the new data frame using the row numbers as the key
grapes_lowrec_gammazero_raw_results <- grapes_lowrec_gammazero_raw_results %>% 
    left_join(df, by = "row_number") %>% 
    select(-data, -row_number) %>%
    as.data.frame

grapes_lowrec_raw_results_wtruth <- grapes_lowrec_gammazero_raw_results %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE]))

grapes_lowrec_raw_results_wclasses <- grapes_lowrec_gammazero_raw_results %>%
  mutate(output = map2(GammaZero.negGmean, GammaZero.negGshape, DFE_proportions_grapes)) %>%
  unnest_wider(output) %>%
  rename(t0 = f0) %>%  #had to rename columns bc i rerun the function
  rename(t1 = f1) %>% 
  rename(t2 = f2) %>% 
  rename(t3 = f3) 

grapes_lowrec_gammazero_summary <- grapes_lowrec_raw_results_wclasses %>% 
    group_by(selfing, DFE) %>%
    summarize(across(where(is.numeric), list(avg = mean, sd = sd)))

grapes_lowrec_gammazero_simple_summary <- 
    grapes_lowrec_gammazero_summary[
        c('DFE','selfing','t0_avg','t1_avg', 't2_avg', 't3_avg',
            't0_sd','t1_sd', 't2_sd', 't3_sd')] %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE])) %>% 
    mutate(output = pmap(list(true_mean, true_shape, 100), DFE_proportions_dfe_alpha)) %>%
    unnest_wider(output) 


dfealpha_lowrec_summary2 <- dfealpha_raw_results_lowrec %>%
  group_by(DFE,matchname) %>%
  select(c(matchname,DFE,gamma,b)) 

grapes_lowrec_gammazero_summmary2 <- grapes_lowrec_gammazero_raw_results %>%
  group_by(DFE) %>%
  select(c(matchname,DFE,GammaZero.negGmean, GammaZero.negGshape)) 

prediction_accuracy_table_lowrec <- dfealpha_lowrec_summary2 %>%
  left_join(grapes_lowrec_gammazero_summmary2) %>% group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE]))

prediction_accuracy_table_lowrec_gamma_melt <- prediction_accuracy_table_lowrec %>% 
  summarize(across(c(true_mean, gamma, GammaZero.negGmean), list(avg = mean, sd = sd))) %>%
  rename("true gamma_avg" = true_mean_avg,
    "DFEalpha gamma_avg" = gamma_avg,
    "Grapes gamma_avg" = GammaZero.negGmean_avg,
    "true gamma_sd" = true_mean_sd,
    "DFEalpha gamma_sd" = gamma_sd,
    "Grapes gamma_sd" = GammaZero.negGmean_sd) %>% 
  pivot_longer(
    cols = c(2:7), 
    names_to = c(".value", "variable"), 
    names_sep = "_", 
    values_drop_na = TRUE
  ) %>%

  melt(variable.name = "gamma_method", id.vars = c("DFE", "variable")) %>%
    pivot_wider(id_cols = c("DFE", "gamma_method"), 
                        names_from = "variable", 
                        values_from = "value")

gammaorder_lowrec <- c("true gamma", "DFEalpha gamma", "Grapes gamma")
ggplot(prediction_accuracy_table_lowrec_gamma_melt, aes(x = factor(gamma_method, levels = gammaorder_lowrec), y = avg, fill = factor(gamma_method, levels = gammaorder))) +
  facet_grid(rows = vars(DFE), scales = "free_y") +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  #stat_summary(fun.data = "mean_sdl", geom = "bar", position = "dodge") +
  geom_errorbar(aes(ymin = avg-sd, ymax = avg+sd), position = position_dodge(width = 0.9), width = 0.2) +
  labs(x = "Gamma estimation method", y = "Gamma", fill = "Gamma estimation method", title = "Gamma estimation accuracy (low recombination, no selfing)") + 
  scale_fill_manual(values=c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(axis.text.x = element_blank(), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15),legend.text = element_text(size=15))


prediction_accuracy_table_beta_melt_lowrec <- prediction_accuracy_table_lowrec %>% 
  summarize(across(c(true_shape, b, GammaZero.negGshape), list(avg = mean, sd = sd))) %>%
  rename("true shape_avg" = true_shape_avg,
    "DFEalpha beta_avg" = b_avg,
    "Grapes beta_avg" = GammaZero.negGshape_avg,
    "true shape_sd" = true_shape_sd,
    "DFEalpha beta_sd" = b_sd,
    "Grapes beta_sd" = GammaZero.negGshape_sd) %>%
  pivot_longer(
    cols = c(2:7), 
    names_to = c(".value", "variable"), 
    names_sep = "_", 
    values_drop_na = TRUE
  ) %>%
   melt(variable.name = "gamma_method", id.vars = c("DFE", "variable")) %>%
    pivot_wider(id_cols = c("DFE", "gamma_method"), 
                        names_from = "variable", 
                        values_from = "value")

ggplot(prediction_accuracy_table_beta_melt_lowrec, aes(x = gamma_method, y = avg, fill = gamma_method)) +
  facet_grid(rows = vars(DFE), scales = "free_y") +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  #stat_summary(fun.data = "mean_sdl", geom = "bar", position = "dodge") +
  geom_errorbar(aes(ymin = avg-sd, ymax = avg+sd), position = position_dodge(width = 0.9), width = 0.2) +
  labs(title = "Beta estimation accuracy (low recombination, no selfing)", x = "Method", y = "Beta", fill = "") +
  scale_fill_manual(values=c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(axis.text.x = element_blank(), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15),legend.text = element_text(size=15))




dfealpha_tidy <- gather(dfealpha_lowrec_summary, 
    key = "generation", value = "value", c(z0,z1,z2,z3)) %>%
    mutate(generation = recode(generation,
    z0 = 'f0', 
    z1 = 'f1', 
    z2 = 'f2', 
    z3 = 'f3')) %>%
    mutate(selfing = as.character(selfing))

dfealpha_for_plotting <- bind_rows(df_true_tidy,dfealpha_tidy) %>%
    select(1:4,32,34,36,38) %>%
    melt() %>% 
    mutate(value = ifelse(is.na(value), 0, value)) %>% 
    filter(variable == "value" | paste0(generation, "_sd") == variable)
    
dfealpha_for_plotting$variable <- ifelse(grepl("_sd", dfealpha_for_plotting$variable), "sd", "value")  

voodoo <- pivot_wider(dfealpha_for_plotting, 
  id_cols = c("generation","DFE","selfing"), names_from = "variable", values_from = "value") %>%
mutate(selfing = case_when(
        selfing == 'True0'  ~ 'truth', 
        selfing == 'true0_recalc' ~ 'true0',
        selfing == '100' ~ '99',
        TRUE ~ selfing 
))

dfealpha_rainbow_plot <- voodoo %>% 
  filter(!grepl("F_adjusted|true", selfing)) 
  #filter(grepl("truth", selfing))
selfing_order <- c("truth", 0)

ggplot(dfealpha_rainbow_plot, aes(x = generation, y = value, fill = factor(selfing, 
    levels = selfing_order))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 1) +
  labs(title = "DFEalpha low recombination test", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))

grapes_lowrec_gammazero_simple_summary <- grapes_lowrec_gammazero_simple_summary %>%
    rename(
        'z0' = f0,
        'z1'= f1 ,
        'z2'= f2 ,
        'z3'= f3 ,
       'f0' = t0_avg , 
       'f1' = t1_avg , 
       'f2' = t2_avg , 
        'f3' = t3_avg ,
        'f0_sd'= t0_sd , 
        'f1_sd'= t1_sd , 
        'f2_sd'= t2_sd , 
        'f3_sd' = t3_sd ) 

df_tidy_lowrec <- gather(grapes_lowrec_gammazero_simple_summary, 
    key = "generation", value = "value", c(f0, f1, f2, f3)) %>%
    mutate(
      selfing =as.character(selfing)
    )
grapes_lowrec_for_plotting <- bind_rows(df_true_tidy,df_tidy_lowrec) %>%
    select(1:8) %>%
    melt() %>% 
    mutate(value = ifelse(is.na(value), 0, value)) %>% 
    filter(variable == "value" | paste0(generation, "_sd") == variable)

grapes_lowrec_for_plotting$variable <- ifelse(grepl("_sd", grapes_lowrec_for_plotting$variable), "sd", "value")  


voodoo_grapes_lowrec <- pivot_wider(grapes_lowrec_for_plotting, id_cols = c("generation","DFE","selfing"), names_from = "variable", values_from = "value")  %>%
    mutate(selfing = case_when(
        selfing == "True0" ~ "truth",
        selfing == "true0_recalc" ~ "true0",
        selfing == 0~ "0_grapes",
        TRUE ~ selfing
    ))
#test plot with all the truth
#selfing_order <- c("truth", "F_adjusted_0","true0", 0, "F_adjusted_50","true50", 50, "F_adjusted_80","true80", 80, "F_adjusted_90","true90", 90, "F_adjusted_95","true95", 95, "F_adjusted_99","true99", 99)

#now we are removing all the rainbows from the plots, so this code removes the adjusted truths 
#also cleaning up names for better legend
grapes_lowrec_rainbow_plot <- voodoo_grapes_lowrec %>%
  filter(!grepl("F_adjusted|true", selfing)) 

selfing_order <- c("truth", "0_grapes")
# Create the grouped bar chart with custom selfing order
ggplot(grapes_lowrec_rainbow_plot, aes(x = generation, y = value, fill = factor(selfing, 
    levels = selfing_order))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 1) +
  labs(title = "Grapes low recombination test", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))


voodoo_grapes2_lowrec_truth <- voodoo_grapes_lowrec %>% filter(grepl("truth", selfing))
voodoo_grapes2_lowrec_not_truth <- voodoo_grapes_lowrec %>% filter(!grepl("truth", selfing))

# Replicate "truth" data frame
df_truth_rep <- voodoo_grapes2_lowrec_truth %>% slice(rep(1:n(), each = 6))

# Assign selfing class to replicated "truth" data frame
df_truth_rep_self <- df_truth_rep %>%
  mutate(
    replicate = rep(1:6, length.out = n()),
    selfing_class = case_when(
      replicate == 1 ~ "0% Selfing",
      #replicate == 2 ~ "50% Selfing",
      #replicate == 3 ~ "80% Selfing",
      #replicate == 4 ~ "90% Selfing",
      #replicate == 5 ~ "95% Selfing",
      #replicate == 6 ~ "99% Selfing",
      TRUE ~ "error"
    )
  ) %>%
  filter(!(replicate %in% c(3, 4, 5))) %>% #leaving this here in case we add more selfing levels
  select(-replicate)

voodoo_grapes2_ <- dom_voodoo_grapes2_not_truth %>%
mutate(
    selfing_class = case_when(
      selfing == "0_grapes" ~ "0% Selfing",
      selfing == "50_grapes" ~ "50% Selfing",
      #selfing == "80_grapes" ~ "80% Selfing",
      #selfing == "90_grapes" ~ "90% Selfing",
      #selfing == "95_grapes" ~ "95% Selfing",
      selfing == "99_grapes" ~ "99% Selfing",
      selfing == paste0("true", selfing_nums[1]) ~ "0% Selfing",
      selfing == paste0("true", selfing_nums[2]) ~ "50% Selfing",
      #selfing == paste0("true", selfing_nums[3]) ~ "80% Selfing",
      #selfing == paste0("true", selfing_nums[4]) ~ "90% Selfing",
      #selfing == paste0("true", selfing_nums[5]) ~ "95% Selfing",
      selfing == paste0("true", selfing_nums[6]) ~ "99% Selfing",
    )
) 


combo_lowrec_plot <- bind_rows(voodoo, grapes_lowrec_rainbow_plot) %>% na.omit() %>%
mutate(
    selfing_class = case_when(
      selfing == "0_grapes" ~ "0% Selfing",
      selfing == 0 ~ "0% Selfing",
      #selfing == "50_grapes" ~ "50% Selfing",
      #selfing == "80_grapes" ~ "80% Selfing",
      #selfing == "90_grapes" ~ "90% Selfing",
      #selfing == "95_grapes" ~ "95% Selfing",
      #selfing == "99_grapes" ~ "99% Selfing",
      selfing == "truth" ~ "0% Selfing"
      #selfing == paste0("true", selfing_nums[2]) ~ "50% Selfing",
      #selfing == paste0("true", selfing_nums[3]) ~ "80% Selfing",
      #selfing == paste0("true", selfing_nums[4]) ~ "90% Selfing",
      #selfing == paste0("true", selfing_nums[5]) ~ "95% Selfing",
      #selfing == paste0("true", selfing_nums[6]) ~ "99% Selfing",
    )
) 


selfing_order <- c("truth", 0, "0_grapes")
ggplot(combo_lowrec_plot, aes(x = generation, y = value, fill = factor(selfing, 
    levels = selfing_order))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "DFE-alpha and Grapes (low recombination)", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values = c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_blank(), legend.text = element_blank())
