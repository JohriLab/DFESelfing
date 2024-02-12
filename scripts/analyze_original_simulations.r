rm(list=ls())
library(RColorBrewer)
library(tidyverse)
library(reshape2)
library(scales)
library(ggpubr)

base_dir <- "/nas/longleaf/home/adaigle/DFESelfing/"
sim_outputs_dir <- "/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/"
results_dir <- paste0(sim_outputs_dir, "original_simulations/")
source(paste0(base_dir, "scripts/calculate_pi.r"))

figures_dir <- paste0(base_dir, "figures_for_publication/")
dfe_results_dir <- paste0(sim_outputs_dir, "original_simulations/dfe_results/")
dfealpha_dir <- paste0(dfe_results_dir, "dfealpha/")
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
  mutate(output = pmap(list(gamma, b, Nw), DFE_proportions_dfe_alpha)) %>%
  unnest_wider(output) 

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

#read in all selfing %, dfes, and experiments
grapes_dir <- paste0(dfe_results_dir, "grapes/")
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
    return(c(f0 = f0, f1 = f1, f2 = f2, f3 = f3))
}
grapes_raw_results_wclasses <- grapes_gammazero_raw_results %>%
  mutate(output = map2(GammaZero.negGmean, GammaZero.negGshape, DFE_proportions_grapes)) %>%
  unnest_wider(output) %>%
  rename(t0 = f0) %>%  #had to rename columns bc i rerun the function
  rename(t1 = f1) %>% 
  rename(t2 = f2) %>% 
  rename(t3 = f3)

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
    empircal_B_selfing_adjust = B/selfing_B, 
    F_adjusted_gamma = selfing_B * true_mean,
    deltaNe = selfing_Ne - empirical_Ne,
    newNE = 5000 - deltaNe,
    newNegamma = (newNE/5000) * true_mean) %>% group_by(DFE,selfing)

#write.csv(prediction_accuracy_table, file="/nas/longleaf/home/adaigle/DFESelfing/pylibseq/gammabeta.csv", quote=F)

## adding adjusted truths to plots
# will just append to end now so I don't break things
F_adjusted_classes <- prediction_accuracy_table %>%
  mutate(output = map2(F_adjusted_gamma, true_shape, DFE_proportions_grapes)) %>%
  unnest_wider(output) %>% 
  #select(DFE, selfing, f0, f1, f2, f3) %>%  # comment out this line for supp. table 3 info. 
  group_by(selfing, DFE) %>%
    summarize(across(where(is.numeric), list(avg = mean, sd = sd))) %>%
    rename(f0 = f0_avg) %>%  #had to rename columns bc i rerun the function
    rename(f1 = f1_avg) %>% 
    rename(f2 = f2_avg) %>% 
    rename(f3 = f3_avg) 

stable03 <- F_adjusted_classes %>% arrange(DFE,selfing) %>% select(DFE,selfing,empirical_Ne_avg,selfing_Ne_avg,selfing_B_avg,B_avg) %>%
  mutate(selfing_adjusted_B = B_avg/selfing_B_avg)


#write.csv(stable03, file="/nas/longleaf/home/adaigle/DFESelfing/stable03.csv", quote=F)

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
    return(c(f0 = f0, f1 = f1, f2 = f2, f3 = f3))
}

#dataframe_of_truth2_adjusted_B <- dataframe_of_truth2 %>% arrange(DFE, matchname)
dataframe_of_truth2_adjusted_B <- dataframe_of_truth2 %>% left_join(prediction_accuracy_table)

B_adjust_csv <- dataframe_of_truth2_adjusted_B %>% group_by(DFE,selfing) %>% summarize(across(where(is.numeric), list(avg = mean, sd = sd)))
write.csv(B_adjust_csv, "/nas/longleaf/home/adaigle/DFESelfing/B_adjust.csv")

dataframe_of_truth2_nopos <- dataframe_of_truth2_adjusted_B %>% filter(!str_detect(fullpath, "pos")) %>%
mutate(output = pmap(list(empircal_B_selfing_adjust, true_mean, true_shape), DFE_proportions_truth)) %>%
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

grapes_for_plotting <- bind_rows(df_true_tidy,df_tidy) %>%
    select(1:8) %>%
    melt() %>% 
    mutate(value = ifelse(is.na(value), 0, value)) %>% 
    filter(variable == "value" | paste0(generation, "_sd") == variable)

grapes_for_plotting$variable <- ifelse(grepl("_sd", grapes_for_plotting$variable), "sd", "value")  


voodoo_grapes <- pivot_wider(grapes_for_plotting, id_cols = c("generation","DFE","selfing"), names_from = "variable", values_from = "value")  
  
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
selfing_order <- c("True0", 0, 50, 80, 90, 95, 100)
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
        selfing == 'True0'  ~ 'Simulated DFE', 
        selfing == 'true0_recalc' ~ 'true0',
        selfing == '100' ~ '99',
        TRUE ~ selfing 
))

#removing adjusted truths for rainbow plots for now
dfealpha_rainbow_plot <- voodoo %>% 
  filter(!grepl("F_adjusted|true", selfing)) 
  #filter(grepl("truth", selfing))
selfing_order <- c("Simulated DFE", 0, 50, 80, 90, 95, 99)

dfealpha_rainbow <- ggplot(dfealpha_rainbow_plot, aes(x = generation, y = value, fill = factor(selfing, 
    levels = selfing_order))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 1) +
  labs(title = "", x = "", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) + 
  scale_fill_manual(values = c("#404040", hue_pal()(6))) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), 
    axis.title.x=element_text(size=15),axis.title.y=element_text(size=15), strip.text = element_text(size=15),
    plot.title= element_text(size=0), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(expression(italic(f[0])), expression(italic(f[1])), expression(italic(f[2])), expression(italic(f[3]))))

# Split data frame into "truth" and "non-truth" data frames
voodoo2_truth <- voodoo %>% filter(grepl("Simulated DFE", selfing))
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
  theme(legend.position="bottom")

### grapes
  grapes_for_plotting <- bind_rows(df_true_tidy,df_tidy,F_adjusted_classes_tidy, truth_tidy) %>%
    select(1:8) %>%
    melt() %>% 
    mutate(value = ifelse(is.na(value), 0, value)) %>% 
    filter(variable == "value" | paste0(generation, "_sd") == variable)

grapes_for_plotting$variable <- ifelse(grepl("_sd", grapes_for_plotting$variable), "sd", "value")  


voodoo_grapes <- pivot_wider(grapes_for_plotting, id_cols = c("generation","DFE","selfing"), names_from = "variable", values_from = "value")  %>%
    mutate(selfing = case_when(
        selfing == "True0" ~ "Simulated DFE",
        selfing == "true0_recalc" ~ "true0",
        TRUE ~ selfing
    ))

grapes_rainbow_plot <- voodoo_grapes %>%
  filter(!grepl("F_adjusted|true", selfing))
selfing_order <- c("Simulated DFE", 0, 50, 80, 90, 95, 99)

# Create the grouped bar chart with custom selfing order
grapes_rainbow <- ggplot(grapes_rainbow_plot, aes(x = generation, y = value, fill = factor(selfing, 
    levels = selfing_order))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 1) +
  labs(title = "", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) + 
  scale_fill_manual(values = c("#404040", hue_pal()(6))) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), 
    axis.title.x=element_text(size=15),axis.title.y=element_text(size=15), strip.text = element_text(size=15),
    plot.title= element_text(size=0), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(expression(italic(f[0])), expression(italic(f[1])), expression(italic(f[2])), expression(italic(f[3]))))

figure1 <- ggarrange(dfealpha_rainbow, grapes_rainbow,
                    labels = c("A", "B"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "bottom")
ggsave(paste0(figures_dir, "figure1.pdf"), plot = figure1, width = 8.5, height = 8.5, dpi = 300)

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

B_value_table <- B_value_table %>%
  mutate(
    selfing_class = case_when(
      selfing == "0" ~ "0% Selfing",
      selfing == "50" ~ "50% Selfing",
      selfing == "80" ~ "80% Selfing",
      selfing == "90" ~ "90% Selfing",
      selfing == "95" ~ "95% Selfing",
      selfing == "99" ~ "99% Selfing",
    ),
    B_avg = B_avg/(1/(1+(.01 *as.numeric(selfing)) / (2-.01*as.numeric(selfing))))
  )
B_value_table$pis_avg <- NULL
B_value_table$empirical_Ne_avg <- NULL
B_value_table$generation <- "f1"
B_value_table$value <- 0.8
B_value_table$selfing <- "GRAPES"
combo_plot <- bind_rows(voodoo3,voodoo3_grapes) %>%
  filter(!grepl("true", selfing)) %>%
  filter(!grepl("F_adjusted", selfing)) %>%
  filter(selfing_class != "80% Selfing") %>%
  filter(selfing_class != "90% Selfing") %>%
  filter(selfing_class != "95% Selfing") 

combo_plot_with_Badjusted <- bind_rows(voodoo3,voodoo3_grapes) %>%
  filter(!grepl("truth", selfing)) %>%
  filter(!grepl("F_adjusted", selfing)) %>%
  #filter(selfing_class != "80% Selfing") %>%
  #filter(selfing_class != "90% Selfing") %>%
  #filter(selfing_class != "95% Selfing") %>% 
  filter(selfing_class != "truth") %>% 
    mutate(selfing = recode(selfing,
     '0_grapes' = 'GRAPES',
     '50_grapes' = 'GRAPES',
     '80_grapes' = 'GRAPES',
     '90_grapes' = 'GRAPES',
     '95_grapes' = 'GRAPES',     
     '99_grapes' = 'GRAPES',
     '0' = 'DFE-alpha',
     '50' = 'DFE-alpha',
     '80' = 'DFE-alpha',
     '90' = 'DFE-alpha',
     '95' = 'DFE-alpha',     
     '99' = 'DFE-alpha',
     'true0' = "Adjusted DFE",
     'true50' = "Adjusted DFE",
     'true80' = "Adjusted DFE",
     'true90' = "Adjusted DFE",
     'true95' = "Adjusted DFE",     
     'true99' = "Adjusted DFE")) 
#B_value_table <- B_value_table %>%
#  mutate(selfing_class = factor(selfing_class, levels = c("0.5x", "0.1x", "0.05x", "0.01x", "0.005x", "0.001x")))


#write.csv(combo_plot, file=paste0(base_dir, "05099selfingbasic.csv"))

Badjusted <- ggplot(combo_plot_with_Badjusted, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", "Adjusted DFE", "DFE-alpha", "GRAPES")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "DFEalpha and Grapes ", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values = c("#404040", rep(c("grey", "#F8766D", "purple"),6))) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), 
    axis.title.x=element_text(size=0),axis.title.y=element_text(size=15), strip.text = element_text(size=13),
    plot.title= element_text(size=0), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(expression(italic(f[0])), expression(italic(f[1])), expression(italic(f[2])), expression(italic(f[3])))) +
  geom_text(data = B_value_table, 
    aes(label = paste("B = ", format(round(B_avg, 2), nsmall = 2, digits = 2))), 
    vjust = -0.5, hjust=0.25, size = 4, position = position_dodge(width = 0.9), fontface="italic") 



ggsave(paste0(figures_dir, "figure1.svg"), plot = Badjusted, width = 8.5, height = 7.5, dpi = 300)

#combo_plot <- combo_plot %>% mutate(dominance=0.5)
combo_plot_int <- combo_plot %>% mutate(intergenic=3000)
write.csv(combo_plot_int, file=paste0(base_dir, "scripts/intergenic_plot/int3000.csv"))

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
  rename("simulated gamma_avg" = true_mean_avg,
    "BGS-adjusted gamma_avg" = adjusted_gamma_avg,
    "F-adjusted gamma_avg" = F_adjusted_gamma_avg,
    "DFE-alpha_avg" = gamma_avg,
    "GRAPES_avg" = GammaZero.negGmean_avg,
    "simulated gamma_sd" = true_mean_sd,
    "BGS-adjusted gamma_sd" = adjusted_gamma_sd,
    "F-adjusted gamma_sd" = F_adjusted_gamma_sd,
    "DFE-alpha_sd" = gamma_sd,
    "GRAPES_sd" = GammaZero.negGmean_sd) %>% 
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
gammaorder <- c("simulated value", "F-adjusted gamma", "BGS-adjusted gamma", "DFE-alpha", "GRAPES")

prediction_accuracy_table_gamma_melt <- prediction_accuracy_table_gamma_melt %>%
mutate(
    selfing_class = case_when(
      selfing == "0" ~ "0% Selfing",
      selfing == "50" ~ "50% Selfing",
      selfing == "80" ~ "80% Selfing",
      selfing == "90" ~ "90% Selfing",
      selfing == "95" ~ "95% Selfing",
      selfing == "99" ~ "99% Selfing"
    ),

  gamma_method = case_when(gamma_method == "simulated gamma" ~ "simulated value", TRUE ~ gamma_method)) %>%
  arrange(factor(gamma_method, levels = gammaorder)) %>%
  filter(!grepl("F-adjusted", gamma_method)) %>%
  filter(!grepl("BGS", gamma_method))

#colorblind palete?
gamma_accuracy <- ggplot(prediction_accuracy_table_gamma_melt, aes(x = factor(gamma_method, levels = gammaorder), y = avg, fill = factor(gamma_method, levels = gammaorder))) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class), scales = "free_y") +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = avg-sd, ymax = avg+sd), position = position_dodge(width = 0.9), width = 0.2) +
  labs(x = "", y = expression(2 * N[e] * bar(s[d])), fill = "", title = "") + 
  scale_fill_manual(values=c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(axis.text.x=element_text(size=0), axis.text.y=element_text(size=15), 
    axis.title.x=element_text(size=0),axis.title.y=element_text(size=15), strip.text = element_text(size=15),
    plot.title= element_text(size=0), legend.position = "bottom", legend.text = element_text(size=12))


#beta stuff
prediction_accuracy_table_beta_melt <- prediction_accuracy_table %>% 
  summarize(across(c(true_shape, b, GammaZero.negGshape), list(avg = mean, sd = sd))) %>%
  rename("true shape_avg" = true_shape_avg,
    "DFE-alpha_avg" = b_avg,
    "GRAPES_avg" = GammaZero.negGshape_avg,
    "true shape_sd" = true_shape_sd,
    "DFE-alpha_sd" = b_sd,
    "GRAPES_sd" = GammaZero.negGshape_sd) %>%
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

beta_accuracy <- ggplot(prediction_accuracy_table_beta_melt, aes(x = gamma_method, y = avg, fill = gamma_method)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class), scales = "free_y") +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = avg-sd, ymax = avg+sd), position = position_dodge(width = 0.9), width = 0.2) +
  labs(title = "", x = "", y = expression(beta), fill = "") +
  scale_fill_manual(values=c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(axis.text.x=element_text(size=0), axis.text.y=element_text(size=15), 
    axis.title.x=element_text(size=0),axis.title.y=element_text(size=15), strip.text = element_text(size=15),
    plot.title= element_text(size=0), legend.position = "bottom", legend.text = element_text(size=12))

sfigure01 <- ggarrange(gamma_accuracy, beta_accuracy,
                    labels = c("A", "B"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "bottom")
ggsave(paste0(figures_dir, "sfigure01.svg"), plot = sfigure01, width = 8.5, height = 9.5, dpi = 150)
