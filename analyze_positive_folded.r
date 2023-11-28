rm(list=ls())
library(RColorBrewer)
library(tidyverse)
source("/nas/longleaf/home/adaigle/DFESelfing/calculate_pi.r")
library(reshape2)
library(scales)
library(ggpubr)
grapes_dir <- "/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/low_positive/dfe_results/grapes/"
pos_output_dirs <- paste(grapes_dir, dir(grapes_dir, pattern = "grapes"), sep = "")
pos_output_dirs <- pos_output_dirs[grepl("grapes_folded_output", pos_output_dirs)]

selfing_nums <- c(0, 50, 80, 90, 95, 99)

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
}

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

#table of results with columns signifying selfing%, DFE, and output
pos_grapes_raw_results <- tibble(
    name = list.files(path = pos_output_dirs, pattern = ".csv$"),
    matchname = sub(".txt.csv","",name), #will be useful for comparing results later, or combining with other programs
    DFE = str_extract(matchname, "(DFE)\\d+"),
    fullpath = paste(pos_output_dirs, "/", name,sep=""),
    selfing = str_extract(fullpath, "(?<=grapes_folded_output_)\\d+"),
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

alphas <- read.table( file = "/nas/longleaf/home/adaigle/DFESelfing/alpha_list_v3.txt", header = TRUE)
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

#pos_dataframe_of_truth <- #must run calculate_pi.r
#pos_summary_table %>% 
#    mutate(row_names = row.names(pos_summary_table),
#    DFE = str_extract(row_names, "(DFE)\\d+"),
#    selfing = str_extract(row_names, "(?<=selfing)\\d+"),) %>% 
#    tibble() %>%
#    group_by(DFE) %>%
#    mutate(true_mean = unlist(true_gammas[DFE])) %>%
#    mutate(true_shape = unlist(true_betas[DFE]))
#
#pos_dataframe_of_truth <- pos_dataframe_of_truth %>% 
#mutate(output = pmap(list(B, true_mean, true_shape), DFE_proportions_truth)) %>%
#    unnest_wider(output) 

##new tidy version with multiple replicates. Can remove other two versions above if this is working out
dataframe_of_truth2 <-
  tidy_summary_table %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE]))
dataframe_of_truth2_pos <- dataframe_of_truth2 %>% #filter(str_detect(fullpath, "pos")) %>%
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
pos_voodoo_grapes <- filter(pos_voodoo_grapes, !(selfing %in% c('true90', 'true80', 'true95'))) %>% #leaving filtering for last so that I can toggle what I want in the chart
  filter(!grepl("F_adjusted|true", selfing)) 
#pos_voodoo_grapes <- filter(pos_voodoo_grapes, !(DFE %in% c('DFE3'))) # was used before DFE2 run was complete
selfing_order <- c("truth", "true0", "0_grapes", "true50", "50_grapes", "true99", "99_grapes")

# Create the grouped bar chart with custom selfing order
ggplot(pos_voodoo_grapes, aes(x = generation, y = value, fill = factor(selfing, levels = c("truth", "true0", "0_grapes", "true50", "50_grapes", "true99", "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 2) +
  labs(title = "Grapes positive selection deleterious DFEs", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1))
ggplot(pos_voodoo_grapes, aes(x = generation, y = value, fill = factor(selfing, 
    levels = selfing_order))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 1) +
  labs(title = "Grapes", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) + 
  scale_fill_manual(values = c("#404040", hue_pal()(6))) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))

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
  #filter(!(replicate %in% c(3, 4, 5))) %>% #leaving this here in case we add more selfing levels
  select(-replicate)
  


pos_voodoo_grapes2 <- pos_voodoo_grapes2_not_truth %>%
mutate(
    selfing_class = case_when(
      selfing == "0_grapes" ~ "0% Selfing",
      selfing == "50_grapes" ~ "50% Selfing",
      selfing == "80_grapes" ~ "80% Selfing",
      selfing == "90_grapes" ~ "90% Selfing",
      selfing == "95_grapes" ~ "95% Selfing",
      selfing == "99_grapes" ~ "99% Selfing",
      selfing == paste0("true", selfing_nums[1]) ~ "0% Selfing",
      selfing == paste0("true", selfing_nums[2]) ~ "50% Selfing",
      selfing == paste0("true", selfing_nums[3]) ~ "80% Selfing",
      selfing == paste0("true", selfing_nums[4]) ~ "90% Selfing",
      selfing == paste0("true", selfing_nums[5]) ~ "95% Selfing",
      selfing == paste0("true", selfing_nums[6]) ~ "99% Selfing",
    )
  ) 
posgrapesplot <- rbind(pos_voodoo_grapes2, df_truth_rep_self) %>%
  filter(!grepl("true", selfing)) %>% 
  filter(!grepl("error", selfing_class)) %>%
    mutate(selfing = recode(selfing,
     '0_grapes' = 'GRAPES',
     '50_grapes' = 'GRAPES',
     '80_grapes' = 'GRAPES',
     '90_grapes' = 'GRAPES',
     '95_grapes' = 'GRAPES',     
     '99_grapes' = 'GRAPES'))


#plot gammaexpo results
positive_grapes_plot <- ggplot(posgrapesplot, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("truth", "true0", "GRAPES", "true50", "50_grapes", "true80", "80_grapes", "true90", "90_grapes", "true95", "95_grapes", "true99", "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "Grapes (with beneficial mutations), folded SFS", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values = c("#404040", rep(c("purple"),6))) + 
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
    mutate(alpha_label = case_when(
        selfing == 'true_alpha_0'  ~ 'True alpha, 0% Selfing', 
        selfing == 'true_alpha_50' ~ 'True alpha, 50% Selfing', 
        selfing == 'true_alpha_99' ~ 'True alpha, 99% Selfing',
        selfing == '0' ~ 'Estimated alpha, 0% Selfing', 
        selfing == '50' ~ 'Estimated alpha, 50% Selfing', 
        selfing == '99' ~ 'Estimated alpha, 99% Selfing', 
        TRUE ~ selfing
    )) %>% 
    mutate(selfing_class = case_when(
        selfing == 'true_alpha_0'  ~ '0% Selfing', 
        selfing == 'true_alpha_50' ~ '50% Selfing', 
        selfing == 'true_alpha_99' ~ '99% Selfing',
        selfing == 'true_alpha_80'  ~ '80% Selfing', 
        selfing == 'true_alpha_90' ~ '90% Selfing', 
        selfing == 'true_alpha_95' ~ '95% Selfing',
        selfing == '80' ~ '80% Selfing', 
        selfing == '90' ~ '90% Selfing', 
        selfing == '95' ~ '95% Selfing', 
        selfing == '0' ~ '0% Selfing', 
        selfing == '50' ~ '50% Selfing', 
        selfing == '99' ~ '99% Selfing', 
        TRUE ~ selfing)) %>% 
    mutate(truth_pred = case_when(
        selfing == 'true_alpha_0'  ~ 'Truth', 
        selfing == 'true_alpha_50' ~ 'Truth', 
        selfing == 'true_alpha_99' ~ 'Truth',
        selfing == 'true_alpha_80'  ~ 'Truth', 
        selfing == 'true_alpha_90' ~ 'Truth', 
        selfing == 'true_alpha_95' ~ 'Truth',
        selfing == '0' ~ 'Prediction', 
        selfing == '50' ~ 'Prediction', 
        selfing == '99' ~ 'Prediction', 
        selfing == '80' ~ 'Prediction', 
        selfing == '90' ~ 'Prediction', 
        selfing == '95' ~ 'Prediction', 
        TRUE ~ selfing)) %>% 
    arrange(factor(truth_pred, levels=c("Truth", "Prediction")))

alpha_order <- c('True alpha, 0% Selfing', 'Estimated alpha, 0% Selfing', 
  'True alpha, 50% Selfing', 'Estimated alpha, 50% Selfing', 
  'True alpha, 99% Selfing', 'Estimated alpha, 99% Selfing')
positive_grapes_plot_alpha <- ggplot(alphaplot, aes(x=factor(truth_pred, levels=c("Truth", "Prediction")),y = alpha_avg, fill = factor(truth_pred, levels=c("Truth", "Prediction")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = alpha_avg - alpha_sd,
                    ymax = alpha_avg + alpha_sd),
                position = position_dodge(width = 0.9), width = 0.2) +
  labs(title = "Alpha prediction", x = "DFE", y = "alpha", fill = "") + 
  scale_fill_manual(values = c(rep(c("#404040", "purple"),6))) + 
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  expand_limits(y=c(0,1)) +
  theme(axis.text.x = element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15),legend.text = element_text(size=15))


figure5 <- ggarrange(positive_grapes_plot, positive_grapes_plot_alpha,
                    labels = c("A", "B"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "right")



dfealpha_dir <- "/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/low_positive/dfe_results/dfealpha/"
dfealpha_output_dirs <- paste(paste(dfealpha_dir, dir(dfealpha_dir, pattern = "DFE_alpha_autofold_output"), sep = ""),
    "/", dir(
    paste(dfealpha_dir, dir(dfealpha_dir, pattern = "DFE_alpha_autofold_output"), sep = ""), pattern = "selected")
    , sep = "") 
    
#dfealpha_output_dirs <- dfealpha_output_dirs[!grepl("DFE_alpha_autofold_output_50/DFE1", dfealpha_output_dirs)]
#dfealpha_output_dirs <- dfealpha_output_dirs[!grepl("99", dfealpha_output_dirs) | !grepl("output2", dfealpha_output_dirs)]

#table of results with columns signifying selfing%, DFE, and output
dfealpha_colnames <- as.data.frame(c("N1", "N2", "t2", "Nw", "b", "Es", "f0","L"))
dfealpha_raw_results <- tibble(
    name = list.files(path = dfealpha_output_dirs, pattern = "est_dfe.out"),
    fullpath = paste(dfealpha_output_dirs, "/est_dfe.out",sep=""),
    selfing = str_extract(fullpath, "(?<=DFE_alpha_autofold_output_)\\d+") %>% if_else(is.na(.), "0", .), #0's were empty
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
    mutate(gamma = -2*Nw*Es) %>% #MODIFIED TO ACCOUNT FOR H=0.1
    rename(f0_fromoutput = f0) #f0 is a dfealpha output, but I also need to name a class f0 later

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

dfealpha_for_plotting <- bind_rows(posgrapes_df_true_tidy,dfealpha_tidy) %>%
    select(1:4,32,34,36,38) %>%
    melt() %>% 
    mutate(value = ifelse(is.na(value), 0, value))

# filter sd's
dfealpha_for_plotting <- dfealpha_for_plotting %>% 
  filter(variable == "value" | paste0(generation, "_sd") == variable)

dfealpha_for_plotting$variable <- ifelse(grepl("_sd", dfealpha_for_plotting$variable), "sd", "value")  

voodoo <- pivot_wider(dfealpha_for_plotting, id_cols = c("generation","DFE","selfing"), names_from = "variable", values_from = "value")

dfealpha_for_plotting <- bind_rows(posgrapes_df_true_tidy,dfealpha_tidy,pos_truth_tidy) %>%
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

voodoo3 <- voodoo3 %>% 
  filter(!grepl("true", selfing)) %>%
  filter(!grepl("F_adjusted", selfing))

ggplot(voodoo3, aes(x = generation, y = value, fill = factor(selfing, 
    levels = selfingc_order))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "DFEalpha", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values = c("#404040", rep(c("#F8766D"),6)))+ 
  theme( axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), plot.title= element_text(size=25))

voodoo4 <- voodoo3 %>%
    mutate(selfing = recode(selfing,
     '0' = 'DFE-alpha',
     '50' = 'DFE-alpha',
     '80' = 'DFE-alpha',
     '90' = 'DFE-alpha',
     '95' = 'DFE-alpha',     
     '99' = 'DFE-alpha')) %>% 
     select(-replicate)

comboplot <- ggplot(rbind(voodoo4, posgrapesplot), aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("truth", "DFE-alpha", "GRAPES")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "DFEalpha", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values = c("#404040", rep(c("#F8766D", "purple"),6)))+ 
  theme( axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), plot.title= element_text(size=25))

ggarrange(comboplot, positive_grapes_plot_alpha,
                    labels = c("A", "B"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "right")