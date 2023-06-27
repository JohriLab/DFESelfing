# An R script to convert slim simulation SFS output to DFE alpha's input format
#TODO fix many hard coded directory paths. 
#Make sure input and output directories go where you want them to

#initialize. Eventually can make the path an argument or at least relative. 
rm(list=ls())
library(tidyverse)


#plot row as DF
#barplot(tapply(unlist(eqm_selfing99_DFE1_neu_100_m1[1, 3:102]), rep(1:10, each = 10), sum))
#
#barplot(tapply(unlist(eqm_selfing99_DFE1_sel_100_m2[1, 3:102]), rep(1:10, each = 10), sum))


#not using 0 and fixed class for now for visualization
#summarize_outputs_neutral <- tibble(
#    means = apply(eqm_selfing50_DFE2_neu_100_m1[3:101], 2, mean),
#    sd = apply(eqm_selfing50_DFE2_neu_100_m1[3:101], 2, sd),
#    entry_number = 1:length(means),
#    SFS = "neutral"
#    )
#summarize_outputs_selected <- tibble(
#    means = apply(eqm_selfing50_DFE2_sel_100_m2[3:101], 2, mean),
#    sd = apply(eqm_selfing50_DFE2_sel_100_m2[3:101], 2, sd),
#    entry_number = 1:length(means),
#    SFS = "selected"
#    )
#
#
#summarize_outputs <- rbind(summarize_outputs_selected, summarize_outputs_neutral)
#gg <- ggplot(summarize_outputs, aes(x = entry_number, y = means, fill = factor(SFS, levels = c("selected", "neutral")))) +
#  geom_bar(stat = "identity", position = "dodge", colour = "black") +
#  geom_errorbar(aes(ymin = means - sd, ymax = means + sd), position = position_dodge(width = 0.9)) +
#  labs(x = "Derived allele frequency", y = "Number of polymorphisms", title = "Site frequency spectra from deleterious-only simulations", fill = "SFS Type") +
#  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
#  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
#  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))
#ggsave("/nas/longleaf/home/adaigle/DFESelfing/sample_plot.png", gg, width=29.7, height=21, units="cm")


#plot attempt 2: split into 10 sections
#neutral_10 <- t(apply(eqm_selfing50_DFE2_neu_100_m1, 1, function(x) {
#  tapply(as.numeric(x[3:102]), rep(1:10, each = 10), sum)
#}))
#selected_10 <- t(apply(eqm_selfing50_DFE2_sel_100_m2, 1, function(x) {
#  tapply(as.numeric(x[3:102]), rep(1:10, each = 10), sum)
#}))
#
#summarize_outputs_neutral_10 <- tibble(
#    mean10 = apply(neutral_10, 2, mean),
#    sd10 = apply(neutral_10, 2, sd),
#    entry_number = 1:length(mean10),
#    SFS = "neutral"
#    )
#summarize_outputs_selected_10 <- tibble(
#    mean10 = apply(selected_10, 2, mean),
#    sd10 = apply(selected_10, 2, sd),
#    entry_number = 1:length(mean10),
#    SFS = "selected"
#    )
#
#summarize_outputs <- rbind(summarize_outputs_selected_10, summarize_outputs_neutral_10)
#ggplot(summarize_outputs, aes(x = entry_number, y = mean10, fill = factor(SFS, levels = c("selected", "neutral")))) +
#  geom_bar(stat = "identity", position = "dodge", colour = "black") +
#  geom_errorbar(aes(ymin = mean10 - sd10, ymax = mean10 + sd10), position = position_dodge(width = 0.9)) +
#  labs(x = "Derived allele frequency", y = "Number of polymorphisms", title = "Site frequency spectra from deleterious-only simulations", fill = "SFS Type") +
#  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
#  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
#  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))

#plot type 3: join last 75
#ignoring fixed class

summarize_experiment_SFS <- function(selfing) {
path_to_files <- paste0("/nas/longleaf/home/adaigle/SFS/", selfing, "/")
#path_to_DFESelfing <- paste0("/nas/longleaf/home/adaigle/rerun_dfealpha/DFE_alpha_input_", selfing, "/")
#path_to_dfe_alpha_output <- paste0("/nas/longleaf/home/adaigle/rerun_dfealpha/DFE_alpha_output_", selfing, "/")
#path_to_grapes_current_input <- paste"/nas/longleaf/home/adaigle/work/dominance_inputsandoutputs/grapes_input_50/"


#total neutral sites is 187500
neutral_sites <- 187500
#total selected sites are 562500
selected_sites <- 562500

#read in slim sfs and fixed counts to a list of dataframes
count_sfs_names <- list.files(path = path_to_files, pattern = "count.sfs")
count_sfs_df_list <- lapply(
    paste(path_to_files,count_sfs_names,sep=""), 
    function(x) read.table(x, header=TRUE))

#new table method. Not currently in use but will change to this eventually
sfs_table <- tibble(
    name = list.files(path = path_to_files, pattern = "count.sfs"),
    matchname = sub("_count.sfs","",name),
    fullpath = paste(path_to_files, "/", name,sep=""), 
    data = lapply(fullpath,read.table)
)

fixed_table <- tibble(
    name = list.files(path = path_to_files, pattern = "count.fixed"),
    matchname = sub("_count.fixed","",name),
    fullpath = paste(path_to_files, "/", name,sep=""), 
    data = lapply(fullpath,read.table)
)

#magic join witchcraft 
join_table <- inner_join(sfs_table, fixed_table, by="matchname")

#try mutate, make columns with names for differetn thing
#can use grepl or regex
join_table <- join_table %>% mutate(DFE = 
str_extract(matchname, "(?<=DFE)\\d+")) %>% 
mutate(across(where(is.character), as.factor)) %>% #make characters factors
mutate(data = map2(data.x, data.y, ~cbind(.x, .y[2]))) #join fixed num to df

#--------------------------------
#stripping _count.sfs from file names to match fixed_sfs header, to name list
names(count_sfs_df_list) <- lapply(count_sfs_names,
    function(x)
        sub("_count.sfs.*","",x)
)

fixed_sfs_names <- list.files(path = path_to_files, pattern = "count.fixed")
fixed_sfs_df_list <- lapply(
    paste(path_to_files,fixed_sfs_names,sep=""), 
    function(x) read.table(x, header=TRUE))
names(fixed_sfs_df_list) <- lapply(fixed_sfs_names,
    function(x)
        sub("_count.fixed.*","",x)
)

#next: for each df in count_sfs, find corresponding df in fixed, and bind X100 column 
#creates new object for each as well. 
#for loop could be vectorized eventually 
combined_df_names_list <- c()
for(x in names(count_sfs_df_list)) {
    combined_df_names_list <- append(combined_df_names_list, x)
    assign(x, cbind(
    count_sfs_df_list[grepl(x, names(count_sfs_df_list))][[1]],
    fixed_sfs_df_list[grepl(x, names(fixed_sfs_df_list))][[1]]$X100
))
}

#Get neu dfs, subtract sum of all rows from num neu sites
for(x in combined_df_names_list[grepl("neu", combined_df_names_list)]) {
    assign(x, 
    add_column(get(x), neutral_sites - rowSums(get(x)[2:101]), .before = 2))
}

#Get sel dfs, subtract sum of all rows from num sel sites
for(x in combined_df_names_list[grepl("sel_", combined_df_names_list)]) {
    assign(x, 
    add_column(get(x), selected_sites - rowSums(get(x)[2:101]), .before = 2))
}


f <- rep(c(1:25,26), times = c(rep(1, each =25),74))
DFE_list <- c("DFE1", "DFE2", "DFE3")
summarize_outputs<-list()
for (DFE in 1:length(DFE_list)) {
    neutral_last75 <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_neu_100_m1")), 1, function(x) {
      tapply(as.numeric(x[3:101]), f, sum)
    }))
    #neutral_last75_prop <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_neu_100_m1")), 1, function(x) {
    #  tapply(as.numeric(x[3:101]), f, function(x) as.numeric(x)/sum(as.numeric(x)))
    #}))
    neutral_last75_prop <- t(apply(neutral_last75, 1, function(x) x/sum(x)))
    selected_last75 <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_sel_100_m2")), 1, function(x) {
      tapply(as.numeric(x[3:101]), f, sum)
    }))
    #selected_last75_prop <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_sel_100_m2")), 1, function(x) {
    #  tapply(as.numeric(x[3:101]), f, function(x) as.numeric(x)/sum(as.numeric(x)))
    #}))
    selected_last75_prop <- t(apply(selected_last75, 1, function(x) x/sum(x)))

    summarize_outputs_neutral_last75 <- tibble(
        mean = apply(neutral_last75, 2, mean),
        sd = apply(neutral_last75, 2, sd),
        prop = apply(neutral_last75_prop, 2, mean),
        propsd = apply(neutral_last75_prop, 2, sd), 
        entry_number = 1:length(mean),
        SFS = "neutral"
        ) #%>% 
        #mutate(
        #    prop = mean/sum(mean),
        #    propsd = prop, 2, sd)
    #summarize_outputs_neutral_last75$proportion <- summarize_outputs_neutral_last75$mean / sum(summarize_outputs_neutral_last75$mean) 
    
    #summarize_outputs_neutral_last75$proportionsd <- apply(summarize_outputs_neutral_last75$proportion, 2, sd)

    summarize_outputs_selected_last75 <- tibble(
        mean = apply(selected_last75, 2, mean),
        sd = apply(selected_last75, 2, sd),
        prop = apply(selected_last75_prop, 2, mean),
        propsd = apply(selected_last75_prop, 2, sd), 
        entry_number = 1:length(mean),
        SFS = "selected"
        )# %>% 
        #mutate(
        #    prop = mean/sum(mean),
        #    propsd = prop, 2, sd)
    #summarize_outputs_selected_last75$mean <- summarize_outputs_selected_last75$mean / sum(summarize_outputs_selected_last75$mean) 
    #summarize_outputs_selected_last75$proportionsd <- apply(summarize_outputs_selected_last75$proportion, 2, sd)

    df <- as.data.frame(rbind(summarize_outputs_selected_last75, summarize_outputs_neutral_last75))
    df$DFE <- DFE_list[[DFE]]
    df$selfing <- selfing
    summarize_outputs$DFE_list[[DFE]] <- df
    }
return(summarize_outputs)
}


results <- lapply(c(0, 50, 80, 90, 95, 99), summarize_experiment_SFS)
plotting_df <- bind_rows(flatten(results))  
#total number
ggplot(plotting_df, aes(x = entry_number, y = mean, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  labs(x = "Derived allele frequency", y = "Number of polymorphisms", title = "Site frequency spectra from deleterious-only simulations", fill = "SFS Type") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.9)) +
  #expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing)) +
  #scale_fill_manual(values = c("#404040", rep(c("purple"),6))) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))
#proportion
ggplot(plotting_df, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  labs(x = "Derived allele frequency", y = "Proportion of polymorphisms", title = "Site frequency spectra from deleterious-only simulations", fill = "SFS Type") +
  geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  #expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing)) +
  #scale_fill_manual(values = c("#404040", rep(c("purple"),6))) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))

plotting_df_less <- plotting_df %>% filter(selfing != 80 & selfing != 90 & selfing != 95)

#total number
ggplot(plotting_df_less, aes(x = entry_number, y = mean, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  labs(x = "Derived allele frequency", y = "Number of polymorphisms", title = "Site frequency spectra from deleterious-only simulations", fill = "SFS Type") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.9)) +
  #expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing)) +
  #scale_fill_manual(values = c("#404040", rep(c("purple"),6))) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))

#simple SFS example
plotting_df_less<- plotting_df_less %>% filter(selfing == 0 & DFE == 'DFE3')
ggplot(plotting_df_less, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  labs(x = "Derived allele frequency", y = "Proportion of polymorphisms", title = "Site Frequency Spectrum (SFS)", fill = "SFS Type") +
  #geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  #expand_limits(y=c(0,1)) +
  #facet_grid(rows = vars(DFE), cols = vars(selfing)) +
  #scale_fill_manual(values = c("#404040", rep(c("purple"),6))) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()+
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15)) +
  scale_x_continuous(breaks = c(1, seq(5, 20, by = 5), 26),
                       labels = c(1, seq(5, 20, by = 5), "26+"), expand = c(0, 0))




plotting_df_809095 <- plotting_df %>% filter(selfing != 0 & selfing != 50 & selfing != 99)
#total number
ggplot(plotting_df_809095, aes(x = entry_number, y = mean, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  labs(x = "Derived allele frequency", y = "Number of polymorphisms", title = "Site frequency spectra from deleterious-only simulations", fill = "SFS Type") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.9)) +
  #expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing)) +
  #scale_fill_manual(values = c("#404040", rep(c("purple"),6))) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))
#proportion
ggplot(plotting_df_809095, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  labs(x = "Derived allele frequency", y = "Proportion of polymorphisms", title = "Site frequency spectra from deleterious-only simulations", fill = "SFS Type") +
  geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  #expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing)) +
  #scale_fill_manual(values = c("#404040", rep(c("purple"),6))) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))


plotting_df_0909599 <- plotting_df %>% filter(selfing != 80 & selfing != 50)
#total number
ggplot(plotting_df_0909599, aes(x = entry_number, y = mean, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  labs(x = "Derived allele frequency", y = "Number of polymorphisms", title = "Site frequency spectra from deleterious-only simulations", fill = "SFS Type") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.9)) +
  #expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing)) +
  #scale_fill_manual(values = c("#404040", rep(c("purple"),6))) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))
#proportion
ggplot(plotting_df_0909599, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  labs(x = "Derived allele frequency", y = "Proportion of polymorphisms", title = "Site frequency spectra from deleterious-only simulations", fill = "SFS Type") +
  geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  #expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing)) +
  #scale_fill_manual(values = c("#404040", rep(c("purple"),6))) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))

ggplot(plotting_df_0909599, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  labs(x = "Derived allele frequency", y = "Proportion of polymorphisms", title = "Site frequency spectra from deleterious-only simulations", fill = "SFS Type") +
  geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing)) +
  #scale_x_continuous(breaks = seq(1, max(plotting_df_0909599$entry_number), by = 9), 
  #                   labels = seq(1, max(plotting_df_0909599$entry_number), by = 9)) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25), 
        strip.text = element_text(size=15), plot.title= element_text(size=25), 
        legend.title = element_text(size=15), legend.text = element_text(size=15)) +
  #important addition to make x axis more readable
  scale_x_continuous(breaks = c(1, seq(5, 20, by = 5), 26),
                       labels = c(1, seq(5, 20, by = 5), "26+"))
#code to plot one SFS
#ggplot(summarize_outputs, aes(x = entry_number, y = mean10, fill = factor(SFS, levels = c("selected", "neutral")))) +
#  geom_bar(stat = "identity", position = "dodge", colour = "black") +
#  geom_errorbar(aes(ymin = mean10 - sd10, ymax = mean10 + sd10), position = position_dodge(width = 0.9)) +
#  labs(x = "Derived allele frequency", y = "Number of polymorphisms", title = "Site frequency spectra from deleterious-only simulations", fill = "SFS Type") +
#  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
#  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
#  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))

elegans_sfs <- read.table("/nas/longleaf/home/adaigle/DFESelfing/elegans_sfs.txt")
elegans_invariant_sfs <- read.table("/nas/longleaf/home/adaigle/DFESelfing/elegans_invariant_sfs.txt")
elegans_diverged_sfs <- read.table("/nas/longleaf/home/adaigle/DFESelfing/elegans_diverged_sfs.txt")
elegans_diverged_invariant_sfs <- read.table("/nas/longleaf/home/adaigle/DFESelfing/elegans_diverged_invariant_sfs.txt")

f <- rep(c(1:25,26), times = c(rep(1, each =25),26))
neutral_last75 <- t(apply(elegans_diverged_sfs[2,][2:52], 1, function(x) {
      tapply(as.numeric(x), f, sum)
    }))
    #neutral_last75_prop <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_neu_100_m1")), 1, function(x) {
    #  tapply(as.numeric(x[3:101]), f, function(x) as.numeric(x)/sum(as.numeric(x)))
    #}))
    neutral_last75_prop <- t(apply(neutral_last75, 1, function(x) x/sum(x)))
    selected_last75 <- t(apply(elegans_diverged_sfs[1,][2:52], 1, function(x) {
      tapply(as.numeric(x), f, sum)
    }))
    #selected_last75_prop <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_sel_100_m2")), 1, function(x) {
    #  tapply(as.numeric(x[3:101]), f, function(x) as.numeric(x)/sum(as.numeric(x)))
    #}))
    selected_last75_prop <- t(apply(selected_last75, 1, function(x) x/sum(x)))

    summarize_outputs_neutral_last75 <- tibble(
        mean = apply(neutral_last75, 2, mean),
        sd = apply(neutral_last75, 2, sd),
        prop = apply(neutral_last75_prop, 2, mean),
        propsd = apply(neutral_last75_prop, 2, sd), 
        entry_number = 1:length(mean),
        SFS = "neutral"
        ) #%>% 
        #mutate(
        #    prop = mean/sum(mean),
        #    propsd = prop, 2, sd)
    #summarize_outputs_neutral_last75$proportion <- summarize_outputs_neutral_last75$mean / sum(summarize_outputs_neutral_last75$mean) 
    
    #summarize_outputs_neutral_last75$proportionsd <- apply(summarize_outputs_neutral_last75$proportion, 2, sd)

    summarize_outputs_selected_last75 <- tibble(
        mean = apply(selected_last75, 2, mean),
        sd = apply(selected_last75, 2, sd),
        prop = apply(selected_last75_prop, 2, mean),
        propsd = apply(selected_last75_prop, 2, sd), 
        entry_number = 1:length(mean),
        SFS = "selected"
        )
df <- as.data.frame(rbind(summarize_outputs_selected_last75, summarize_outputs_neutral_last75))

#code to plot one SFS

ggplot(df, aes(x = entry_number, y = mean, fill = factor(SFS, levels = c("selected", "neutral")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.9)) +
  labs(x = "Derived allele frequency", y = "Number of polymorphisms", title = "C. elegans_diverged natural population SFS", fill = "SFS Type") +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))

ggplot(df, aes(x = entry_number, y = prop, fill = factor(SFS, levels = c("selected", "neutral")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.9)) +
  labs(x = "Derived allele frequency", y = "Number of polymorphisms", title = "C. elegans_diverged natural population SFS", fill = "SFS Type") +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))


f <- rep(c(1:25,26), times = c(rep(1, each =25),26))
neutral_last75 <- t(apply(elegans_diverged_invariant_sfs[2,][2:52], 1, function(x) {
      tapply(as.numeric(x), f, sum)
    }))
    #neutral_last75_prop <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_neu_100_m1")), 1, function(x) {
    #  tapply(as.numeric(x[3:101]), f, function(x) as.numeric(x)/sum(as.numeric(x)))
    #}))
    neutral_last75_prop <- t(apply(neutral_last75, 1, function(x) x/sum(x)))
    selected_last75 <- t(apply(elegans_diverged_invariant_sfs[1,][2:52], 1, function(x) {
      tapply(as.numeric(x), f, sum)
    }))
    #selected_last75_prop <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_sel_100_m2")), 1, function(x) {
    #  tapply(as.numeric(x[3:101]), f, function(x) as.numeric(x)/sum(as.numeric(x)))
    #}))
    selected_last75_prop <- t(apply(selected_last75, 1, function(x) x/sum(x)))

    summarize_outputs_neutral_last75 <- tibble(
        mean = apply(neutral_last75, 2, mean),
        sd = apply(neutral_last75, 2, sd),
        prop = apply(neutral_last75_prop, 2, mean),
        propsd = apply(neutral_last75_prop, 2, sd), 
        entry_number = 1:length(mean),
        SFS = "neutral"
        ) #%>% 
        #mutate(
        #    prop = mean/sum(mean),
        #    propsd = prop, 2, sd)
    #summarize_outputs_neutral_last75$proportion <- summarize_outputs_neutral_last75$mean / sum(summarize_outputs_neutral_last75$mean) 
    
    #summarize_outputs_neutral_last75$proportionsd <- apply(summarize_outputs_neutral_last75$proportion, 2, sd)

    summarize_outputs_selected_last75 <- tibble(
        mean = apply(selected_last75, 2, mean),
        sd = apply(selected_last75, 2, sd),
        prop = apply(selected_last75_prop, 2, mean),
        propsd = apply(selected_last75_prop, 2, sd), 
        entry_number = 1:length(mean),
        SFS = "selected"
        )
df <- as.data.frame(rbind(summarize_outputs_selected_last75, summarize_outputs_neutral_last75))
ggplot(df, aes(x = entry_number, y = mean, fill = factor(SFS, levels = c("selected", "neutral")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.9)) +
  labs(x = "Derived allele frequency", y = "Number of polymorphisms", title = "C. elegans_diverged natural population SFS", fill = "SFS Type") +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))

ggplot(df, aes(x = entry_number, y = prop, fill = factor(SFS, levels = c("selected", "neutral")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.9)) +
  labs(x = "Derived allele frequency", y = "Number of polymorphisms", title = "C. elegans_diverged natural population SFS", fill = "SFS Type") +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))