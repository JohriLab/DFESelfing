# A script to plot SFS's obtained from forward simulations
# In general the last 75 entries of the SFS are added together into one bar
# for easier visualization

# This script plots the SFS from the original experiments 
# and sims with varying dominance coefficients
rm(list=ls())
library(tidyverse)
library(ggpubr)
base_dir <- "/nas/longleaf/home/adaigle/DFESelfing/"
figures_dir <- paste0(base_dir, "figures_for_publication/")
sim_outputs_dir <- "/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/"

SFS_dir <- paste0(sim_outputs_dir, "original_simulations/SFS_plot/")
#SFS_dir <- paste0(sim_outputs_dir, "scaling/SFS_plot/")

summarize_experiment_SFS <- function(selfing, base_dir) {
path_to_files <- paste0(base_dir, selfing, "/")


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


#f <- rep(c(1:99))
f <- rep(c(1:10,11), times = c(rep(1, each =10),89))

DFE_list <- c("DFE1", "DFE2", "DFE3")
#DFE_list <- c("DFE3")
#DFE_list <- c("DFE2", "DFE3")

summarize_outputs<-list()
for (DFE in 1:length(DFE_list)) {
    neutral_last75 <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_neu_100_m1")), 1, function(x) {
      tapply(as.numeric(x[3:101]), f, sum)
    }))
    neutral_last75_prop <- t(apply(neutral_last75, 1, function(x) x/sum(x)))
    selected_last75 <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_sel_100_m2")), 1, function(x) {
      tapply(as.numeric(x[3:101]), f, sum)
    }))

    selected_last75_prop <- t(apply(selected_last75, 1, function(x) x/sum(x)))

    summarize_outputs_neutral_last75 <- tibble(
        mean = apply(neutral_last75, 2, mean),
        sd = apply(neutral_last75, 2, sd),
        prop = apply(neutral_last75_prop, 2, mean),
        propsd = apply(neutral_last75_prop, 2, sd), 
        entry_number = 1:length(mean),
        SFS = "neutral"
        ) 
    summarize_outputs_selected_last75 <- tibble(
        mean = apply(selected_last75, 2, mean),
        sd = apply(selected_last75, 2, sd),
        prop = apply(selected_last75_prop, 2, mean),
        propsd = apply(selected_last75_prop, 2, sd), 
        entry_number = 1:length(mean),
        SFS = "selected"
        )
    df <- as.data.frame(rbind(summarize_outputs_selected_last75, summarize_outputs_neutral_last75))
    df$DFE <- DFE_list[[DFE]]
    df$selfing <- selfing
    summarize_outputs$DFE_list[[DFE]] <- df
    }
return(summarize_outputs)
}

selfings <- c(0,50,80,90,95,99)
#selfings <- c(80,90,95,99)
#selfings <- c(99)

input_dir <- rep(SFS_dir,6)
results <- mapply(summarize_experiment_SFS, selfings, input_dir)
plotting_df <- bind_rows(flatten(results))  

plotting_df_0909599 <- plotting_df #%>% filter(selfing != 80 & selfing != 50)

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

#FIGURE 3
plotting_df_0909599 <- plotting_df_0909599 %>% 
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

#investigating second half of the sfs for signs of U shape
#plotting_df_0909599 <- plotting_df_0909599[plotting_df_0909599$entry_number >= 50 & plotting_df_0909599$entry_number <= 99, ]

ggplot(plotting_df_0909599, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  #geom_smooth(method = "auto", se = T) +
  #scale_y_log10() +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", fill = "Site type") +
  geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  #scale_x_continuous(breaks = seq(1, max(plotting_df_0909599$entry_number), by = 9), 
  #                   labels = seq(1, max(plotting_df_0909599$entry_number), by = 9)) +
  scale_fill_manual(values=c("#619CFF", "#F8766D")) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), 
        strip.text = element_text(size=15), plot.title= element_text(size=0), 
        legend.title = element_text(size=12), legend.text = element_text(size=12),
        legend.position = "bottom") +
  #important addition to make x axis more readable
  scale_x_continuous(breaks = c(1, seq(5, 9, by = 5), 11),
                       labels = c(1, seq(5, 9, by = 5), "11+"))

#ggsave(paste0(figures_dir, "figure3.svg"), plot = figure3, width = 8.5, height = 8, dpi = 300)
#ggsave(paste0(figures_dir, "figure3_poster.png"), plot = figure3, width = 8.5, height = 9, dpi = 300, units="in")


summarize_experiment_SFS_noclip <- function(selfing, base_dir) {
path_to_files <- paste0(base_dir, selfing, "/")


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


f <- rep(c(1:99))
#f <- rep(c(1:10,11), times = c(rep(1, each =10),89))

DFE_list <- c("DFE1", "DFE2", "DFE3")
#DFE_list <- c("DFE3")
#DFE_list <- c("DFE2", "DFE3")

summarize_outputs<-list()
for (DFE in 1:length(DFE_list)) {
    neutral_last75 <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_neu_100_m1")), 1, function(x) {
      tapply(as.numeric(x[3:101]), f, sum)
    }))
    neutral_last75_prop <- t(apply(neutral_last75, 1, function(x) x/sum(x)))
    selected_last75 <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_sel_100_m2")), 1, function(x) {
      tapply(as.numeric(x[3:101]), f, sum)
    }))

    selected_last75_prop <- t(apply(selected_last75, 1, function(x) x/sum(x)))

    summarize_outputs_neutral_last75 <- tibble(
        mean = apply(neutral_last75, 2, mean),
        sd = apply(neutral_last75, 2, sd),
        prop = apply(neutral_last75_prop, 2, mean),
        propsd = apply(neutral_last75_prop, 2, sd), 
        entry_number = 1:length(mean),
        SFS = "neutral"
        ) 
    summarize_outputs_selected_last75 <- tibble(
        mean = apply(selected_last75, 2, mean),
        sd = apply(selected_last75, 2, sd),
        prop = apply(selected_last75_prop, 2, mean),
        propsd = apply(selected_last75_prop, 2, sd), 
        entry_number = 1:length(mean),
        SFS = "selected"
        )
    df <- as.data.frame(rbind(summarize_outputs_selected_last75, summarize_outputs_neutral_last75))
    df$DFE <- DFE_list[[DFE]]
    df$selfing <- selfing
    summarize_outputs$DFE_list[[DFE]] <- df
    }
return(summarize_outputs)
}
selfings <- c(0,50,80,90,95,99)
#selfings <- c(80,90,95,99)
#selfings <- c(99)

input_dir <- rep(SFS_dir,6)
results2 <- mapply(summarize_experiment_SFS_noclip, selfings, input_dir)
plotting_df2 <- bind_rows(flatten(results2))  

plotting_df_0909599_2 <- plotting_df2 #%>% filter(selfing != 80 & selfing != 50)


#FIGURE 3
plotting_df_0909599_2 <- plotting_df_0909599_2 %>% 
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

ggplot(plotting_df_0909599_2, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  #geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  geom_smooth(method = "auto", se = T) +
  scale_y_log10() +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", fill = "Site type") +
  #geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  #scale_x_continuous(breaks = seq(1, max(plotting_df_0909599$entry_number), by = 9), 
  #                   labels = seq(1, max(plotting_df_0909599$entry_number), by = 9)) +
  scale_fill_manual(values=c("#619CFF", "#F8766D")) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), 
        strip.text = element_text(size=15), plot.title= element_text(size=0), 
        legend.title = element_text(size=12), legend.text = element_text(size=12),
        legend.position = "bottom") #+
  #important addition to make x axis more readable
  #scale_x_continuous(breaks = c(1, seq(5, 9, by = 5), 11),
  #                     labels = c(1, seq(5, 9, by = 5), "11+"))



plotting_df_0909599_DFE2 <- plotting_df_0909599 %>% filter(DFE=="DFE2") %>% 
    filter(selfing != 80 & selfing != 50)

plotting_df_0909599_2_DFE2 <- plotting_df_0909599_2 %>% filter(DFE=="DFE2") %>% 
    filter(selfing != 80 & selfing != 50)

clipped_DFEs <- ggplot(plotting_df_0909599_DFE2, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  #geom_smooth(method = "auto", se = T) +
  #scale_y_log10() +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", fill = "Site type") +
  geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  #scale_x_continuous(breaks = seq(1, max(plotting_df_0909599$entry_number), by = 9), 
  #                   labels = seq(1, max(plotting_df_0909599$entry_number), by = 9)) +
  scale_fill_manual(values=c("#619CFF", "#F8766D")) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), 
        strip.text = element_text(size=15), plot.title= element_text(size=0), 
        legend.title = element_text(size=12), legend.text = element_text(size=12),
        legend.position = "bottom") +
  #important addition to make x axis more readable
  scale_x_continuous(breaks = c(1, seq(5, 9, by = 5), 11),
                       labels = c(1, seq(5, 9, by = 5), "11+"))

log_dfes <- ggplot(plotting_df_0909599_2_DFE2, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_smooth(method = "auto", se = T, color = "black") +
  scale_y_log10() +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", fill = "Site type") +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values=c("#619CFF", "#F8766D")) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), 
        strip.text = element_text(size=15), plot.title= element_text(size=0), 
        legend.title = element_text(size=12), legend.text = element_text(size=12),
        legend.position = "bottom") +
  #important addition to make x axis more readable
  scale_x_continuous(breaks = c(1, 50, 99),
                       labels = c(1, 50, 99))


figure <- ggarrange(clipped_DFEs, log_dfes,
                    labels = c("A", "B"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "bottom", vjust=1) 
ggsave(paste0(figures_dir, "figure3.svg"), plot = figure, width = 8.5, height = 8, dpi = 300)

plotting_df_0909599_DFE1 <- plotting_df_0909599 %>% filter(DFE=="DFE1") %>% 
    filter(selfing != 80 & selfing != 50)

plotting_df_0909599_2_DFE1 <- plotting_df_0909599_2 %>% filter(DFE=="DFE1") %>% 
    filter(selfing != 80 & selfing != 50)

clipped_DFE1 <- ggplot(plotting_df_0909599_DFE1, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  #geom_smooth(method = "auto", se = T) +
  #scale_y_log10() +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", fill = "Site type") +
  geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  #scale_x_continuous(breaks = seq(1, max(plotting_df_0909599$entry_number), by = 9), 
  #                   labels = seq(1, max(plotting_df_0909599$entry_number), by = 9)) +
  scale_fill_manual(values=c("#619CFF", "#F8766D")) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), 
        strip.text = element_text(size=15), plot.title= element_text(size=0), 
        legend.title = element_text(size=12), legend.text = element_text(size=12),
        legend.position = "bottom") +
  #important addition to make x axis more readable
  scale_x_continuous(breaks = c(1, seq(5, 9, by = 5), 11),
                       labels = c(1, seq(5, 9, by = 5), "11+"))

log_dfe1 <- ggplot(plotting_df_0909599_2_DFE1, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_smooth(method = "auto", se = T, color = "black") +
  scale_y_log10() +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", fill = "Site type") +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values=c("#619CFF", "#F8766D")) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), 
        strip.text = element_text(size=15), plot.title= element_text(size=0), 
        legend.title = element_text(size=12), legend.text = element_text(size=12),
        legend.position = "bottom") +
  #important addition to make x axis more readable
  scale_x_continuous(breaks = c(1, 50, 99),
                       labels = c(1, 50, 99))


sfigure1 <- ggarrange(clipped_DFE1, log_dfe1,
                    labels = c("A", "B"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "bottom", vjust=1) 
ggsave(paste0(figures_dir, "sfigure5_new.svg"), plot = sfigure1, width = 8.5, height = 8, dpi = 300)

plotting_df_0909599_DFE3 <- plotting_df_0909599 %>% filter(DFE=="DFE3") %>% 
    filter(selfing != 80 & selfing != 50)

plotting_df_0909599_2_DFE3 <- plotting_df_0909599_2 %>% filter(DFE=="DFE3") %>% 
    filter(selfing != 80 & selfing != 50)

clipped_DFE3 <- ggplot(plotting_df_0909599_DFE3, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  #geom_smooth(method = "auto", se = T) +
  #scale_y_log10() +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", fill = "Site type") +
  geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  #scale_x_continuous(breaks = seq(1, max(plotting_df_0909599$entry_number), by = 9), 
  #                   labels = seq(1, max(plotting_df_0909599$entry_number), by = 9)) +
  scale_fill_manual(values=c("#619CFF", "#F8766D")) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), 
        strip.text = element_text(size=15), plot.title= element_text(size=0), 
        legend.title = element_text(size=12), legend.text = element_text(size=12),
        legend.position = "bottom") +
  #important addition to make x axis more readable
  scale_x_continuous(breaks = c(1, seq(5, 9, by = 5), 11),
                       labels = c(1, seq(5, 9, by = 5), "11+"))

log_DFE3 <- ggplot(plotting_df_0909599_2_DFE3, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_smooth(method = "auto", se = T, color = "black") +
  scale_y_log10() +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", fill = "Site type") +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values=c("#619CFF", "#F8766D")) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), 
        strip.text = element_text(size=15), plot.title= element_text(size=0), 
        legend.title = element_text(size=12), legend.text = element_text(size=12),
        legend.position = "bottom") +
  #important addition to make x axis more readable
  scale_x_continuous(breaks = c(1, 50, 99),
                       labels = c(1, 50, 99))


sfigure2 <- ggarrange(clipped_DFE3, log_DFE3,
                    labels = c("A", "B"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "bottom", vjust=1) 
ggsave(paste0(figures_dir, "sfigure6_new.svg"), plot = sfigure2, width = 8.5, height = 8, dpi = 300)
