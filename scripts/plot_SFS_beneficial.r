# Plots the SFS figures for sims with beneficial mutations.
# To reproduce, change paths to sim outputs as needed.

rm(list=ls())
library(tidyverse)

summarize_experiment_SFS <- function(selfing) {
path_to_files <- paste0("/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/positive/SFS/", selfing, "/")

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


f <- rep(c(1:10,11), times = c(rep(1, each =10),89))
DFE_list <- c("DFE1", "DFE2", "DFE3")
summarize_outputs<-list()
for (DFE in 1:length(DFE_list)) {
    neutral_last75 <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_neu_100_m1")), 1, function(x) {
      tapply(as.numeric(x[3:101]), f, sum)
    }))

    neutral_last75_prop <- t(apply(neutral_last75, 1, function(x) x/sum(x)))
    selected_last75 <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_sel_100_m2,m4")), 1, function(x) {
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


results <- lapply(c(0, 50, 80, 90, 95, 99), summarize_experiment_SFS)
plotting_df <- bind_rows(flatten(results))  

plotting_df_0909599 <- plotting_df %>% filter(selfing != 80 & selfing != 50& selfing != 90 & selfing != 95)

beneficial_1percent <- ggplot(plotting_df_0909599, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", fill = "SFS Type") +
  geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing)) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), 
        strip.text = element_text(size=15), plot.title= element_text(size=0), 
        legend.title = element_text(size=15), legend.text = element_text(size=15)) +
  #important addition to make x axis more readable
  scale_x_continuous(breaks = c(1, seq(5, 9, by = 5), 11),
                       labels = c(1, seq(5, 9, by = 5), "11+")) +
  expand_limits(y=c(0,0.7))


summarize_experiment_SFS_low <- function(selfing) {
path_to_files <- paste0("/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/low_positive/SFS/", selfing, "/")

#total neutral sites is 187500
neutral_sites <- 187500
#total selected sites are 562500
selected_sites <- 562500

#read in slim sfs and fixed counts to a list of dataframes
count_sfs_names <- list.files(path = path_to_files, pattern = "count.sfs")
count_sfs_df_list <- lapply(
    paste(path_to_files,count_sfs_names,sep=""), 
    function(x) read.table(x, header=TRUE))

#new table method. Not currently in use.
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


f <- rep(c(1:10,11), times = c(rep(1, each =10),89))
DFE_list <- c("DFE1", "DFE2", "DFE3")
summarize_outputs<-list()
for (DFE in 1:length(DFE_list)) {
    neutral_last75 <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_neu_100_m1")), 1, function(x) {
      tapply(as.numeric(x[3:101]), f, sum)
    }))

    neutral_last75_prop <- t(apply(neutral_last75, 1, function(x) x/sum(x)))
    selected_last75 <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_sel_100_m2,m4")), 1, function(x) {
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


results <- lapply(c(0, 50, 80, 90, 95, 99), summarize_experiment_SFS_low)
plotting_df <- bind_rows(flatten(results))  

plotting_df_0909599 <- plotting_df %>% filter(selfing != 80 & selfing != 50 & selfing != 90 & selfing != 95)
#FIGURE 2
beneficial_01percent <- ggplot(plotting_df_0909599, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", fill = "SFS Type") +
  geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing)) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), 
        strip.text = element_text(size=15), plot.title= element_text(size=0), 
        legend.title = element_text(size=15), legend.text = element_text(size=15)) +
  #important addition to make x axis more readable
  scale_x_continuous(breaks = c(1, seq(5, 9, by = 5), 11),
                       labels = c(1, seq(5, 9, by = 5), "11+")) +
  expand_limits(y=c(0,0.7))

sfigure10 <- ggarrange(beneficial_01percent, beneficial_1percent,
                    labels = c("A", "B"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "bottom")
ggsave("/nas/longleaf/home/adaigle/DFESelfing/figures_for_publication/sfigure10.svg", plot = sfigure10, width = 8.5, height = 10, dpi = 150)
