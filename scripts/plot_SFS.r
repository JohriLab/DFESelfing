# A script to plot SFS's obtained from forward simulations
# In general the last 75 entries of the SFS are added together into one bar
# for easier visualization

# This script plots the SFS from the original experiments 
# and sims with varying dominance coefficients
rm(list=ls())
library(tidyverse)

base_dir <- "/nas/longleaf/home/adaigle/DFESelfing/"
figures_dir <- paste0(base_dir, "figures_for_publication/")
sim_outputs_dir <- "/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/"

SFS_dir <- paste0(sim_outputs_dir, "original_simulations/SFS_plot/")
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


f <- rep(c(1:10,11), times = c(rep(1, each =10),89))
DFE_list <- c("DFE1", "DFE2", "DFE3")
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

summarize_dom_experiment_SFS <- function(selfing, base_dir) {
path_to_files <- paste0(base_dir, selfing, "/")

#total num neutral sites is 187500
neutral_sites <- 187500
#total num selected sites is 562500
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
    neutral_last75 <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_hdel_0_25_neu_100_m1")), 1, function(x) {
      tapply(as.numeric(x[3:101]), f, sum)
    }))

    neutral_last75_prop <- t(apply(neutral_last75, 1, function(x) x/sum(x)))
    selected_last75 <- t(apply(get(paste0("eqm_selfing", selfing, "_", DFE_list[DFE], "_hdel_0_25_sel_100_m2")), 1, function(x) {
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
input_dir <- rep(SFS_dir,6)
results <- mapply(summarize_experiment_SFS, selfings, input_dir)
plotting_df <- bind_rows(flatten(results))  

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

figure3 <- ggplot(plotting_df_0909599, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
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

ggsave(paste0(figures_dir, "figure3.svg"), plot = figure3, width = 8.5, height = 8, dpi = 300)
ggsave(paste0(figures_dir, "figure3_poster.png"), plot = figure3, width = 8.5, height = 9, dpi = 300, units="in")

#dominance experiment: h = 0.25
domh025_dir <- paste0(sim_outputs_dir, "dom/hdel_0_25/SFS/")
selfings_h025 <- c(0,50,99)
input_dir <- rep(domh025_dir,3)
results_h025 <- mapply(summarize_dom_experiment_SFS, selfings_h025, input_dir)
plotting_df_h025 <- bind_rows(flatten(results_h025))

plotting_df_h025 <- plotting_df_h025 %>% 
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

hdel025_linked_figure <- ggplot(plotting_df_h025, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", fill = "Site type") +
  geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25), 
        strip.text = element_text(size=15), plot.title= element_text(size=25), 
        legend.title = element_text(size=15), legend.text = element_text(size=15)) +
  expand_limits(y=c(0,0.7)) +
  #important addition to make x axis more readable
  scale_x_continuous(breaks = c(1, seq(5, 9, by = 5), 11),
                       labels = c(1, seq(5, 9, by = 5), "11+"))
#ggsave("/nas/longleaf/home/adaigle/DFESelfing/figures_for_publication/hdel025_linked_figure.svg", plot = hdel025_linked_figure, width = 8.5, height = 7, dpi = 600)


plotting_df_05099 <- plotting_df %>% filter(selfing != 80 & selfing != 90 & selfing != 95)
plotting_df_05099 <- plotting_df_05099 %>% 
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

hdel05_linked_figure <- ggplot(plotting_df_05099, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", fill = "Site type") +
  geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25), 
        strip.text = element_text(size=15), plot.title= element_text(size=25), 
        legend.title = element_text(size=15), legend.text = element_text(size=15)) +
  #important addition to make x axis more readable
  expand_limits(y=c(0,0.7)) +
  scale_x_continuous(breaks = c(1, seq(5, 9, by = 5), 11),
                       labels = c(1, seq(5, 9, by = 5), "11+"))

#unlinked dominance experiment: h = 0.5
domh05_unlinked_dir <- paste0(sim_outputs_dir, "nolinkage/final_50k_outputs/SFS_for_plotting/")
selfings_h05 <- c(0,50,99)
input_unlinkedh05_dir <- rep(domh05_unlinked_dir,3)
results_h05_unlinked <- mapply(summarize_experiment_SFS, selfings_h05, input_unlinkedh05_dir)
plotting_df_h05_unlinked <- bind_rows(flatten(results_h05_unlinked))

plotting_df_h05_unlinked <- plotting_df_h05_unlinked %>% 
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

plotting_df_h05_unlinked_new <- plotting_df_h05_unlinked %>% 
   filter(selfing != 50 | SFS != "neutral") %>% #filter(selfing != 99 | SFS != "neutral") %>%
    mutate(selfing = case_when(
      selfing == 0 & SFS == "neutral" ~ "Neutral, 0% selfing",
      selfing == 99 & SFS == "neutral" ~ "Neutral, 99% selfing",
      TRUE ~ as.character(selfing)))


plotting_df_h05_unlinked_neutral <- plotting_df_h05_unlinked %>% 
  filter(selfing == 99 | SFS == "neutral") %>% filter(selfing == 50 | SFS == "neutral") #%>%
    #mutate(selfing = case_when(
    #  selfing == 0 & SFS == "neutral" ~ "neutral",
    #  TRUE ~ as.character(selfing)))

# Reorder the 'selfing' column in your data frame
plotting_df_h05_unlinked_new$selfing <- factor(plotting_df_h05_unlinked_new$selfing, levels = c("Neutral, 0% selfing", "Neutral, 99% selfing", "0", "50", "99"))

# Plot with reordered factor levels
hdel05_unlinked_figure <- ggplot(plotting_df_h05_unlinked_new, aes(x = entry_number, y = prop, fill = selfing)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", fill = "Site type") +
  geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(DFE)) +
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), 
        strip.text = element_text(size = 15), plot.title = element_text(size = 25), 
        legend.title = element_text(size = 15), legend.text = element_text(size = 15),
        legend.position = "bottom") +
  expand_limits(y = c(0, 0.7)) +
  scale_x_continuous(breaks = c(1, seq(5, 9, by = 5), 11),
                     labels = c(1, seq(5, 9, by = 5), "11+"))
hdel05_unlinked_neutral <- ggplot(plotting_df_h05_unlinked_neutral, aes(x = entry_number, y = prop, fill = selfing)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", fill = "Site type") +
  geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  #facet_grid(rows = vars(DFE)) +
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 25), 
        strip.text = element_text(size = 15), plot.title = element_text(size = 15), 
        legend.title = element_text(size = 15), legend.text = element_text(size = 15),
        legend.position = "bottom") +
  expand_limits(y = c(0, 0.7)) +
  scale_x_continuous(breaks = c(1, seq(5, 9, by = 5), 11),
                     labels = c(1, seq(5, 9, by = 5), "11+"))

ggsave(paste0(figures_dir,"sfigure03.svg"), plot = hdel05_unlinked_figure, width = 8.5, height = 8.5, dpi = 150)


plotting_df_h05_linked <- plotting_df_05099 %>% 
mutate(
    SFS = case_when(
      SFS == "selected" ~ "selected, h=0.5",
      SFS == "neutral" ~ "neutral, h=0.5",
      TRUE ~ SFS
    ))
plotting_df_h025_linked <- plotting_df_h025 %>% 
mutate(
    SFS = case_when(
      SFS == "selected" ~ "selected, h=0.25",
      SFS == "neutral" ~ "neutral, h=0.25",
      TRUE ~ SFS
    ))
plotting_linked_dom_comparison <- 
  bind_rows(plotting_df_h05_linked,plotting_df_h025_linked) %>%
  distinct()  



dom_comp_fig_linked <- ggplot(plotting_linked_dom_comparison, aes(x = entry_number, y = prop, fill = factor(SFS))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", aes(group = interaction(SFS, selfing))) +
  labs(x = "Derived allele count", y = "Proportion of polymorphisms", fill = "Site type") +
  geom_errorbar(aes(ymin = prop - propsd, ymax = prop + propsd), position = position_dodge(width = 0.9)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +

  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
        axis.title.x=element_text(size=15), axis.title.y=element_text(size=15), 
        strip.text = element_text(size=15), plot.title= element_text(size=25), 
        legend.title = element_text(size=15), legend.text = element_text(size=15),
        legend.position = "bottom" ) +
  expand_limits(y=c(0,0.7)) +
  #important addition to make x axis more readable
  scale_x_continuous(breaks = c(1, seq(5, 9, by = 5), 11),
                       labels = c(1, seq(5, 9, by = 5), "11+")) +
  guides(fill = guide_legend(nrow = 2))
ggsave(paste0(figures_dir,"sfigure08.svg"), plot = dom_comp_fig_linked, width = 8.5, height = 10, dpi = 150)
