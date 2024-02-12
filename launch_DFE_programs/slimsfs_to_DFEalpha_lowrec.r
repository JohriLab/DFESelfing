# An R script to convert slim simulation SFS output to DFE alpha and GRAPES input formats
# Prepares slurm scripts for launching from specific directories. 

rm(list=ls())
library(tidyverse)
path_to_files <- "/nas/longleaf/home/adaigle/work/lowrecombination_inputsandoutputs/SFS/"
path_to_DFESelfing <- "/nas/longleaf/home/adaigle/work/lowrecombination_inputsandoutputs/DFE_alpha_lowrecom/"
path_to_dfe_alpha_output <- "/nas/longleaf/home/adaigle/work/lowrecombination_inputsandoutputs/DFE_alpha_lowrecom_output/"
path_to_grapes_current_input <- "/nas/longleaf/home/adaigle/work/lowrecombination_inputsandoutputs/grapes_lowrecom_input/"

#total neutral sites is 187500
neutral_sites <- 187500
#total selected sites are 562500
selected_sites <- 562500

#read in slim sfs and fixed counts to a list of dataframes
count_sfs_names <- list.files(path = path_to_files, pattern = "count.sfs")
count_sfs_df_list <- lapply(
    paste(path_to_files,count_sfs_names,sep=""), 
    function(x) read.table(x, header=TRUE))

#new table method. Not currently in use 
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

DFE_list <- c("DFE1", "DFE2", "DFE3")

#this assumes all files in SFS folder have same number and name of replicates
#if this assumption is violated the code will need to get more complex

replicates <- get(combined_df_names_list[1])$filename

#clears out script to run dfe alpha 
#write("", file = paste(path_to_DFESelfing, "run_dfealpha", sep = ""))



dfealpha_sfs <- function(x) {
    for(y in replicates) {
    df <- data.frame(Map(c,
        get(combined_df_names_list[grepl(paste(x, "_sel", sep = ""), 
            combined_df_names_list)])[get(combined_df_names_list[grepl(paste(x, "_sel", sep = ""), 
            combined_df_names_list)])$filename == y, ],
        get(combined_df_names_list[grepl(paste(x, "_neu", sep = ""), 
            combined_df_names_list)])[get(combined_df_names_list[grepl(paste(x, "_neu", sep = ""), 
            combined_df_names_list)])$filename == y, ]))
    #path and name of final file
    filepath <- paste(path_to_DFESelfing, x, y, sep = "")
    neuconfigpath <- paste(path_to_DFESelfing, "neu_config_", x, y, sep = "")
    selconfigpath <- paste(path_to_DFESelfing, "sel_config_", x, y, sep = "")
    df_stripped <- df[2:102]
    names(df_stripped) <- NULL
    #write to fle with proper header structure. Assumes 100 alleles were chosen
    write(1, file = filepath)
    write(100, file = filepath, append = TRUE)
    write.table(df_stripped, row.names = FALSE, quote = FALSE, 
       file = filepath, append = TRUE)#
    #write neutral and selected config files
    write("data_path_1 /nas/longleaf/home/adaigle/work/johri_elegans/data
site_class 0
fold 0
epochs 2
search_n2 1
t2_variable 1
t2 50", file = neuconfigpath )
    write(paste("sfs_input_file", filepath), 
        file = neuconfigpath, append = TRUE)
    write(
        paste("est_dfe_results_dir ", path_to_dfe_alpha_output, 
            x, y, "_neutral", sep = ""), 
        file = neuconfigpath, append = TRUE)
    write("data_path_1 /nas/longleaf/home/adaigle/work/johri_elegans/data
site_class 1
fold 0
epochs 2
mean_s_variable 1
mean_s -0.01
beta_variable 1
beta 0.5
p_additional 0
s_additional 0 ", file = selconfigpath)
    write(paste("sfs_input_file ", filepath, sep = ""), 
        file = selconfigpath, append = TRUE)
    write(
        paste("est_dfe_results_dir ", path_to_dfe_alpha_output, 
            x, y, "_selected", sep = ""), 
        file = selconfigpath, append = TRUE)
    write(
        paste("est_dfe_demography_results_file ", path_to_dfe_alpha_output, 
            x, y, "_neutral", "/est_dfe.out", sep = ""), 
        file = selconfigpath, append = TRUE)

    #create a script to run dfe_alpha on all configs on the command line
    #to run: "bash run_dfealpha"
    write(paste("./est_dfe -c ", neuconfigpath, sep = ""), 
        file = paste(path_to_DFESelfing, "run_dfealpha", sep = ""), append = TRUE)
    write(paste("./est_dfe -c ", selconfigpath, sep = ""), 
        file = paste(path_to_DFESelfing, "run_dfealpha", sep = ""), append = TRUE)
#return(assign(df_stripped, paste0(x,y)))
}}

lapply(DFE_list, dfealpha_sfs)


# this function creates input sfs files for grapes
#will add commands in a bit
grapes_sfs <- function(x) {
for(y in replicates) {
    neudf <- data.frame(Map(c,
        get(combined_df_names_list[grepl(paste(x, "_neu", sep = ""), 
            combined_df_names_list)])[get(combined_df_names_list[grepl(paste(x, "_neu", sep = ""), 
            combined_df_names_list)])$filename == y, ]))
    seldf <- data.frame(Map(c,
        get(combined_df_names_list[grepl(paste(x, "_sel", sep = ""), 
            combined_df_names_list)])[get(combined_df_names_list[grepl(paste(x, "_sel", sep = ""), 
            combined_df_names_list)])$filename == y, ]))
    #path and name of final file
    filepath <- paste(path_to_grapes_current_input, x, y, sep = "")
    neudf_stripped <- neudf[3:101]
    seldf_stripped <- cbind(seldf[1],100,selected_sites,seldf[3:101])
    #write to fle with proper header structure. Assumes 100 alleles were chosen
    grapes_format <- as.data.frame(c(seldf_stripped,neutral_sites,neudf_stripped, selected_sites, seldf[102], neutral_sites, neudf[102]))
    names(grapes_format) <- NULL
    write(paste(x, y), file = filepath)
    write("#unfolded", file = filepath, append = TRUE)
    write.table(grapes_format, row.names = FALSE, quote = FALSE, sep = "\t",
        file = filepath, append = TRUE)
}}

lapply(DFE_list, grapes_sfs)
