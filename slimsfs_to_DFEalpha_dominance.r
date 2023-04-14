# An R script to convert slim simulation SFS output to DFE alpha's input format
#TODO fix many hard coded directory paths. 
#Make sure input and output directories go where you want them to

#initialize. Eventually can make the path an argument or at least relative. 
rm(list=ls())
library(tidyverse)
path_to_files <- "/nas/longleaf/home/adaigle/work/dominance_inputsandoutputs/SFS/0/"
path_to_DFESelfing <- "/nas/longleaf/home/adaigle/work/dominance_inputsandoutputs/DFE_alpha_dom_input_0/"
path_to_dfe_alpha_output <- "/nas/longleaf/home/adaigle/work/dominance_inputsandoutputs/DFE_alpha_dom_output_0/"
path_to_grapes_current_input <- "/nas/longleaf/home/adaigle/work/dominance_inputsandoutputs/grapes_dom_input_0/"

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

#plot row as DF
#barplot(unlist(eqm_selfing99_DFE1_neu_100_m1[1, 3:101]))
#barplot(unlist(eqm_selfing99_DFE1_neu_100_m1[2, 3:101]))
#barplot(unlist(eqm_selfing99_DFE1_neu_100_m1[3, 3:101]))
#barplot(unlist(eqm_selfing99_DFE1_neu_100_m1[4, 3:101]))
#barplot(unlist(eqm_selfing99_DFE1_neu_100_m1[5, 3:101]))
#
#
#barplot(unlist(eqm_selfing99_DFE1_sel_100_m2[1, 3:101]))
#barplot(unlist(eqm_selfing99_DFE1_sel_100_m2[2, 3:101]))
#barplot(unlist(eqm_selfing99_DFE1_sel_100_m2[3, 3:101]))
#barplot(unlist(eqm_selfing99_DFE1_sel_100_m2[4, 3:101]))
#barplot(unlist(eqm_selfing99_DFE1_sel_100_m2[5, 3:101]))
#
#barplot(unlist(eqm_selfing99_DFE1_neu_100_m1[1, 50:102]))
#barplot(unlist(eqm_selfing99_DFE1_neu_100_m1[2, 50:102]))
#barplot(unlist(eqm_selfing99_DFE1_neu_100_m1[3, 50:102]))
#barplot(unlist(eqm_selfing99_DFE1_neu_100_m1[4, 50:102]))
#barplot(unlist(eqm_selfing99_DFE1_neu_100_m1[5, 50:102]))
#
#barplot(unlist(eqm_selfing99_DFE1_sel_100_m2[1, 50:102]))
#barplot(unlist(eqm_selfing99_DFE1_sel_100_m2[2, 50:102]))
#barplot(unlist(eqm_selfing99_DFE1_sel_100_m2[3, 50:102]))
#barplot(unlist(eqm_selfing99_DFE1_sel_100_m2[4, 50:102]))
#barplot(unlist(eqm_selfing99_DFE1_sel_100_m2[5, 50:102]))

#complex fraction change stuff
#cell to frac
#neu2frac <- unlist(eqm_selfing99_DFE1_neu_100_m1[3, 2:102])/sum(unlist(eqm_selfing99_DFE1_neu_100_m1[3, 2:102]))
#sel2frac <- unlist(eqm_selfing99_DFE1_sel_100_m2[3, 2:102])/sum(unlist(eqm_selfing99_DFE1_sel_100_m2[3, 2:102]))
#barplot(sel2frac-neu2frac)
#div <- sel2frac/neu2frac
#div[is.infinite(div)] <- 0
#barplot(div)
#
#neu2frac <- unlist(eqm_selfing99_DFE1_neu_100_m1[4, 2:102])/sum(unlist(eqm_selfing99_DFE1_neu_100_m1[4, 2:102]))
#sel2frac <- unlist(eqm_selfing99_DFE1_sel_100_m2[4, 2:102])/sum(unlist(eqm_selfing99_DFE1_sel_100_m2[4, 2:102]))
#barplot(sel2frac-neu2frac)
#div <- sel2frac/neu2frac
#div[is.infinite(div)] <- 0
#barplot(div)
#
#neu2frac <- unlist(eqm_selfing99_DFE1_neu_100_m1[5, 2:102])/sum(unlist(eqm_selfing99_DFE1_neu_100_m1[5, 2:102]))
#sel2frac <- unlist(eqm_selfing99_DFE1_sel_100_m2[5, 2:102])/sum(unlist(eqm_selfing99_DFE1_sel_100_m2[5, 2:102]))
#barplot(sel2frac-neu2frac)
#div <- sel2frac/neu2frac
#div[is.infinite(div)] <- 0
#barplot(div)
#
##dfe2
#barplot(unlist(eqm_selfing99_DFE2_neu_100_m1[1, 3:101]))
#barplot(unlist(eqm_selfing99_DFE2_neu_100_m1[2, 3:101]))
#barplot(unlist(eqm_selfing99_DFE2_neu_100_m1[3, 3:101]))
#barplot(unlist(eqm_selfing99_DFE2_neu_100_m1[4, 3:101]))
#barplot(unlist(eqm_selfing99_DFE2_neu_100_m1[5, 3:101]))
#
#barplot(unlist(eqm_selfing99_DFE2_sel_100_m2[1, 3:101]))
#barplot(unlist(eqm_selfing99_DFE2_sel_100_m2[2, 3:101]))
#barplot(unlist(eqm_selfing99_DFE2_sel_100_m2[3, 3:101]))
#barplot(unlist(eqm_selfing99_DFE2_sel_100_m2[4, 3:101]))
#barplot(unlist(eqm_selfing99_DFE2_sel_100_m2[5, 3:101]))
#
## set up the plot layout
#par(mfrow = c(2, 3), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0), xpd = TRUE)
#for (bp in list(bp1, bp2, bp3, bp4, bp5)) {
#  axis(1, at = bp, labels = FALSE, tick = FALSE)
#  axis(2, ylim = c(0, max(bp)), las = 1)
#  box()
#}
##now for each experiment, I need to loop through each output, 
##pasting neu and sel together in a file
#
##code to get output1 from DFE1_nue
#
DFE_list <- c("DFE1", "DFE2", "DFE3")

#this assumes all files in SFS folder have same number and name of replicates
#if this assumption is violated the code will need to get more complex

replicates <- get(combined_df_names_list[1])$filename

#clears out script to run dfe alpha 
#write("", file = paste(path_to_DFESelfing, "run_dfealpha", sep = ""))


dfealpha_sfs <- function(x) {
    for(y in replicates) {
    df <- data.frame(Map(c,
        get(combined_df_names_list[grepl(paste(x, "_hdel_0_25_sel", sep = ""), 
            combined_df_names_list)])[get(combined_df_names_list[grepl(paste(x, "_hdel_0_25_sel", sep = ""), 
            combined_df_names_list)])$filename == y, ],
        get(combined_df_names_list[grepl(paste(x, "_hdel_0_25_neu", sep = ""), 
            combined_df_names_list)])[get(combined_df_names_list[grepl(paste(x, "_hdel_0_25_neu", sep = ""), 
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


#path_to_polyDFE_current_input <- "/nas/longleaf/home/adaigle/DFESelfing/polyDFE_input/"

# this function creates input sfs files for polydfe 
# we won't be using it so its commented out for now
#will add commands in a bit
#polydfe_sfs <- function(x) {
#for(y in replicates) {
#    df <- data.frame(Map(c,
#        get(combined_df_names_list[grepl(paste(x, "_neu", sep = ""), 
#            combined_df_names_list)])[get(combined_df_names_list[grepl(paste(x, "_neu", sep = ""), 
#            combined_df_names_list)])$filename == y, ],
#        get(combined_df_names_list[grepl(paste(x, "_sel", sep = ""), 
#            combined_df_names_list)])[get(combined_df_names_list[grepl(paste(x, "_sel", sep = ""), 
#            combined_df_names_list)])$filename == y, ]))
#    #path and name of final file
#    filepath <- paste(path_to_DFESelfing, x, y, sep = "")
#    df_stripped <- df[2:102]
#    names(df_stripped) <- NULL
#    #write to fle with proper header structure. Assumes 100 alleles were chosen
#    write(1 1 100, file = filepath)
#    write(100, file = filepath, append = TRUE)
#    write.table(df_stripped, row.names = FALSE, quote = FALSE, 
#        file = filepath, append = TRUE)
#
#}}


# this function creates input sfs files for grapes
#will add commands in a bit
grapes_sfs <- function(x) {
for(y in replicates) {
    neudf <- data.frame(Map(c,
        get(combined_df_names_list[grepl(paste(x, "_hdel_0_25_neu", sep = ""), 
            combined_df_names_list)])[get(combined_df_names_list[grepl(paste(x, "_hdel_0_25_neu", sep = ""), 
            combined_df_names_list)])$filename == y, ]))
    seldf <- data.frame(Map(c,
        get(combined_df_names_list[grepl(paste(x, "_hdel_0_25_sel", sep = ""), 
            combined_df_names_list)])[get(combined_df_names_list[grepl(paste(x, "_hdel_0_25_sel", sep = ""), 
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
