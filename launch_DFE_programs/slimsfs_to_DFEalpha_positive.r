# An R script to convert slim simulation SFS output to DFE alpha and GRAPES input formats
# Prepares slurm scripts for launching from specific directories. 
# A bit more complex than the others as it also sets up folded SFS analyses. 

rm(list=ls())
library(tidyverse)
selfing_levels <- c("0", "50", "80", "90", "95", "99")
#selfing_levels <- c("0")

fold_sfs <- function(unfolded_df) {
    # set up to fold a SFS in the format of this script for 100 individuals.
    # can be made generic if needed
    folded_df <- data.frame(Column1 = rep(NA, 5))
    for (i in 1:49) {
        #fold SFS
        folded_df[1] <- unfolded_df[1]
        folded_df[2] <- unfolded_df[2]
        folded_df[i+2] <- unfolded_df[i+2] + unfolded_df[100-i+2]
    }
        # entry 50
        folded_df[52] <- unfolded_df[52]
        # fixed stuff
        folded_df[53] <- unfolded_df[102]
        return(folded_df) 
}


for(selfing_level in selfing_levels) {
print(selfing_level)
main_dir <- "/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/low_positive"
path_to_files <- paste0(main_dir, "/SFS/")
path_to_DFESelfing <- paste0(main_dir, "/dfe_results/dfealpha/DFE_alpha_input_", selfing_level, "/")
print(path_to_DFESelfing)
path_to_DFEalpha_autofold <- paste0(main_dir, "/dfe_results/dfealpha/DFE_alpha_autofold_input_", selfing_level, "/")
print(path_to_DFEalpha_autofold)
#path_to_DFEalpha_manualfold <- paste0(main_dir, "/dfe_results/dfealpha/DFE_alpha_manualfold_input_", selfing_level, "/")
#print(path_to_DFEalpha_manualfold)

path_to_dfe_alpha_output <- paste0(main_dir, "/dfe_results/dfealpha/DFE_alpha_output_", selfing_level, "/")
path_to_dfe_alpha_autofold_output <- paste0(main_dir, "/dfe_results/dfealpha/DFE_alpha_autofold_output_", selfing_level, "/")
#path_to_dfe_alpha_manualfold_output <- paste0(main_dir, "/dfe_results/dfealpha/DFE_alpha_manualfold_output_", selfing_level, "/")

path_to_grapes_current_input <- paste0(main_dir, "/dfe_results/grapes/grapes_input_", selfing_level, "/")
grapes_output <- paste0(main_dir, "/dfe_results/grapes/grapes_output_", selfing_level, "/")

path_to_grapes_nodivparam_current_input <- paste0(main_dir, "/dfe_results/grapes/grapes_nodivparam_input_", selfing_level, "/")
grapes_nodivparam_output <- paste0(main_dir, "/dfe_results/grapes/grapes_nodivparam_output_", selfing_level, "/")

path_to_grapes_folded_current_input <- paste0(main_dir, "/dfe_results/grapes/grapes_folded_input_", selfing_level, "/")
grapes_folded_output <- paste0(main_dir, "/dfe_results/grapes/grapes_folded_output_", selfing_level, "/")

dir.create(file.path(paste0(main_dir, "/dfe_results/")))
dir.create(file.path(paste0(main_dir, "/dfe_results/dfealpha/")))
dir.create(file.path(paste0(main_dir, "/dfe_results/grapes")))

dir.create(file.path(path_to_DFESelfing))
dir.create(file.path(path_to_dfe_alpha_output))
dir.create(file.path(path_to_DFEalpha_autofold))
dir.create(file.path(path_to_dfe_alpha_autofold_output))
#dir.create(file.path(path_to_DFEalpha_manualfold))
#dir.create(file.path(path_to_dfe_alpha_manualfold_output))

dir.create(file.path(path_to_grapes_current_input))
dir.create(file.path(grapes_output))

dir.create(file.path(path_to_grapes_nodivparam_current_input))
dir.create(file.path(grapes_nodivparam_output))

dir.create(file.path(path_to_grapes_folded_current_input))
dir.create(file.path(grapes_folded_output))

neutral_sites <- 187500
selected_sites <- 562500

eqm <- paste0("eqm_selfing", selfing_level)
#read in slim sfs and fixed counts to a list of dataframes

count_sfs_names <- list.files(path = path_to_files, pattern = "count.sfs")
count_sfs_names <- count_sfs_names[grepl(eqm, count_sfs_names) & grepl("count.sfs", count_sfs_names)]
count_sfs_df_list <- lapply(
    paste(path_to_files,count_sfs_names,sep=""), 
    function(x) read.table(x, header=TRUE))


#--------------------------------
#stripping _count.sfs from file names to match fixed_sfs header, to name list
names(count_sfs_df_list) <- lapply(count_sfs_names,
    function(x)
        sub("_count.sfs.*","",x)
)

fixed_sfs_names <- list.files(path = path_to_files, pattern = "count.fixed")
fixed_sfs_names <- fixed_sfs_names[grepl(eqm, fixed_sfs_names) & grepl("count.fixed", fixed_sfs_names)]

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

#make folded SFS 
for(x in combined_df_names_list[grepl("neu", combined_df_names_list)]) {
    assign(paste0(x, "_folded"), 
      fold_sfs(get(x)))
}

#Get sel dfs, subtract sum of all rows from num sel sites
for(x in combined_df_names_list[grepl("sel_", combined_df_names_list)]) {
    assign(x, 
    add_column(get(x), selected_sites - rowSums(get(x)[2:101]), .before = 2))
}

for(x in combined_df_names_list[grepl("sel_", combined_df_names_list)]) {
    assign(paste0(x, "_folded"), 
      fold_sfs(get(x)))
}


DFE_list <- c("DFE1", "DFE2", "DFE3")

#this assumes all files in SFS folder have same number and name of replicates
#if this assumption is violated the code will need to get more complex

replicates <- get(combined_df_names_list[1])$filename
replicates <- replicates[!grepl("m1", replicates)]



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

    slurmpath <- paste0(path_to_DFESelfing, "../../dfealpha_slurm", selfing_level, ".sh")
    write("#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 00-02:00:00
#SBATCH -n 1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adaigle@email.unc.edu

", file = slurmpath)
    write("cd ~/dfe-alpha-release-2.16/", file = slurmpath, append = TRUE)
    write(paste0("bash ", path_to_DFESelfing, "run_dfealpha"), file = slurmpath, append = TRUE)

}}

lapply(DFE_list, dfealpha_sfs)

dfealpha_sfs_autofold <- function(x) {
#key difference from other cases: Fold = 1 instead of 0
    for(y in replicates) {
    df <- data.frame(Map(c,
        get(combined_df_names_list[grepl(paste(x, "_sel", sep = ""), 
            combined_df_names_list)])[get(combined_df_names_list[grepl(paste(x, "_sel", sep = ""), 
            combined_df_names_list)])$filename == y, ],
        get(combined_df_names_list[grepl(paste(x, "_neu", sep = ""), 
            combined_df_names_list)])[get(combined_df_names_list[grepl(paste(x, "_neu", sep = ""), 
            combined_df_names_list)])$filename == y, ]))
    #path and name of final file
    filepath <- paste(path_to_DFEalpha_autofold, x, y, sep = "")
    neuconfigpath <- paste(path_to_DFEalpha_autofold, "neu_config_", x, y, sep = "")
    selconfigpath <- paste(path_to_DFEalpha_autofold, "sel_config_", x, y, sep = "")
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
fold 1
epochs 2
search_n2 1
t2_variable 1
t2 50", file = neuconfigpath )
    write(paste("sfs_input_file", filepath), 
        file = neuconfigpath, append = TRUE)
    write(
        paste("est_dfe_results_dir ", path_to_dfe_alpha_autofold_output, 
            x, y, "_neutral", sep = ""), 
        file = neuconfigpath, append = TRUE)
    write("data_path_1 /nas/longleaf/home/adaigle/work/johri_elegans/data
site_class 1
fold 1
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
        paste("est_dfe_results_dir ", path_to_dfe_alpha_autofold_output, 
            x, y, "_selected", sep = ""), 
        file = selconfigpath, append = TRUE)
    write(
        paste("est_dfe_demography_results_file ", path_to_dfe_alpha_autofold_output, 
            x, y, "_neutral", "/est_dfe.out", sep = ""), 
        file = selconfigpath, append = TRUE)

    #create a script to run dfe_alpha on all configs on the command line
    #to run: "bash run_dfealpha"
    write(paste("./est_dfe -c ", neuconfigpath, sep = ""), 
        file = paste(path_to_DFEalpha_autofold, "run_dfealpha", sep = ""), append = TRUE)
    write(paste("./est_dfe -c ", selconfigpath, sep = ""), 
        file = paste(path_to_DFEalpha_autofold, "run_dfealpha", sep = ""), append = TRUE)

    slurmpath <- paste0(path_to_DFEalpha_autofold, "../../dfealpha_autofold_slurm", selfing_level, ".sh")
    write("#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 00-02:00:00
#SBATCH -n 1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adaigle@email.unc.edu

", file = slurmpath)
    write("cd ~/dfe-alpha-release-2.16/", file = slurmpath, append = TRUE)
    write(paste0("bash ", path_to_DFEalpha_autofold, "run_dfealpha"), file = slurmpath, append = TRUE)

}}

lapply(DFE_list, dfealpha_sfs_autofold)

dfealpha_sfs_manualfold <- function(x) {
    for(y in replicates) {
    df <- data.frame(Map(c,
        get(combined_df_names_list[grepl(paste(x, "_sel", sep = ""), 
            combined_df_names_list)])[get(combined_df_names_list[grepl(paste(x, "_sel", sep = ""), 
            combined_df_names_list)])$filename == y, ],
        get(combined_df_names_list[grepl(paste(x, "_neu", sep = ""), 
            combined_df_names_list)])[get(combined_df_names_list[grepl(paste(x, "_neu", sep = ""), 
            combined_df_names_list)])$filename == y, ]))
    #path and name of final file
    filepath <- paste(path_to_DFEalpha_manualfold, x, y, sep = "")
    neuconfigpath <- paste(path_to_DFEalpha_manualfold, "neu_config_", x, y, sep = "")
    selconfigpath <- paste(path_to_DFEalpha_manualfold, "sel_config_", x, y, sep = "")
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
        paste("est_dfe_results_dir ", path_to_dfe_alpha_manualfold_output, 
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
        paste("est_dfe_results_dir ", path_to_dfe_alpha_manualfold_output, 
            x, y, "_selected", sep = ""), 
        file = selconfigpath, append = TRUE)
    write(
        paste("est_dfe_demography_results_file ", path_to_dfe_alpha_manualfold_output, 
            x, y, "_neutral", "/est_dfe.out", sep = ""), 
        file = selconfigpath, append = TRUE)

    #create a script to run dfe_alpha on all configs on the command line
    #to run: "bash run_dfealpha"
    write(paste("./est_dfe -c ", neuconfigpath, sep = ""), 
        file = paste(path_to_DFEalpha_manualfold, "run_dfealpha", sep = ""), append = TRUE)
    write(paste("./est_dfe -c ", selconfigpath, sep = ""), 
        file = paste(path_to_DFEalpha_manualfold, "run_dfealpha", sep = ""), append = TRUE)

    slurmpath <- paste0(path_to_DFEalpha_manualfold, "../../dfealpha_slurm", selfing_level, ".sh")
    write("#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 00-02:00:00
#SBATCH -n 1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adaigle@email.unc.edu

", file = slurmpath)
    write("cd ~/dfe-alpha-release-2.16/", file = slurmpath, append = TRUE)
    write(paste0("bash ", path_to_DFESelfing, "run_dfealpha"), file = slurmpath, append = TRUE)

}}

#lapply(DFE_list, dfealpha_sfs_manualfold)

grapes_sfs <- function(x) {
    #Grapes GammaExpo, default parameters, fullDFE
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

    slurmpath <- paste0(path_to_DFESelfing, "../../grapes_slurm", selfing_level, ".sh")
    write("#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 00-02:00:00
#SBATCH --mem=16g
#SBATCH -n 15
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adaigle@email.unc.edu

", file = slurmpath)
    write("cd ~/grapes/", file = slurmpath, append = TRUE)
    write(paste0("ls ", path_to_grapes_current_input, " | parallel ./grapes -in ", 
        path_to_grapes_current_input, "{1} -out ", grapes_output, "{1}.csv -model GammaExpo "), 
        file = slurmpath, append = TRUE)
}}

lapply(DFE_list, grapes_sfs)


grapes_sfs_nodivparam <- function(x) {
    #Grapes GammaExpo, no_div_param , fullDFE
    #does this account for adaptive polymorphisms in the data?
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
    filepath <- paste(path_to_grapes_nodivparam_current_input, x, y, sep = "")
    neudf_stripped <- neudf[3:101]
    seldf_stripped <- cbind(seldf[1],100,selected_sites,seldf[3:101])
    #write to fle with proper header structure. Assumes 100 alleles were chosen
    grapes_format <- as.data.frame(c(seldf_stripped,neutral_sites,neudf_stripped, selected_sites, seldf[102], neutral_sites, neudf[102]))
    names(grapes_format) <- NULL
    write(paste(x, y), file = filepath)
    write("#unfolded", file = filepath, append = TRUE)
    write.table(grapes_format, row.names = FALSE, quote = FALSE, sep = "\t",
        file = filepath, append = TRUE)

    slurmpath <- paste0(path_to_DFESelfing, "../../grapes_slurm_nodivparam", selfing_level, ".sh")
    write("#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 00-02:00:00
#SBATCH --mem=16g
#SBATCH -n 15
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adaigle@email.unc.edu

", file = slurmpath)
    write("cd ~/grapes/", file = slurmpath, append = TRUE)
    write(paste0("ls ", path_to_grapes_nodivparam_current_input, " | parallel ./grapes -in ", 
        path_to_grapes_nodivparam_current_input, "{1} -out ", grapes_nodivparam_output, "{1}.csv -model GammaExpo -no_div_param "), 
        file = slurmpath, append = TRUE)
}}

lapply(DFE_list, grapes_sfs_nodivparam)

combined_df_names_list_folded <- c()
for (i in 1:length(combined_df_names_list)) {
    combined_df_names_list_folded[i] <- paste0(combined_df_names_list[i], "_folded")
}
grapes_sfs_folded <- function(x) {
    #grapes GammaExpo on a folded SFS 
for(y in replicates) {
    neudf <- data.frame(Map(c,
        get(combined_df_names_list_folded[grepl(paste(x, "_neu", sep = ""), 
            combined_df_names_list_folded)])[get(combined_df_names_list_folded[grepl(paste(x, "_neu", sep = ""), 
            combined_df_names_list_folded)])[[1]] == y, ]))
    seldf <- data.frame(Map(c,
        get(combined_df_names_list_folded[grepl(paste(x, "_sel", sep = ""), 
            combined_df_names_list_folded)])[get(combined_df_names_list_folded[grepl(paste(x, "_sel", sep = ""), 
            combined_df_names_list_folded)])[[1]] == y, ]))
    #path and name of final file
    filepath <- paste(path_to_grapes_folded_current_input, x, y, sep = "")
    neudf_stripped <- neudf[3:52]
    seldf_stripped <- cbind(seldf[1],100,selected_sites,seldf[3:52])
    #write to fle with proper header structure. Assumes 100 alleles were chosen
    grapes_format <- as.data.frame(c(seldf_stripped,neutral_sites,neudf_stripped, selected_sites, seldf[53], neutral_sites, neudf[53]))
    names(grapes_format) <- NULL
    write(paste(x, y), file = filepath)
    write.table(grapes_format, row.names = FALSE, quote = FALSE, sep = "\t",
        file = filepath, append = TRUE)

    slurmpath <- paste0(path_to_DFESelfing, "../../grapes_slurm_folded", selfing_level, ".sh")
    write("#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 00-02:00:00
#SBATCH --mem=16g
#SBATCH -n 15
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adaigle@email.unc.edu

", file = slurmpath)
    write("cd ~/grapes/", file = slurmpath, append = TRUE)
    write(paste0("ls ", path_to_grapes_folded_current_input, " | parallel ./grapes -in ", 
        path_to_grapes_folded_current_input, "{1} -out ", grapes_folded_output, "{1}.csv -model GammaExpo "), 
        file = slurmpath, append = TRUE)
}}

lapply(DFE_list, grapes_sfs_folded)

}
