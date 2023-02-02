# An R script to convert slim simulation SFS output to DFE alpha's input format

#initialize. Eventually can make the path an argument or at least relative. 
rm(list=ls())
library(tidyverse)
path_to_files <- "/nas/longleaf/home/adaigle/DFESelfing/SFS/"
#total neutral sites is 187500
neutral_sites <- 187500
#total selected sites are 562500
selected_sites <- 562500

#read in slim sfs and fixed counts to a list of dataframes
count_sfs_names <- list.files(path = path_to_files, pattern = "count.sfs")
count_sfs_df_list <- lapply(
    paste(path_to_files,count_sfs_names,sep=""), 
    function(x) read.table(x, header=TRUE))

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

#now for each experiment, I need to loop through each output, 
#pasting neu and sel together in a file

#code to get output1 from DFE1_nue
# can make loops to move through DFEs, then through outputs 
# can hard code lists for now?
neu <- c("DFE1_neu", "DFE2_neu")
sel <- c("DFE1_sel", "DFE2_sel")

DFE_list <- c("DFE1", "DFE2")

paste(x, "_neu", sep = "")
paste(x, "_sel", sep = "")

paste("/DFE_alpha_SFS_format/", x, y, sep='')
replicates <- get(combined_df_names_list[grepl(paste(x, "_neu", sep = ""), combined_df_names_list)])$filename
dfealpha_sfs <- function(x) {
for(y in replicates) {
    df <- data.frame(Map(c,
        get(combined_df_names_list[grepl(paste(x, "_neu", sep = ""), combined_df_names_list)])[get(combined_df_names_list[grepl(paste(x, "_neu", sep = ""), combined_df_names_list)])$filename == y, ],
        get(combined_df_names_list[grepl(paste(x, "_sel", sep = ""), combined_df_names_list)])[get(combined_df_names_list[grepl(paste(x, "_sel", sep = ""), combined_df_names_list)])$filename == y, ]))
    #path and name of final file
    filepath <- paste(path_to_files, paste("/DFE_alpha_SFS_format/", x, y, sep=""), sep = "")
    df_stripped <- df[2:102]
    names(df_stripped) <- NULL
    #write to fle with proper header structure. Assumes 100 alleles were chosen
    write(1, file = filepath)
    write(100, file = filepath, append = TRUE)
    write.table(df_stripped, row.names = FALSE, quote = FALSE, 
        file = filepath, append = TRUE)
}

}

lapply(DFE_list, dfealpha_sfs)

c(1,100,dfealphaformat[3][[1]])

#need to write code to make unique config files for each run, so that DFE alpha does not overwrite them

df <- data.frame(Map(c,
    get(combined_df_names_list[grepl("DFE1_neu", combined_df_names_list)])[get(combined_df_names_list[grepl("DFE1_neu", combined_df_names_list)])$filename == 'output1.txt', ],
    get(combined_df_names_list[grepl("DFE1_sel", combined_df_names_list)])[get(combined_df_names_list[grepl("DFE1_sel", combined_df_names_list)])$filename == 'output1.txt', ]))


path <- paste(path_to_files, "/DFE_alpha_SFS_format/DFE1SFS_dfealpha.txt", sep = "")
df_stripped <- df[2:102]
names(df_stripped) <- NULL


write(1, file = path)
write(100, file = path, append = TRUE)
write.table(df_stripped, row.names = FALSE, quote = FALSE, 
    file = path, append = TRUE)