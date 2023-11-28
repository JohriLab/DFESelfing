# the purpose of this script is to combine the outputs of the five 
#subpopulations into one output
library(tidyverse)
dir <- "/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/popstructure_uneven_2/"
outputs <- rep(paste0("output", 1:5),each=5)
pops <- paste0("p", 1:5)
string_pattern <- c("m1","m2")

root_dir <-"/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/popstructure_uneven_2/"
dirs <- paste0(root_dir, dir(root_dir))
dirs <- dirs[24:45]
#dirs <- dirs[!grepl("Nem05|Nem10", dirs)]
#dirs <- dir()
#dirs <- dir
for (dir in dirs) {
your_tibble <- tibble(
    filenames = paste0(outputs,pops,".txt"),
    outputs = outputs,
    pop = gsub("output(\\d+)p(\\d+)\\.txt", "\\2", filenames),
    fullpath = paste0(dir,"/", filenames),
    data = lapply(fullpath, readLines), # must do this as rows have different column numbers
    filtered_lines = lapply(data, function(x) x[grep(paste(string_pattern, collapse="|"), x)]),
    df = lapply(filtered_lines, function(x) read.table(text=x, sep = " "))
) %>% group_by(outputs, pop)


# Define a list of groups
groups <- unique(your_tibble$outputs)

# Function to update pop1 using other pops for a specific group
# For each pop, add up mutations of same label in other pops
# Then uniquify to ensure no mutations are duplicated
# Mutation freq should add up the same way in each iteration
update_pop <- function(group) {
  # Filter the data for the current group and pop1
  pop_df_list <- list()
  for (i in 2:5) {
    #print(i)
    p1 <- your_tibble %>%
      filter(outputs == group, pop == i) %>% 
      select(df)
  
      p1_df <- p1$df[[1]]
      unique_id <- p1_df$V2
  
    # Extract the filenames for other pops in the current group
    update <- your_tibble %>%
      filter(outputs == group, pop != i) %>% 
      select(df)
  
  
    list_of_dataframes <- update$df
    for (df in list_of_dataframes) {
    # Check if the unique ID exists in the current dataframe
    matching_rows <- df$V2 %in% unique_id
    #matching_rows <- df$V2 %in% unique_id
  
    # If there are matches, update column 9 in df1
    if (any(matching_rows)) {
      p1_df$V9[unique_id %in% df$V2] <- p1_df$V9[unique_id %in% df$V2] + df$V9[matching_rows]
  }}
    #print(length(p1_df))]
    #append()
    pop_df_list <- c(pop_df_list, list(p1_df))

}
  unique_df <- do.call(rbind, pop_df_list) %>% 
      arrange(V2) %>% mutate(V1=V2) %>% distinct()

  #ensure no mutations are duplicated
  #um(duplicated(unique_df$V2))
  return(unique_df)
}

update_pop1 <- function(group) {
  df <- your_tibble %>%
      filter(outputs == group) %>% 
      ungroup %>%
      select(df) %>% unnest(df)

  result <- df %>% 
    group_by(V2) %>% 
    summarize(
    V1 = first(V1),
    V3 = first(V3),
    V4 = first(V4),  # Add additional columns as needed
    V5 = first(V5),
    V6 = first(V6),
    V7 = first(V7),
    V8 = first(V8),
    V9 = sum(V9)
  )
}

#run on outputs 1-5
new_outputs <- lapply(groups, update_pop1)
names(new_outputs) <- groups

#save to output#.txt
lapply(names(new_outputs), function(name) {
  file_name <- paste0(dir, "/", name, ".txt")
  write.table(new_outputs[[name]], file = file_name, sep = " ",quote=FALSE, row.names = FALSE, col.names = FALSE)
})
}
