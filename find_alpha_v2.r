library(tidyverse)
dirs <- list.dirs(path = "/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/low_positive")
dirs <- dirs[ grepl("eqm", dirs) ]
string_pattern <- c("m4", "m2", "m1")

tidy_summary_table <- tibble(
    name = list.files(path = dirs, pattern = "^output..fixed"),
    fullpath = paste0(dirs, "/", name),
    selfing = str_extract(fullpath, "(?<=eqm_selfing)\\d+"),
    DFE = str_extract(fullpath, "(DFE)\\d+"),
    output = str_extract(fullpath, "(?<=output)\\d+"),
    data = lapply(fullpath, readLines), # must do this as rows have different column numbers
    filtered_lines = lapply(data, function(x) x[grep(paste(string_pattern, collapse="|"), x)])
    #pis = unlist(lapply(data, function(z) sum(2*(z$V9/100) *(1-(z$V9/100)))/total_neu_mutations)),
    #B = pis/theta,
    #empirical_Ne = B*5000
) %>% 
mutate(data = NULL, m1 = NULL, m2 = NULL, m4 = NULL)


# have to do some gymnastics now to turn my lists to df's, filter out muts fixed before 
tidy_summary_table <- tidy_summary_table %>% 
    mutate(newdata = lapply(filtered_lines, as.data.frame))

tidy_summary_table <- tidy_summary_table %>% 
    mutate(newdata2= lapply(newdata, function(x) separate(x, 1, into = c("id1", "id2", "mutation_class", "id3", "s", "h", "p1", "start_gen", "fix_gen"), sep = " ")
    ))

tidy_summary_table <- tidy_summary_table %>% mutate(newdata3 = lapply(newdata2, function(x) x %>% filter(as.numeric(fix_gen) > 50000))) %>%
    mutate(m1 = lapply(newdata3, function(df) sum(grepl("m1", df$mutation_class))), # quick way to count m1 rows 
    m2 = lapply(newdata3, function(df) sum(grepl("m2", df$mutation_class))),
    m4 = lapply(newdata3, function(df) sum(grepl("m4", df$mutation_class))),
    m1_count = unlist(lapply(m1, sum)), # quick way to count m1 rows  pt 2
    m2_count = unlist(lapply(m2, sum)),
    m4_count = unlist(lapply(m4, sum)),
    neutral_m4 = lapply(newdata3, function(x) x %>% filter((mutation_class == "m4"), 2*5000*as.numeric(s) < 5 )),
    neutral_m4_count = unlist(lapply(neutral_m4, function(x) nrow(x))),
    alpha = (as.numeric(m4_count) - neutral_m4_count)/ (as.numeric(m2_count) + as.numeric(m4_count))
)

tidy_summary_table <- tidy_summary_table %>% mutate(filtered_lines = NULL, newdata = NULL,newdata1 = NULL,newdata2 = NULL,newdata3 = NULL, m1 = NULL, m2 = NULL, m4 = NULL, neutral_m4 = NULL)

write.table(tidy_summary_table, file = "/nas/longleaf/home/adaigle/DFESelfing/alpha_list_v3.txt", 
    quote = FALSE, row.names = FALSE)
