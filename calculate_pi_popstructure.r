# a script to calculate pi using m1 mutations
# in slim output files
library(tidyverse)

#calculation will be p = last_column/100
#total_neu_mutations = 187500
#sum(2*p*(1-p))/total_neu_mutations
total_neu_mutations <- 187500
total_sel_mutations <- 562500
#need to read in all m1_output.txt files in a way that preserves selfing%_DFE#_replicate structure
#get pi for each output
#average and std dev for each DFE and make table like I did for gamma and beta

#dirs <- dir(path = "/nas/longleaf/home/adaigle/work/johri_elegans/gammaDFE", pattern = "DFE")
getpi <- function(Nem) {
dirs <- list.dirs(path = paste0("/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/popstructure/", Nem), recursive = F)
dirs <- dirs[ grepl("eqm", dirs) ]
#dirs <- dirs[!grepl("Nem2_eqm_selfing99_DFE2", dirs)]

string_pattern <- c("m1")

tidy_summary_table <- tibble(
    name = list.files(path = dirs, pattern = "^output..txt"),
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
mutate(data = NULL, m2 = NULL, m4 = NULL)

# have to do some gymnastics now to turn my lists to df's, filter out muts fixed before 
tidy_summary_table <- tidy_summary_table %>% 
    mutate(newdata = lapply(filtered_lines, as.data.frame))

tidy_summary_table <- tidy_summary_table %>% 
    mutate(newdata2= lapply(newdata, function(x) separate(x, 1, into = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"), sep = " ")
    ))

tidy_summary_table_final <- tidy_summary_table %>% #mutate(newdata3 = lapply(newdata2, function(x) x %>% filter(as.numeric(fix_gen) > 50000))) %>%
    mutate(#m1 = lapply(newdata3, function(df) sum(grepl("m1", df$mutation_class))), # quick way to count m1 rows 
    pis = unlist(lapply(newdata2, function(z) sum(2*(as.numeric(z$V9)/100) *(1-(as.numeric(z$V9)/100)))/total_neu_mutations))
) %>% 
group_by(selfing,DFE) %>% 
summarize(across(where(is.numeric), list(avg = mean))) 

tidy_summary_table_final$Nem <- Nem
return(tidy_summary_table_final)}

selfing_Ne <- function(selfing_prop, census_size) {
    census_size/(1+(selfing_prop/(2-selfing_prop)))
}

theta <- 4*5000*3.3e-9*100
theta_50 <- 4*selfing_Ne(0.5, 5000)*3.3e-9*100
theta_99 <- 4*selfing_Ne(0.99, 5000)*3.3e-9*100

theta_nem2 <- 4*4629.63*3.3e-9*100
theta_50_nem2 <- 4*selfing_Ne(0.5, 4629.63)*3.3e-9*100
theta_99_nem2 <- 4*selfing_Ne(0.99, 4629.63)*3.3e-9*100

theta_nem1 <- 4*4310.345*3.3e-9*100
theta_50_nem1 <- 4*selfing_Ne(0.5, 4310.345)*3.3e-9*100
theta_99_nem1 <- 4*selfing_Ne(0.99, 4310.345)*3.3e-9*100

theta_nem01 <- 4*1923.077*3.3e-9*100
theta_50_nem01 <- 4*selfing_Ne(0.5, 1923.077)*3.3e-9*100
theta_99_nem01 <- 4*selfing_Ne(0.99, 1923.077)*3.3e-9*100
#tidy_summary_table <- tidy_summary_table %>% mutate(filtered_lines = NULL, newdata = NULL,newdata1 = NULL,newdata2 = NULL,newdata3 = NULL, m1 = NULL, m2 = NULL, m4 = NULL, neutral_m4 = NULL)

pitables <- lapply(c("Nem10", "Nem2", "Nem1", "Nem05", "Nem01"), getpi)
pitables_final <- bind_rows(pitables)
write.csv(pitables_final, file = "/nas/longleaf/home/adaigle/DFESelfing/popstructure_pis.csv", quote=F, row.names=F)
