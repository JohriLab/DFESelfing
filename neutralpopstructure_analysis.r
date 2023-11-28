# the purpose of this script is to make a table comparing the metapopulation 
# and subpopulation nucleotide diversity and FST to theoretical expectations.
library(tidyverse)

data_dir <- "/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/neutralpopstructure2/"

sub_dirs <- dir(data_dir)
#sub_dirs <- c("Nem1_eqm_selfing0", "Nem2_eqm_selfing0", "Nem01_eqm_selfing99")

parameter_table <- tibble(
    Nem= c(10,2,1,0.5,0.1),
    m = c(0.01016, 0.00216, 0.00116, 0.00066, 0.00026),
    Ne_sub = c(984, 926, 862, 785, 385),    
) 

summary_table <- tibble(output_paths = paste0(data_dir, sub_dirs, "/output", rep(1:5, each = length(sub_dirs))), 
    selfing = as.numeric(sub(".*selfing([0-9]+).*", "\\1", output_paths)),
    Nem = as.numeric(sub(".*Nem([[:digit:]]+).*", "\\1", output_paths)),
    Nem01 = as.numeric(sub(".*Nem(01).*", "\\1", output_paths)),
    Nem05 = as.numeric(sub(".*Nem(05).*", "\\1", output_paths)),
    meta_pi = unlist(lapply(output_paths, function(x) as.numeric(readLines(paste0(x, "_pi.txt"))))),
    subpop_pi = unlist(lapply(output_paths, function(x) as.numeric(readLines(paste0(x, "_pi_p1.txt"))))),
    fst = 1-(subpop_pi/meta_pi)
    ) %>%
    replace_na(list(Nem01 = 0, Nem05 = 0)) %>%  # Replace NA with 0 for specified columns
    # fix Nem values for Nem01 and Nem05
    mutate(Nem = ifelse(Nem01 == 1, 0.1, Nem),
        Nem = ifelse(Nem05 == 5, 0.5, Nem)) %>%
    select(!c(Nem01,Nem05)) 

#summary_table$Nem2 <- ifelse(summary_table$Nem01 == 1, 0.1, "test")

selfing_Ne <- function(S) {
    F = S*0.01/(2-S*0.01)
    Ne = 5000/(1+F)
    return(Ne)
}

csv <- right_join(summary_table, parameter_table, by="Nem") %>% 
    mutate(
        metapopsize_self = ((2*5*(Ne_sub*(selfing_Ne(selfing)/5000))) + (16/(10*m))) / 2,
        meta_pi_exp = 4*5000*3.3E-7 * (metapopsize_self)/5000,
        subpop_pi_exp = 4*5*Ne_sub*3.3E-7 * (selfing_Ne(selfing)/5000)
    ) %>% mutate(fst_exp = 1-(subpop_pi_exp/meta_pi_exp)) %>%
    group_by(Nem,selfing) %>%
    summarize(across(where(is.numeric), list(avg = mean))) %>%
    select(Nem, selfing, meta_pi_exp_avg, meta_pi_avg, subpop_pi_exp_avg, subpop_pi_avg, fst_exp_avg, fst_avg)

write.csv(csv, file = "/nas/longleaf/home/adaigle/DFESelfing/neutralpopstructure.csv", quote = F, row.names=F)
