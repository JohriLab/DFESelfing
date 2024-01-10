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
dirs <- list.dirs(path = "/nas/longleaf/home/adaigle/work/johri_elegans/gammaDFE")
dirs <- dirs[ grepl("eqm", dirs) ]
#dirs <- dirs[!grepl("eqm_lowrec50_DFE1", dirs)]

string_pattern <- c("m2")

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
    pis = unlist(lapply(newdata2, function(z) sum(2*(as.numeric(z$V9)/100) *(1-(as.numeric(z$V9)/100)))/total_sel_mutations)),
    p = unlist(lapply(newdata2, function(z) sum((as.numeric(z$V9)/100) /length(z$V9)))),
    #B = pis/theta,
    #empirical_Ne = B*5000
    #m2 = lapply(newdata3, function(df) sum(grepl("m2", df$mutation_class))),
    #m4 = lapply(newdata3, function(df) sum(grepl("m4", df$mutation_class))),
    #m1_count = unlist(lapply(m1, sum)), # quick way to count m1 rows  pt 2
    #m2_count = unlist(lapply(m2, sum)),
    #m4_count = unlist(lapply(m4, sum)),
    #neutral_m4 = lapply(newdata3, function(x) x %>% filter((mutation_class == "m4"), 2*5000*as.numeric(s) < 5 )),
    #neutral_m4_count = unlist(lapply(neutral_m4, function(x) nrow(x))),
    #alpha = (as.numeric(m4_count) - neutral_m4_count)/ (as.numeric(m2_count) + as.numeric(m4_count))
) %>% 
group_by(selfing,DFE) %>% 
summarize(across(where(is.numeric), list(avg = mean))) 
write.csv(tidy_summary_table_final, file="/nas/longleaf/home/adaigle/DFESelfing/base_exp_averagep_m2.csv")

#above this is all I need for now. 

selfing_Ne <- function(selfing_prop) {
    5000/(1+(selfing_prop/(2-selfing_prop)))
}
theta <- 4*5000*3.3e-9*100
theta_50 <- 4*selfing_Ne(0.5)*3.3e-9*100
theta_99 <- 4*selfing_Ne(0.99)*3.3e-9*100
#tidy_summary_table <- tidy_summary_table %>% mutate(filtered_lines = NULL, newdata = NULL,newdata1 = NULL,newdata2 = NULL,newdata3 = NULL, m1 = NULL, m2 = NULL, m4 = NULL, neutral_m4 = NULL)

allele_freqs <- tidy_summary_table %>% #mutate(newdata3 = lapply(newdata2, function(x) x %>% filter(as.numeric(fix_gen) > 50000))) %>%
    mutate(#m1 = lapply(newdata3, function(df) sum(grepl("m1", df$mutation_class))), # quick way to count m1 rows 
    pis = unlist(lapply(newdata2, function(z) sum(2*(as.numeric(z$V9)/100) *(1-(as.numeric(z$V9)/100)))/total_neu_mutations)),
    allele_freqs = lapply(newdata2, function(z) as.numeric(z$V9)),
    mean_allele_freqs = unlist(lapply(allele_freqs, mean)),
    sd_allele_freqs = unlist(lapply(allele_freqs, sd))
    )

allele_freqs_summary <- allele_freqs[c(3:4,9,11:12)]
#write.csv(allele_freqs_summary, file = "/nas/longleaf/home/adaigle/DFESelfing/nolinkage_scripts/h05_allelefreqs.csv")

#quick_comparison_plot_noself <-  allele_freqs$allele_freqs[[1]]
#quick_comparison_plot_50self <-  allele_freqs$allele_freqs[[4]]
#quick_comparison_plot_99self <-  allele_freqs$allele_freqs[[7]]
#proportions1 <- table(quick_comparison_plot_noself) / length(quick_comparison_plot_noself)
#proportions2 <- table(quick_comparison_plot_50self) / length(quick_comparison_plot_50self)
#proportions3 <- table(quick_comparison_plot_99self) / length(quick_comparison_plot_99self)
#
## Create a data frame for the histograms
#columns <- c("Group", "Value", "Value.Freq")
#df1 <- data.frame(Group = factor(rep("s = 0%", length(proportions1))), Value = proportions1)
#df2 <- data.frame(Group = factor(rep("s = 50%", length(proportions2))), Value = proportions2)
#df3 <- data.frame(Group = factor(rep("s = 99%", length(proportions3))), Value = proportions3)
#
#colnames(df1) <- columns
#colnames(df2) <- columns
#colnames(df3) <- columns
#
#combined_df <- rbind(df1, df2,df3)
#
#ggplot(combined_df, aes(x = Value, y = Value.Freq, fill = Group)) +
#  geom_bar(stat = "identity", position = "dodge", color = "black") +
#  scale_fill_manual(values = c("s = 0%" = "blue", "s = 50%" = "red", "s = 99%" = "green")) +
#  labs(title = "DFE1, 50% Selfing", x = "Value", y = "Proportion") +
# theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
#  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
#  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))
  
##quick_comparison_ploth05 <-  allele_freqs$allele_freqs[[9]]
proportions1 <- table(quick_comparison_ploth025) / length(quick_comparison_ploth025)
proportions2 <- table(quick_comparison_ploth05) / length(quick_comparison_ploth05)
# Create a data frame for the histograms
columns <- c("Group", "Value", "Value.Freq")
df1 <- data.frame(Group = factor(rep("h = 0.25", length(proportions1))), Value = proportions1)
df2 <- data.frame(Group = factor(rep("h = 0.5", length(proportions2))), Value = proportions2)
colnames(df1) <- columns
colnames(df2) <- columns

combined_df <- rbind(df1, df2)

ggplot(combined_df, aes(x = Value, y = Value.Freq, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("h = 0.25" = "blue", "h = 0.5" = "red")) +
  labs(title = "DFE3, 99% Selfing", x = "Value", y = "Proportion") +
 theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))
  