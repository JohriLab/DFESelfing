# a script to calculate pi using m1 mutations
# in slim output files
library(tidyverse)

total_neu_mutations <- 187500 * 6
total_sel_mutations <- 562500 * 6
theta <- 4*5000*3.3e-9*100

#need to read in all m1_output.txt files in a way that preserves selfing%_DFE#_replicate structure
#get pi for each output
#average and std dev for each DFE and make table like I did for gamma and beta

#dirs <- dir(path = "/nas/longleaf/home/adaigle/work/johri_elegans/gammaDFE", pattern = "DFE")

dirs <- list.dirs(path = results_dir)
dirs <- dirs[ grepl("eqm", dirs) ]
dirs <- dirs[ !grepl("subpops", dirs) ]

# has code to load in files with nice names, but i realized I didn't need that 
#makes a vector of avg pi for each experiment 
#outputs are averaged and std deviation put into table
#summary_table <- data.frame()
#
#for(x in dirs) {
#    replicate_pis <- c()
#    replicate_seg_sites <- c()
#    for(y in list.files(path = x, pattern = "m1")) {
#        #assign( 
#        #paste(sub(
#        #    sub("[^/]+$", "", x),
#        #   "", x), y, sep=""), read.table(paste(x, "/", y, sep = "")))
#        z <- read.table(paste(x, "/", y, sep = ""))
#        replicate_pis <- c(replicate_pis, sum(2*(z$V9/100) *(1-(z$V9/100)))/total_neu_mutations)
#        replicate_seg_sites <- c(replicate_seg_sites,sum(z$V9/z$V9))
#    }
#    summary_table <- rbind(summary_table, c(sub(sub("[^/]+$", "", x),"", x),
#        mean(replicate_pis), sd(replicate_pis), mean(replicate_seg_sites), sd(replicate_seg_sites)))
#}
#
#row.names(summary_table) <- summary_table[[1]]
#summary_table <- summary_table[-1]
#colnames(summary_table) <- c("pi", "pi_sd", "segregating_sites", "segregating_sites_sd")
##relative biases due to selfing/dfe
##E(pi) is theta in a neutral model
##aka B, background selection
#summary_table$B <- as.numeric(summary_table$pi)/(theta)
#summary_table$Bsd <- as.numeric(summary_table$pi_sd)/(theta)
##empircal Ne = B*Ne
#summary_table$empirical_Ne <- summary_table$B*5000
##write.csv(summary_table, file="/nas/longleaf/home/adaigle/DFESelfing/pi_summary.csv")
#
#summary_table$Nempirical <- as.numeric(summary_table$pi)/(4*3.3e-9*100)
#
#
#pos_dirs <- list.dirs(path = "/nas/longleaf/home/adaigle/work/johri_elegans/gammaDFE/pos")
#pos_dirs <- pos_dirs[ grepl("eqm", pos_dirs) ]

# has code to load in files with nice names, but i realized I didn't need that 
#makes a vector of avg pi for each experiment 
#outputs are averaged and std deviation put into table
#pos_summary_table <- data.frame()
#for(x in pos_dirs) {
#    replicate_pis <- c()
#    for(y in list.files(path = x, pattern = "m1")) {
#        #assign( 
#        #paste(sub(
#        #    sub("[^/]+$", "", x),
#        #   "", x), y, sep=""), read.table(paste(x, "/", y, sep = "")))
#        z <- read.table(paste(x, "/", y, sep = ""))
#        replicate_pis <- c(replicate_pis, sum(2*(z$V9/100) *(1-(z$V9/100)))/total_neu_mutations)
#        #avg <- mean z$V10
#    }
#    pos_summary_table <- rbind(pos_summary_table, c(sub(sub("[^/]+$", "", x),"", x),
#        mean(replicate_pis), sd(replicate_pis)))
#}
#
#row.names(pos_summary_table) <- pos_summary_table[[1]]
#pos_summary_table <- pos_summary_table[-1]
#colnames(pos_summary_table) <- c("pi", "pi_sd")
#relative biases due to selfing/dfe
#E(pi) is theta in a neutral model
#aka B, background selection
#pos_summary_table$B <- as.numeric(pos_summary_table$pi)/(theta)
#empircal Ne = B*Ne
#pos_summary_table$empirical_Ne <- pos_summary_table$B*5000



#commenting this out because it takes forever to run...

#seg_site_table <- data.frame()
#for(x in dirs) {
#    replicate_m1_seg_sites <- c()
#    replicate_m2_seg_sites <- c()
#    replicate_m3_seg_sites <- c()
#    replicate_m4_seg_sites <- c()
#    m2_pi <- c()
#    for(y in list.files(path = x, pattern = "^output\\d.txt")) {
#        # Read the file as a character vector
#        lines <- readLines(paste(x, "/", y, sep = ""))
#
#        # Select lines containing "m1", "m2", "m3", or "m4"
#        selected_lines <- lines[grep("m1|m2|m3|m4", lines)]
#        # Remove the first three lines
#        #selected_lines <- selected_lines[-(1:3)]
#        # Read the selected lines into a data frame
#        mydata <- read.table(text = selected_lines)
#        #z <- read.table(paste(x, "/", y, sep = ""))
#        #replicate_pis <- c(replicate_pis, sum(2*(z$V9/100) *(1-(z$V9/100)))/total_neu_mutations)
#        #replicate_seg_sites <- c(replicate_seg_sites,sum(z$V9/z$V9))
#        m1_rows <- subset(mydata, grepl("m1", mydata[[3]]))
#        replicate_m1_seg_sites <- c(replicate_m1_seg_sites,length(m1_rows$V9))
#
#        m2_rows <- subset(mydata, grepl("m2", mydata[[3]]))
#        replicate_m2_seg_sites <- c(replicate_m2_seg_sites,length(m2_rows$V9))
#        m2_pi <- c(m2_pi, sum(2*(m2_rows$V9/100) *(1-(m2_rows$V9/100)))/total_sel_mutations)
#
#        m3_rows <- subset(mydata, grepl("m3", mydata[[3]]))
#        replicate_m3_seg_sites <- c(replicate_m3_seg_sites,length(m3_rows$V9))
#        m4_rows <- subset(mydata, grepl("m4", mydata[[3]]))
#        replicate_m4_seg_sites <- c(replicate_m4_seg_sites,length(m4_rows$V9))    
#    }
#    seg_site_table <- rbind(seg_site_table, c(sub(sub("[^/]+$", "", x),"", x),
#        mean(replicate_m1_seg_sites), sd(replicate_m1_seg_sites), 
#        mean(replicate_m2_seg_sites), sd(replicate_m2_seg_sites),
#        mean(m2_pi), sd(m2_pi),
#        mean(replicate_m3_seg_sites), sd(replicate_m3_seg_sites), 
#        mean(replicate_m4_seg_sites), sd(replicate_m4_seg_sites)))
#}
#row.names(seg_site_table) <- seg_site_table[[1]]
#seg_site_table <- seg_site_table[-1]
#colnames(seg_site_table) <- c("m1_mean", "m1_sd", 
#    "m2_mean", "m2_sd","m2_pi_mean", "m2_sd_mean", "m3_mean", "m3_sd","m4_mean", "m4_sd")
#pi_segsites <- merge(summary_table, seg_site_table, by = 0)
#pi_segsites2 <- data.frame(pi_segsites, DFE = sapply(strsplit(pi_segsites$Row.names, "_"), `[`, 3),
#                 selfing = paste0(gsub("^eqm_selfing([0-9]+).*$", "\\1", pi_segsites$Row.names), "%"),
#                 pos = grepl("_pos$", pi_segsites$Row.names))
#pi_segsites2 <- pi_segsites2[order(pi_segsites2$selfing), ]
#pi_segsites2 <- pi_segsites2[order(pi_segsites2$DFE), ]
#pi_segsites2 <- pi_segsites2[order(pi_segsites2$pos), ]
#
#write.csv(pi_segsites2, file="/nas/longleaf/home/adaigle/DFESelfing/pi_segsites_summary.csv")

####new tidy way of doing this:

#tidy_summary_table <- tibble(
#    name = list.files(path = dirs, pattern = "m1_output..txt"),
#    fullpath = paste0(dirs, "/", name),
#    selfing = str_extract(fullpath, "(?<=eqm_selfing)\\d+"),
#    DFE = str_extract(fullpath, "(DFE)\\d+"),
#    output = str_extract(fullpath, "(?<=output)\\d+"),
#    data = lapply(fullpath, read.table),
#    pis = unlist(lapply(data, function(z) sum(2*(z$V9/100) *(1-(z$V9/100)))/total_neu_mutations)),
#    B = pis/theta,
#    empirical_Ne = B*5000
#) %>% 
#mutate(data = NULL)#don't need the big dataframes anymore

string_pattern <- c("m1")
#missing two replicates from 0% selfing
#due to excessive sim run time and high similarity of results these weren't rerun
dirs <- dirs[!grepl("_0/DFE1output2.txt", dirs) & !grepl("_0/DFE1output3.txt", dirs)]
dirs <- dirs[!grepl("_0/DFE2output2.txt", dirs) & !grepl("_0/DFE2output3.txt", dirs)]
dirs <- dirs[!grepl("_0/DFE3output2.txt", dirs) & !grepl("_0/DFE3output3.txt", dirs)]

tidy_summary_table <- tibble(
    name = list.files(path = dirs, pattern = "^output..txt"),
    fullpath = paste0(dirs, "/", name)) %>% 
  filter(!grepl("g0_DFE1/output2|g0_DFE2/output2|g0_DFE3/output2", fullpath))

tidy_summary_table <- tidy_summary_table %>% mutate(
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

tidy_summary_table <- tidy_summary_table %>% #mutate(newdata3 = lapply(newdata2, function(x) x %>% filter(as.numeric(fix_gen) > 50000))) %>%
    mutate(#m1 = lapply(newdata3, function(df) sum(grepl("m1", df$mutation_class))), # quick way to count m1 rows 
    pis = unlist(lapply(newdata2, function(z) sum(2*(as.numeric(z$V9)/100) *(1-(as.numeric(z$V9)/100)))/total_neu_mutations)),
    B = pis/theta,
    empirical_Ne = B*5000,
    filtered_lines = NULL,
    newdata = NULL,
    newdata2 = NULL
    #m2 = lapply(newdata3, function(df) sum(grepl("m2", df$mutation_class))),
    #m4 = lapply(newdata3, function(df) sum(grepl("m4", df$mutation_class))),
    #m1_count = unlist(lapply(m1, sum)), # quick way to count m1 rows  pt 2
    #m2_count = unlist(lapply(m2, sum)),
    #m4_count = unlist(lapply(m4, sum)),
    #neutral_m4 = lapply(newdata3, function(x) x %>% filter((mutation_class == "m4"), 2*5000*as.numeric(s) < 5 )),
    #neutral_m4_count = unlist(lapply(neutral_m4, function(x) nrow(x))),
    #alpha = (as.numeric(m4_count) - neutral_m4_count)/ (as.numeric(m2_count) + as.numeric(m4_count))
) 
#%>% 
#group_by(selfing,DFE) %>% 
#summarize(across(where(is.numeric), list(avg = mean))) 

B_value_table <-tidy_summary_table %>% 
group_by(selfing,DFE) %>% 
summarize(across(where(is.numeric), list(avg = mean))) 