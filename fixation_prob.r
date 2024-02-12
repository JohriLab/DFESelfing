library(tidyverse)
library(reshape2)
library(ggpubr)
library(scales)
base_dir <- "/nas/longleaf/home/adaigle/DFESelfing/"

source(paste0(base_dir, "scripts/calculate_pi.r"))
sim_outputs_dir <- "/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/"
tidy_summary_table_pi <- tidy_summary_table
figures_dir <- paste0(base_dir, "figures_for_publication/")
dirs <- list.dirs(path = paste0(sim_outputs_dir, "/original_simulations"))
dirs <- dirs[ grepl("eqm", dirs) ]
#remove sims with beneficial muts
#dirs <- dirs[ !grepl("pos", dirs) ]
string_pattern <- c("m2", "m1")

#2 x N x mu x no. selected sites x no. generations
total_neu_sites <- 187500
total_sel_sites <- 562500
mu_scaled <- 3.3e-9*100
N_scaled <- 5000
generations <- 20000
theta <- 4*5000*3.3e-9*100

#kind of expected fixations except I'm not accounting for s yet
expected_fixations <- 2 * N_scaled * mu_scaled * total_sel_sites * 20000



tidy_summary_table <- tibble(
    name = list.files(path = dirs, pattern = "^output..fixed"),
    fullpath = paste0(dirs, "/", name),
    selfing = str_extract(fullpath, "(?<=eqm_selfing)\\d+"),
    DFE = str_extract(fullpath, "(DFE)\\d+"),
    output = str_extract(fullpath, "(?<=output)\\d+"),
    data = lapply(fullpath, readLines), # must do this as rows have different column numbers
    filtered_lines = lapply(data, function(x) x[grep(paste(string_pattern, collapse="|"), x)])
    
) %>% 
mutate(data = NULL, m1 = NULL, m2 = NULL, m4 = NULL)


# have to do some gymnastics now to turn my lists to df's, filter out muts fixed before 
tidy_summary_table <- tidy_summary_table %>% 
    mutate(newdata = lapply(filtered_lines, as.data.frame))

tidy_summary_table <- tidy_summary_table %>% 
    mutate(newdata2= lapply(newdata, function(x) separate(x, 1, into = c("id1", "id2", "mutation_class", "id3", "s", "h", "p1", "start_gen", "fix_gen"), sep = " ")
    ))

#have to change table name here bc I have an identically named table I'm pulling from another script
tidy_summary_table_fixed <- tidy_summary_table %>% mutate(newdata3 = lapply(newdata2, function(x) x %>% filter(as.numeric(fix_gen) > 50000))) %>%
    mutate(m1 = lapply(newdata3, function(df) sum(grepl("m1", df$mutation_class))), # quick way to count m1 rows 
    m2 = lapply(newdata3, function(df) sum(grepl("m2", df$mutation_class))),
    m1_count = unlist(lapply(m1, sum)), # quick way to count m1 rows  pt 2
    m2_count = unlist(lapply(m2, sum)),
    m2_filtered = lapply(newdata3, function(df) df[grepl("m2", df$mutation_class), ]),
    m1_filtered = lapply(newdata3, function(df) df[grepl("m1", df$mutation_class), ]),
) %>%
  mutate(beta = ifelse(DFE == "DFE1", 0.9,
                       ifelse(DFE == "DFE2", 0.5,
                              ifelse(DFE == "DFE3", 0.3, NA)))) %>%
  mutate(gamma = ifelse(DFE == "DFE1", 5,
                       ifelse(DFE == "DFE2", 50,
                              ifelse(DFE == "DFE3", 1000, NA))))

#source("/nas/longleaf/home/adaigle/DFESelfing/calculate_pi.r")

selfing_Ne <- function (selfing) {
  pself <- selfing * .01
  f <- pself / (2 - pself)
  nself <- 5000 / (1 + f)
  return(nself)
}

tidy_summary_table <- tidy_summary_table_pi[ !grepl("pos", tidy_summary_table$fullpath), ] %>% mutate(selfing_Ne = selfing_Ne(as.numeric(selfing)), 
  noselfing_B = empirical_Ne / selfing_Ne)

gamma_prop_calc <- function(B, meanGamma,beta,s1_count, s2_count, s3_count, s4_count, s5_count, s6_count, s7_count, s8_count, s9_count, s10_count) { 
    # calc fixation prob for classes of muts in given range for given dfe
    Nw <- 100 
    meanS <- B * (meanGamma/(2*Nw))

    # this code adapted from get_dfe_class_proportions
    s_shape <- beta
    s_rate <- s_shape/abs(meanS)
    x1 <- 1/(2*Nw)
    x2 <- 2/(2*Nw)
    x3 <- 3/(2*Nw)
    x4 <- 4/(2*Nw)
    x5 <- 5/(2*Nw)
    x6 <- 6/(2*Nw)
    x7 <- 7/(2*Nw)
    x8 <- 8/(2*Nw)
    x9 <- 9/(2*Nw)
    x10 <- 10/(2*Nw)
    f0 <- pgamma(x1, shape=s_shape, rate=s_rate)
    f1 <- pgamma(x2, shape=s_shape, rate=s_rate) - f0
    f2 <- pgamma(x3, shape=s_shape, rate=s_rate) - f1 - f0
    f3 <- pgamma(x4, shape=s_shape, rate=s_rate) -f2 - f1 - f0
    f4 <- pgamma(x5, shape=s_shape, rate=s_rate) -f2 - f1 - f0 - f3
    f5 <- pgamma(x6, shape=s_shape, rate=s_rate) -f2 - f1 - f0 - f3 - f4
    f6 <- pgamma(x7, shape=s_shape, rate=s_rate)-f2 - f1 - f0 - f3 - f4 - f5
    f7 <- pgamma(x8, shape=s_shape, rate=s_rate) - f1 - f0 - f3 - f4 - f5 - f6
    f8 <- pgamma(x9, shape=s_shape, rate=s_rate) - f1 - f0 - f3 - f4 - f5 - f6 - f7
    f9 <- pgamma(x10, shape=s_shape, rate=s_rate) - f1 - f0 - f3 - f4 - f5 - f6 - f7 - f8
    p0 <- s1_count/(total_sel_sites * 2 * 5000 * mu_scaled * generations * f0)
    p1 <- s2_count/(total_sel_sites * 2 * 5000 * mu_scaled * generations * f1)
    p2 <- s3_count/(total_sel_sites * 2 * 5000 * mu_scaled * generations * f2)
    p3 <- s4_count/(total_sel_sites * 2 * 5000 * mu_scaled * generations * f3 )
    p4 <- s5_count/(total_sel_sites * 2 * 5000 * mu_scaled * generations * f4 )
    p5 <- s6_count/(total_sel_sites * 2 * 5000 * mu_scaled * generations * f5 )
    p6 <- s7_count/(total_sel_sites * 2 * 5000 * mu_scaled * generations * f6)
    p7 <- s8_count/(total_sel_sites * 2 * 5000 * mu_scaled * generations * f7 )
    p8 <- s9_count/(total_sel_sites * 2 * 5000 * mu_scaled * generations * f8 )
    p9 <- s10_count/(2 * 5000 * mu_scaled * generations * f9)
    #print(paste(f0, f1, f2, f3, sep=", "))
    #return(c(list(f0), list(f1), list(f2), list(f3)))
    return(c(p0 = p0, p1 = p1, p2 = p2, p3 = p3, p4=p4, p5=p5, p6=p6, p7=p7, p8=p8, p9=p9))
}



tidy_summary_table_fixed <- tidy_summary_table_fixed %>% 
  mutate(s1 = map(m2_filtered, ~ filter(.x, -(as.numeric(s)) <= 1/(5000*2))),
    s1_count = lapply(s1, function(s1) length(s1$s)),
    s1_count = unlist(lapply(s1_count, sum)),
    s1_prop = s1_count/m2_count,
    s2 = map(m2_filtered, ~ filter(.x, -(as.numeric(s)) >= 1/(5000*2) & -(as.numeric(s)) <= 2/(5000*2))),
    s2_count = lapply(s2, function(s2) length(s2$s)),
    s2_count = unlist(lapply(s2_count, sum)),
    s2_prop = s2_count/m2_count,
    s3 = map(m2_filtered, ~ filter(.x, -(as.numeric(s)) >= 2/(5000*2) & -(as.numeric(s)) <= 3/(5000*2))),
    s3_count = lapply(s3, function(s3) length(s3$s)),
    s3_count = unlist(lapply(s3_count, sum)),
    s3_prop = s3_count/m2_count,
    s4 = map(m2_filtered, ~ filter(.x, -(as.numeric(s)) >= 3/(5000*2) & -(as.numeric(s)) <= 4/(5000*2))),
    s4_count = lapply(s4, function(s4) length(s4$s)),
    s4_count = unlist(lapply(s4_count, sum)),
    s4_prop = s4_count/m2_count,
    s5 = map(m2_filtered, ~ filter(.x, -(as.numeric(s)) >= 4/(5000*2) & -(as.numeric(s)) <= 5/(5000*2))),
    s5_count = lapply(s5, function(s5) length(s5$s)),
    s5_count = unlist(lapply(s5_count, sum)),
    s5_prop = s5_count/m2_count,
    s6 = map(m2_filtered, ~ filter(.x, -(as.numeric(s)) >= 5/(5000*2) & -(as.numeric(s)) <= 6/(5000*2))),
    s6_count = lapply(s6, function(s6) length(s6$s)),
    s6_count = unlist(lapply(s6_count, sum)),
    s6_prop = s6_count/m2_count,
    s7 = map(m2_filtered, ~ filter(.x, -(as.numeric(s)) >= 6/(5000*2) & -(as.numeric(s)) <= 7/(5000*2))),
    s7_count = lapply(s7, function(s7) length(s7$s)),
    s7_count = unlist(lapply(s7_count, sum)),
    s7_prop = s7_count/m2_count,
    s8 = map(m2_filtered, ~ filter(.x, -(as.numeric(s)) >= 7/(5000*2) & -(as.numeric(s)) <= 8/(5000*2))),
    s8_count = lapply(s8, function(s8) length(s8$s)),
    s8_count = unlist(lapply(s8_count, sum)),
    s8_prop = s8_count/m2_count,
    s9 = map(m2_filtered, ~ filter(.x, -(as.numeric(s)) >= 8/(5000*2) & -(as.numeric(s)) <= 9/(5000*2))),
    s9_count = lapply(s9, function(s9) length(s9$s)),
    s9_count = unlist(lapply(s9_count, sum)),
    s9_prop = s9_count/m2_count,
    s10 = map(m2_filtered, ~ filter(.x, -(as.numeric(s)) >= 9/(5000*2) & -(as.numeric(s)) <= 10/(5000*2))),
    s10_count = lapply(s10, function(s10) length(s10$s)),
    s10_count = unlist(lapply(s10_count, sum)),
    s10_prop = s10_count/m2_count,
    ) 


#binding pi/B table to fixation table 
tidy_summary_table <- tidy_summary_table %>% select(3:4, 10) %>%
  group_by(selfing,DFE) %>% 
  summarize(across(where(is.numeric), list(avg = mean, sd = sd))) 

tidy_summary_table_fixed <- merge(tidy_summary_table_fixed, tidy_summary_table, by = c("DFE", "selfing"), all.x = TRUE)

tidy_summary_table_fixed <- tidy_summary_table_fixed %>% 
  mutate(output = pmap(list(1, gamma, beta, s1_count, s2_count, s3_count, s4_count, s5_count, s6_count, s7_count, s8_count, s9_count,s10_count), gamma_prop_calc)) %>%
  unnest_wider(output) 


quick_test <- tidy_summary_table_fixed %>% 
  mutate(across(5:14, ~ ./ (1/10000))) %>% 
  select(1:14,57) %>%
  group_by(selfing,DFE) %>% 
  summarize(across(where(is.numeric), list(avg = mean, sd = sd))) 
  #cbind(tidy_summary_table)
  #%>%
  #melt() %>% 
  #mutate(value = ifelse(is.na(value), 0, value)) %>% 
  #filter(variable == "value" | paste0(generation, "_sd") == variable) 
  #pivot_wider(tidy_summary_table, 
  #id_cols = c("DFE","selfing"), names_from = "variable", values_from = "value") %>%



ggplot(quick_test, aes(x = selfing, y = p2_avg, fill = selfing)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 1) +
  #labs(title = "DFEalpha", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = p2_avg - p2_sd, ymax = p2_avg + p2_sd), position = position_dodge(width = 0.9)) +
  #expand_limits(y=c(0,1)) + 
  #scale_fill_manual(values = c("#404040", hue_pal()(6))) +
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), 
  plot.title= element_text(size=25), legend.title = element_text(size=15), legend.text = element_text(size=15))



DFE_proportions_truth <- function(B, meanGamma,beta) { 
    # modified function to calc vector of dfe classes given gamma and beta
    Nw <- 100 
    meanS <- B * (meanGamma/(2*Nw))

    # this code adapted from get_dfe_class_proportions
    s_shape <- beta
    s_rate <- s_shape/abs(meanS)
    x1 <- 1/(2*Nw)
    x10 <- 10/(2*Nw)
    x100 <- 100/(2*Nw)
    f0 <- pgamma(x1, shape=s_shape, rate=s_rate)
    f1 <- pgamma(x10, shape=s_shape, rate=s_rate) - f0
    f2 <- pgamma(x100, shape=s_shape, rate=s_rate) - f1 - f0
    f3 <- 1.0 - f2 - f1 - f0
    #print(paste(f0, f1, f2, f3, sep=", "))
    #return(c(list(f0), list(f1), list(f2), list(f3)))
    return(c(f0 = f0, f1 = f1, f2 = f2, f3 = f3))
}



dataframe_of_truth <- tidy_summary_table_fixed  %>%
mutate(output = pmap(list(1, gamma, beta), DFE_proportions_truth)) %>%
    unnest_wider(output) 



# Example usage
meanGamma <- 1000
beta <- 0.3
Ns <- 100


### Attempt 3

# Mutation class density function
density_func <- function(x) {
  meanGamma <- 1000
  shape <- beta  # Shape parameter
  N <- 5000
  meanS <- (meanGamma)
  rate <- shape/abs(meanS)
  dgamma(x, shape = beta, rate = rate)
  
}

# Parameters for Kimura's equation
N <- 5000  # Effective population size
t <- 20000   # Time in generations

# Mutation class ranges
#class_ranges <- list(
#  class1 = c(0, 1),     # Range for class 1 mutations
#  class2 = c(1, 2),   # Range for class 2 mutations
#  class3 = c(2, 3)    # Range for class 3 mutations
#)
#
## Calculate the expected proportion of fixations in each mutation class
#proportions <- numeric(length(class_ranges))
#num_fixations <- numeric(length(class_ranges))
#
#for (i in seq_along(class_ranges)) {
#  range <- class_ranges[[i]]
#  lower <- range[1]
#  upper <- range[2]
#  
#  # Integrate the mutation class density function over the range of the class
#  integral <- integrate(density_func, lower = lower, upper = upper)$value
#  
#  # Apply Kimura's equation to calculate the expected proportion of fixations
#  proportions[i] <- integral * (1 - exp(-2 * N * integral * t))
#
#  # Calculate the expected number of fixations
#  num_fixations[i] <- proportions[i] * mu_scaled * total_sel_sites * t
#
#}
#
## Normalize the proportions to sum up to 1
#proportions <- proportions / sum(proportions)
#
## Print the expected proportions of fixations in each mutation class
#for (i in seq_along(class_ranges)) {
#  cat("Class", i, "Proportion:", proportions[i], "Number of Fixations:", num_fixations[i], "\n")
#}


#then I need to calc the real Fixation prob using prop sub in data

Pfix <- function(s, B) {
  (s*B) / (exp(2 * N * B * s) - 1)
}
meanGamma <- 5
shape <- beta  # Shape parameter
N <- 5000
meanS <- (meanGamma/(2*N))
rate <- shape/abs(meanS)
s0 <- 0 / (2 * N)
s1 <- 1 / (2 * N)
s2 <- 2 / (2 * N)
s3 <- 3 / (2 * N)
s4 <- 4 / (2 * N)
s5 <- 5 / (2 * N)
s6 <- 6 / (2 * N)
s7 <- 7 / (2 * N)
s8 <- 8 / (2 * N)
s9 <- 9 / (2 * N)
s10 <- 10/(2*N)

fixprob <- function (s1,s2,B) { 
# Integrate over the specified range
  Integral <- integrate(Pfix, B=B, lower = s1, upper = s2)
  MeanProb <- Integral$value * (1 / (s2 - s1))
  relMeanProb <- MeanProb / (1 / (2 * N))
  return(relMeanProb)
}

filter <- quick_test  %>%
mutate(p0_truth = unlist(pmap(list(s0,s1, noselfing_B_avg_avg), fixprob))) %>%
mutate(p1_truth = unlist(pmap(list(s1,s2, noselfing_B_avg_avg), fixprob))) %>%
mutate(p2_truth = unlist(pmap(list(s2,s3, noselfing_B_avg_avg), fixprob))) %>%
mutate(p3_truth = unlist(pmap(list(s3,s4, noselfing_B_avg_avg), fixprob)))

adjusted_truths <- filter %>% select(1:2, 25:28) %>%
  mutate(selfing = case_when(
        selfing == '0'  ~ 'true0', 
        selfing == '50' ~ 'true50', 
        selfing == '80' ~ 'true80', 
        selfing == '90' ~ 'true90',
        selfing == '95' ~ 'true95',
        selfing == '99' ~ 'true99',
        TRUE ~ selfing
    ))
colnames(adjusted_truths) <- c("selfing", "DFE", "p0_avg", "p1_avg", "p2_avg", "p3_avg")

truth1 <- data.frame("truth", "DFE1", fixprob(s0,s1,1) , 0, fixprob(s1,s2,1) , 0, fixprob(s2,s3,1) , 
  0, fixprob(s3,s4,1) , 0, fixprob(s4,s5,1) , 0)
truth2 <- data.frame("truth", "DFE2", fixprob(s0,s1,1) , 0, fixprob(s1,s2,1) , 0, fixprob(s2,s3,1) , 
  0, fixprob(s3,s4,1) , 0, fixprob(s4,s5,1) , 0)
truth3 <- data.frame("truth", "DFE3", fixprob(s0,s1,1) , 0, fixprob(s1,s2,1) , 0, fixprob(s2,s3,1) , 
  0, fixprob(s3,s4,1) , 0, fixprob(s4,s5,1) , 0)
colnames(truth1) <- colnames(filter[1:12])
colnames(truth2) <- colnames(filter[1:12])
colnames(truth3) <- colnames(filter[1:12])

plot_fixprob <- filter %>% select(1:12) %>% rbind(truth1, truth2, truth3, adjusted_truths, filter)
#selfing_order <- c("truth", 0, 50, 80, 90, 95, 99)
selfing_order <- c("truth", "true0", "0", "true50", "50", "true80", "80", 
  "true90", "90","true95", "95", "true99", "99")
p0 <- ggplot(plot_fixprob, aes(x = factor(selfing, levels = selfing_order), y = p0_avg, fill = factor(selfing, levels = selfing_order))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 1) +
  labs(title = "Ns=0-1 Fixation Probabilities", x = "Selfing Level", y = "fixation probability", fill = "Selfing %") +
  geom_errorbar(aes(ymin = p0_avg - p0_sd, ymax = p0_avg + p0_sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) + 
  scale_fill_manual(values = c("#404040", hue_pal()(12))) +
  theme(axis.text.x=element_text(size=0), axis.text.y=element_text(size=15), 
    axis.title.x=element_text(size=15),axis.title.y=element_text(size=15), strip.text = element_text(size=15),
    plot.title= element_text(size=0), legend.position = "none", legend.text = element_text(size=12))

p1 <- ggplot(plot_fixprob, aes(x = factor(selfing, levels = selfing_order), y = p1_avg, fill = factor(selfing, levels = selfing_order))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 1) +
  labs(title = "Ns=1-2 Fixation Probabilities", x = "Selfing Level", y = "fixation probability", fill = "Selfing %") +
  geom_errorbar(aes(ymin = p1_avg - p1_sd, ymax = p1_avg + p1_sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) + 
  scale_fill_manual(values = c("#404040", hue_pal()(12))) +
  theme(axis.text.x=element_text(size=0), axis.text.y=element_text(size=15), 
    axis.title.x=element_text(size=15),axis.title.y=element_text(size=15), strip.text = element_text(size=15),
    plot.title= element_text(size=0), legend.position = "none", legend.text = element_text(size=12))

p2 <- ggplot(plot_fixprob, aes(x = factor(selfing, levels = selfing_order), y = p2_avg, fill = factor(selfing, levels = selfing_order))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 1) +
  labs(title = "Ns=2-3 Fixation Probabilities", x = "Selfing Level", y = "fixation probability", fill = "Selfing %") +
  geom_errorbar(aes(ymin = p2_avg - p2_sd, ymax = p2_avg + p2_sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) + 
  scale_fill_manual(values = c("#404040", hue_pal()(12))) +
  theme(axis.text.x=element_text(size=0), axis.text.y=element_text(size=15), 
    axis.title.x=element_text(size=15),axis.title.y=element_text(size=15), strip.text = element_text(size=15),
    plot.title= element_text(size=0), legend.position = "none", legend.text = element_text(size=12))

p3 <- ggplot(plot_fixprob, aes(x = factor(selfing, levels = selfing_order), y = p3_avg, fill = factor(selfing, levels = selfing_order))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 1) +
  labs(title = "Ns=3-4 Fixation Probabilities", x = "Selfing Level", y = "fixation probability", fill = "Selfing %") +
  geom_errorbar(aes(ymin = p3_avg - p3_sd, ymax = p3_avg + p3_sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) + 
  scale_fill_manual(values = c("#404040", hue_pal()(12))) +
  theme(axis.text.x=element_text(size=0), axis.text.y=element_text(size=15), 
    axis.title.x=element_text(size=15),axis.title.y=element_text(size=15), strip.text = element_text(size=15),
    plot.title= element_text(size=0), legend.position = "none", legend.text = element_text(size=12))

ggarrange(p0, p1, p2, p3,
                    labels = c("A", "B", "C", "D"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 2,
                    common.legend = TRUE, legend = "none", vjust=1)
