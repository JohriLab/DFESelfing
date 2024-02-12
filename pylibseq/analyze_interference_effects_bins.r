library(tidyverse)
library(ggplot2)
library(ggpubr)

#df <- read.table("/nas/longleaf/home/adaigle/DFESelfing/pylibseq/allresults_window_5000_1,2,3,4,5BinSize50.stats", header=F, skip=1)
df <- read.table("/nas/longleaf/home/adaigle/DFESelfing/pylibseq/allresults_window_5000_1,2,3,4,5BinSize500.stats", header=F, skip=1)

gammabeta <- read.csv("/nas/longleaf/home/adaigle/DFESelfing/pylibseq/gammabeta.csv")
length_dist_cols <- length(names(df))-15
header <- c("filename","output","site","posn","S","thetapi","thetaw","thetah","hprime","tajimasd","numSing","hapdiv",1:length_dist_cols,"rsq","D","Dprime")
#header <- c("filename","output","site","posn","S","thetapi","thetaw","thetah","hprime","tajimasd","numSing","hapdiv","rsq","D","Dprime")

colnames(df) <- header
# Get the column names to process (excluding the 'ID' column)
columns_to_process <- colnames(df)[-1]

# Iterate over each column and split it into numeric columns
for (col_name in columns_to_process) {
  if ("," %in% unlist(strsplit(as.character(df[1, col_name]), split=""))) {
    df <- df %>%
      separate(col_name, into = paste0(col_name, 1:4), sep = ",", convert = TRUE)
  }
}
#unique(df[c(1:3,270:273)])
df <- df %>% mutate(
    selfing = as.numeric(str_extract(filename, "(?<=eqm_selfing)\\d+")),
    DFE = str_extract(filename, "(DFE)\\d+")
    ) %>% group_by(selfing, DFE)
  
gammabeta <- gammabeta %>% select(B, empirical_Ne, DFE, selfing, gamma, b, GammaZero.negGmean, GammaZero.negGshape, selfing_Ne, newNE) %>%
  group_by(selfing, DFE)

join_df <- left_join(df, gammabeta, multiple="all") 

#summarize_df <- join_df %>% ungroup() %>% group_by(filename,site,DFE,selfing) %>% 
#    summarize(across(where(is.numeric), list(avg = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE)))) %>%
#    ungroup() %>% group_by(DFE)

#first find the average within each replicate, then find the sd using our five replicates
summarize_df <- join_df %>% ungroup() %>% group_by(filename,output,site,DFE,selfing) %>% 
    summarize(across(where(is.numeric), list(avg = ~mean(., na.rm = TRUE)), .names = "{col}")) %>% #.names makes sure columns aren't renamed to _avg yet
    ungroup() %>% group_by(DFE, selfing,site) %>%
    summarize(across(where(is.numeric), list(avg = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE))))

#plot(summarize_df$selfing, summarize_df$Dprime_avg)

DFE1 <- summarize_df %>% filter(DFE=="DFE1") %>% ungroup()

ggplot(DFE1, aes(x = selfing, y = Dprime_avg, fill=site)) +
  geom_errorbar(aes(ymin = Dprime_avg - Dprime_sd, ymax = Dprime_avg + Dprime_sd,color=site), 
    width = 5, position=position_dodge(0.1)) +
  geom_line(aes(color=site)) +
  geom_point() + 
  #facet_wrap(~ filename, ncol = 1) +  # Separate plots by 'filename'
  labs(title = "Dprime_avg with Error Bars",
       x = "Selfing",
       y = "Dprime_avg")

DFE2 <- summarize_df %>% filter(DFE=="DFE2") %>% ungroup()

ggplot(DFE2, aes(x = selfing, y = Dprime_avg, fill=site)) +
  geom_errorbar(aes(ymin = Dprime_avg - Dprime_sd, ymax = Dprime_avg + Dprime_sd,color=site), 
    width = 5, position=position_dodge(0.1)) +
  geom_line(aes(color=site)) +
  geom_point() + 
  #facet_wrap(~ filename, ncol = 1) +  # Separate plots by 'filename'
  labs(title = "Dprime_avg with Error Bars",
       x = "Selfing",
       y = "Dprime_avg")


DFE3 <- summarize_df %>% filter(DFE=="DFE3") %>% ungroup()

ggplot(DFE3, aes(x = selfing, y = Dprime_avg, fill=site)) +
  geom_errorbar(aes(ymin = Dprime_avg - Dprime_sd, ymax = Dprime_avg + Dprime_sd,color=site), 
    width = 5, position=position_dodge(0.1)) +
  geom_line(aes(color=site)) +
  geom_point() + 
  #facet_wrap(~ filename, ncol = 1) +  # Separate plots by 'filename'
  labs(title = "Dprime_avg with Error Bars",
       x = "Selfing",
       y = "Dprime_avg")

quick_summary <- DFE3 %>% select(selfing,site,Dprime_avg)
#write.csv(quick_summary, file = "/nas/longleaf/home/adaigle/DFESelfing/pylibseq/DFE3.csv", quote=F)

plot_stat <- function(DFE, stat) {
  ggplot(DFE, aes(x = selfing, y = get(paste0(stat, "_avg")), fill=site)) +
  geom_errorbar(aes(ymin = get(paste0(stat, "_avg")) - get(paste0(stat, "_sd")), ymax = get(paste0(stat, "_avg")) + get(paste0(stat, "_sd")),color=site), 
    width = 5, position=position_dodge(0.1)) +
  geom_line(aes(color=site)) +
  geom_point() + 
  #facet_wrap(~ filename, ncol = 1) +  # Separate plots by 'filename'
  labs(#title = "Dprime_avg with Error Bars",
       x = "Selfing",
       y = paste0(stat, "_avg"))
}

plot_DFEinference <- function(DFE) {
  ggplot(DFE, aes(x = selfing)) +

  geom_line(aes(y = B_avg), colour = "#F8766D") +
  geom_point(aes(y = B_avg)) + 
  geom_errorbar(aes(ymin = gamma_avg - gamma_sd, ymax = gamma_avg + gamma_sd,color="#F8766D"), 
    width = 5, position=position_dodge(0.1)) 
  #geom_line(aes(y = GammaZero.negGshape_avg), colour = "purple") + 
  #geom_point(aes(y = GammaZero.negGshape_avg)) 
  #geom_errorbar(aes(ymin = GammaZero.negGmean_avg - GammaZero.negGmean_sd, ymax = GammaZero.negGmean_avg + GammaZero.negGmean_sd,color="purple"), 
  #  width = 5, position=position_dodge(0.1))
  #scale_fill_manual(values=c("#F8766D", "purple"))
  #facet_wrap(~ filename, ncol = 1) +  # Separate plots by 'filename'
  #labs(#title = "Dprime_avg with Error Bars",
  #     x = "Selfing",
  #     y = paste0(stat, "_avg"))
}

DFE3_Dprime <- plot_stat(DFE3, "Dprime")
DFE3_D <- plot_stat(DFE3, "D")
DFE3_tajimasd <- plot_stat(DFE3, "tajimasd")
DFE3_hapdiv <- plot_stat(DFE3, "hapdiv")
DFE3_inference <- plot_DFEinference(DFE3)
plot_stat(DFE3, "14")


ggarrange(DFE3_inference, DFE3_Dprime, DFE3_tajimasd, DFE3_hapdiv,
                    labels = c("A", "B", "C", "D"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 2,
                    common.legend = TRUE, legend = "right")

 plot_stat(DFE1, "hprime")
DFE1_Dprime <- plot_stat(DFE1, "Dprime")
DFE1_D <- plot_stat(DFE1, "D")
DFE1_tajimasd <- plot_stat(DFE1, "tajimasd")
DFE1_hapdiv <- plot_stat(DFE1, "hapdiv")
DFE1_inference <- plot_DFEinference(DFE1)


ggarrange(DFE1_inference, DFE1_Dprime, DFE1_tajimasd, DFE1_hapdiv,
                    labels = c("A", "B", "C", "D"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 2,
                    common.legend = TRUE, legend = "right")

DFE2_Dprime <- plot_stat(DFE2, "Dprime")
DFE2_D <- plot_stat(DFE2, "D")
DFE2_tajimasd <- plot_stat(DFE2, "tajimasd")
DFE2_hapdiv <- plot_stat(DFE2, "hapdiv")
DFE2_inference <- plot_DFEinference(DFE2)


ggarrange(DFE2_inference, DFE2_Dprime, DFE2_tajimasd, DFE2_hapdiv,
                    labels = c("A", "B", "C", "D"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 2,
                    common.legend = TRUE, legend = "right")


fig <- ggarrange(DFE1_tajimasd, DFE1_Dprime, DFE2_tajimasd, DFE2_Dprime, DFE3_tajimasd, DFE3_Dprime,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 3,
                    common.legend = TRUE, legend = "right")
DFE1_thetapi <- plot_stat(DFE1, "thetapi")
DFE2_thetapi <- plot_stat(DFE2, "thetapi")
DFE3_thetapi <- plot_stat(DFE3, "thetapi")

ggarrange(DFE1_thetapi, DFE1_Dprime, DFE2_thetapi, DFE2_Dprime, DFE3_thetapi, DFE3_Dprime,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 3,
                    common.legend = TRUE, legend = "right")

annotate_figure(fig, 
                left = "DFE3      DFE2     DFE1")



index <- which(colnames(summarize_df)=="hapdiv_sd") +1
index2 <- which(colnames(summarize_df)=="rsq_avg") -1

distances_indices <- seq(from = index, to=index2 - 7, by=8)
dprime_indices <- seq(from = index + 6, to=index2 - 1, by=8)
dprime_sd_indices <- seq(from = index + 7, to=index2, by=8)
d_indices <- seq(from = index + 4, to=index2 - 3, by=8)
d_sd_indices <- seq(from = index + 5, to=index2-2, by=8)
distances <- summarize_df[distances_indices]
dprime <- summarize_df[dprime_indices]
dprime_sd <- summarize_df[dprime_sd_indices]
d <- summarize_df[d_indices]
d_sd <- summarize_df[d_sd_indices]
condensed_df <- data.frame(Vector = apply(distances, 1, function(row) as.numeric(unlist(row))))


rbind(summarize_df[1:4], dprime, dprime_sd, distances)

plot <- ggplot(df, aes(x = x, y = y, color = interaction(selfing, dfe))) +
  geom_line() +
  geom_point(aes(shape = interaction(selfing, dfe)), size = 3) +
  geom_errorbar(aes(ymin = y - sd, ymax = y + sd), width = 0.1) +
  labs(x = "X Axis Label", y = "Y Axis Label", color = "Legend Label") +
  scale_color_discrete(name = "Legend Label") +
  scale_shape_manual(name = "Legend Label", values = c(16, 17, 18)) +
  theme_minimal()

plot(unlist(distances[1,]), unlist(dprime[1,]),type = "l", col = "blue")
lines(unlist(distances[2,]), unlist(dprime[2,]), col = "red")

lines(unlist(distances[31,]), unlist(dprime[31,]), col = "purple")
lines(unlist(distances[32,]), unlist(dprime[32,]), col = "pink")


plot_dist <- function(DFE_choice, selfing_choice) {
  df1 <- summarize_df %>% filter(DFE == DFE_choice) %>% filter(selfing==selfing_choice)
  distances_neutral <- as.numeric(df1[distances_indices][1,])
  distances_selected <- as.numeric(df1[distances_indices][2,])
  dprime_neutral <- as.numeric(df1[dprime_indices][1,])
  dprime_selected <- as.numeric(df1[dprime_indices][2,])
  dprime_sd_neutral <- as.numeric(df1[dprime_sd_indices][1,])
  dprime_sd_selected <- as.numeric(df1[dprime_sd_indices][2,])
  data <- tibble(type=c(rep("neutral", length(distances_neutral)), rep("selected", length(distances_selected))), 
    x=c(distances_neutral, distances_selected), y=c(dprime_neutral, dprime_selected), 
    sd=c(dprime_sd_neutral, dprime_sd_selected)) %>% group_by(type)

  color <- c(rep("blue", length(distances_neutral)), rep("red", length(distances_selected)))
  # Create a line plot with error bars
  plt <- ggplot(data, aes(x = x, y = y, group=type, fill=factor(type))) +
    geom_line(aes(group=type), colour=color) +
    geom_errorbar(aes(ymin = y - sd, ymax = y + sd), width = 3, position = position_dodge(0.05), colour=color) +
    labs(title = "Line Plot with Error Bars", x = "X-axis label", y = "Y-axis label")
    #scale_fill_manual(values = c("blue", "red"))  # Custom color scheme
  return(plt)
}
plot_dist_d <- function(DFE_choice, selfing_choice) {
  df1 <- summarize_df %>% filter(DFE == DFE_choice) %>% filter(selfing==selfing_choice)
  distances_neutral <- as.numeric(df1[distances_indices][1,])
  distances_selected <- as.numeric(df1[distances_indices][2,])
  d_neutral <- as.numeric(df1[d_indices][1,])
  d_selected <- as.numeric(df1[d_indices][2,])
  d_sd_neutral <- as.numeric(df1[d_sd_indices][1,])
  d_sd_selected <- as.numeric(df1[d_sd_indices][2,])
  data <- tibble(type=c(rep("neutral", length(distances_neutral)), rep("selected", length(distances_selected))), 
    x=c(distances_neutral, distances_selected), y=c(d_neutral, d_selected), 
    sd=c(d_sd_neutral, d_sd_selected)) %>% group_by(type)

  color <- c(rep("blue", length(distances_neutral)), rep("red", length(distances_selected)))
  # Create a line plot with error bars
  plt <- ggplot(data, aes(x = x, y = y, group=type, fill=factor(type))) +
    geom_line(aes(group=type), colour=color) +
    geom_errorbar(aes(ymin = y - sd, ymax = y + sd), width = 3, position = position_dodge(0.05), colour=color) +
    labs(title = "Line Plot with Error Bars", x = "X-axis label", y = "Y-axis label")
    #scale_fill_manual(values = c("blue", "red"))  # Custom color scheme
  return(plt)
}
plot_dist("DFE1", 0)
#plot_dist_d("DFE1", 0)

plot_dist("DFE1", 50)
plot_dist("DFE1", 99)

plot_dist("DFE2", 0)
plot_dist("DFE2", 50)
plot_dist("DFE2", 99)

plot_dist("DFE3", 0)
plot_dist("DFE3", 50)
plot_dist("DFE3", 80)
plot_dist("DFE3", 90)
plot_dist("DFE3", 99)
plot_dist_d("DFE3", 0)
plot_dist_d("DFE3", 50)


plot_dist_d("DFE3", 99)

#write.csv(summarize_df, file="/nas/longleaf/home/adaigle/DFESelfing/pylibseq/freq1_5.csv")
write.csv(summarize_df, file="/nas/longleaf/home/adaigle/DFESelfing/pylibseq/freq1_5_500bin.csv")

