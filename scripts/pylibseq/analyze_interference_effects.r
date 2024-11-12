rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)

#frequency filtered results
df <- read.table("/nas/longleaf/home/adaigle/DFESelfing/scripts/pylibseq/allresults_5000_0.stats", header=T)
#df <- read.table("/nas/longleaf/home/adaigle/DFESelfing/scripts/pylibseq/lowrec/allresults_window_5000_0BinSize5000.stats", header=T)

#non frequency filtered results
#df <- read.table("/nas/longleaf/home/adaigle/DFESelfing/pylibseq/allresults_50000.stats", header=T)

gammabeta <- read.csv("/nas/longleaf/home/adaigle/DFESelfing_unused_analyses/pylibseq/gammabeta.csv")
df <- df %>% mutate(
    selfing = as.numeric(str_extract(filename, "(?<=eqm_selfing)\\d+")),
    DFE = str_extract(filename, "(DFE)\\d+")
    ) %>% group_by(selfing, DFE)
gammabeta <- gammabeta %>% select(B, empirical_Ne, DFE, selfing, gamma, b, GammaZero.negGmean, GammaZero.negGshape, selfing_Ne, newNE) %>%
  group_by(selfing, DFE)

gammabeta_propersd <- gammabeta
gammabeta_propersd$selfing_B <- gammabeta_propersd$B / (gammabeta_propersd$selfing_Ne / 5000)
gammabeta_propersd$s_dfealpha <- -gammabeta_propersd$gamma / (2*gammabeta_propersd$selfing_B * 5000)
gammabeta_propersd$s_grapes <- -gammabeta_propersd$GammaZero.negGmean / (2*gammabeta_propersd$selfing_B * 5000)
gammabeta_propersd <- gammabeta_propersd %>% select(
  DFE, selfing, s_dfealpha, b, 
    s_grapes, GammaZero.negGshape) %>% 
  group_by(DFE, selfing) %>%
  summarize(across(where(is.numeric), 
    list(avg = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE))))

join_df <- left_join(df, gammabeta, multiple="all") 

join_df$thetadelta <- 1 - (join_df$thetapi / join_df$thetaw)

join_df$selfing_B <- join_df$B / (join_df$selfing_Ne / 5000)
join_df$s_dfealpha <- -join_df$gamma / (2*join_df$selfing_B * 5000)
join_df$s_grapes <- -join_df$GammaZero.negGmean / (2*join_df$selfing_B * 5000)
join_df$faywuh <- join_df$thetapi - join_df$thetah

#find the error in gamma and beta
truth_df <- data.frame(DFE = c("DFE1", "DFE2", "DFE3"), true_gamma = c(-5, -50, -1000), true_s = c(-5e-04, -5e-03, -0.1), true_beta = c(0.9, 0.6, 0.3))

join_df <- merge(join_df, truth_df, by = "DFE")
join_df$gammaerror <- abs((-join_df$gamma - join_df$true_gamma ) / (join_df$true_gamma))
join_df$serror <- abs((join_df$s_dfealpha - join_df$true_s ) / (join_df$true_s))
join_df$betaerror <- abs((join_df$b - join_df$true_beta ) / (join_df$true_beta))

summarize_df <- join_df %>% ungroup() %>% group_by(filename,output,site,DFE,selfing) %>% 
    summarize(across(where(is.numeric), list(avg = ~mean(., na.rm = TRUE)), .names = "{col}")) %>% #.names makes sure columns aren't renamed to _avg yet
    ungroup() %>% group_by(DFE, selfing,site) %>%
    summarize(across(where(is.numeric), list(avg = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE))))


lowfreqs <- read.csv("/nas/longleaf/home/adaigle/DFESelfing/scripts/pylibseq/freq1_5.csv")
summarize_df$Dprime_avg <- lowfreqs$Dprime_avg
summarize_df$Dprime_sd <- lowfreqs$Dprime_sd
summarize_df <- summarize_df %>% arrange(desc(site))
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
DFE2$gammaerror_avg <- abs(50 - DFE2$gamma_avg) / 50
DFE2$gammaerror_sd <- 0
DFE2$betaerror_avg <- abs(.5 - DFE2$b_avg) / .5
DFE2$betaerror_sd <- 0

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

quick_summary_all <- summarize_df %>% select(DFE,selfing,site,Dprime_avg,thetadelta_avg) %>% arrange(selfing)
#write.csv(quick_summary_all, file = "/nas/longleaf/home/adaigle/DFESelfing/pylibseq/quicksumall.csv", quote=F)


plot_stat <- function(DFE, stat) {
  ggplot(DFE, aes(x = selfing, y = get(paste0(stat, "_avg")), fill=site)) +
  geom_errorbar(aes(ymin = get(paste0(stat, "_avg")) - get(paste0(stat, "_sd")), ymax = get(paste0(stat, "_avg")) + get(paste0(stat, "_sd")),color=site), 
    width = 5, position=position_dodge(0.1)) +
  geom_line(aes(color=site)) +
  geom_point() + 
  labs(x = "% Selfing",
       y = paste0(stat, "_avg")) +
  scale_color_manual(values=c("#619CFF", "#F8766D")) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
  axis.title.x=element_text(size=12),axis.title.y=element_text(size=12), strip.text = element_text(size=15), 
  plot.title= element_text(size=0), legend.title = element_text(size=15), legend.text = element_text(size=15))
}

plot_stat_dfealpha <- function(DFE, stat) {
  ggplot(DFE, aes(x = selfing, y = get(paste0(stat, "_avg")), fill=site)) +
  geom_errorbar(aes(ymin = get(paste0(stat, "_avg")) - get(paste0(stat, "_sd")), ymax = get(paste0(stat, "_avg")) + get(paste0(stat, "_sd")),color=site), 
    width = 5, position=position_dodge(0.1)) +
  geom_line(aes(color=site)) +
  geom_point() + 
  labs(x = "% Selfing",
       y = paste0(stat, "_avg")) +
  scale_color_manual(values=c("#F8766D", "purple")) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
  axis.title.x=element_text(size=12),axis.title.y=element_text(size=12), strip.text = element_text(size=15), 
  plot.title= element_text(size=0), legend.title = element_text(size=15), legend.text = element_text(size=15))
}

plot_DFEinference <- function(DFE) {
  ggplot(DFE, aes(x = selfing)) +

  geom_line(aes(y = B_avg), colour = "#F8766D") +
  geom_point(aes(y = B_avg)) 
  #geom_errorbar(aes(ymin = B_avg - B_sd, ymax = B_avg + B_sd,color="#F8766D"), 
  #  width = 5, position=position_dodge(0.1)) 
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
DFE3_Dprime <- plot_stat(DFE3, "Dprime") + labs(y = expression(italic("D'"))) + ylim(c(-0.975, -0.850))
DFE3_D <- plot_stat(DFE3, "D")
DFE3_tajimasd <- plot_stat(DFE3, "tajimasd")
DFE3_hapdiv <- plot_stat(DFE3, "hapdiv")
DFE3_inference <- plot_DFEinference(DFE3)
DFE3_brian <- plot_stat(DFE3, "thetadelta") + labs(y = expression(paste(Delta,theta)[w])) + ylim(c(0, 0.75))
DFE3_rsq <- plot_stat(DFE3, "rsq")
DFE3_thetah <- plot_stat(DFE3, "thetah")
DFE3_gamma <- plot_stat(DFE3, "b")


ggarrange(DFE3_thetah, DFE3_Dprime, DFE3_tajimasd, DFE3_hapdiv,
                    labels = c("A", "B", "C", "D"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 2,
                    common.legend = TRUE, legend = "right")

 plot_stat(DFE1, "hprime")
DFE1_Dprime <- plot_stat(DFE1, "Dprime") + labs(y = expression(italic("D'"))) + ylim(c(-0.975, -0.850))
DFE1_tajimasd <- plot_stat(DFE1, "tajimasd")
DFE1_hapdiv <- plot_stat(DFE1, "hapdiv")
DFE1_inference <- plot_DFEinference(DFE1)
DFE1_brian <- plot_stat(DFE1, "thetadelta") + labs(y = expression(paste(Delta,theta)[w])) + ylim(c(0, 0.75))
DFE1_rsq <- plot_stat(DFE1, "rsq")

#DFE1$faywuh_avg <- DFE1$thetapi_avg - DFE1$thetah_avg
#DFE1$faywuh_sd <- 0
#DFE1_faywuh <- plot_stat(DFE1, "faywuh")
#DFE1$faywud_avg <- DFE1$thetaw_avg - DFE1$thetapi_avg
#DFE1$faywud_sd <- 0
#DFE1_faywud <- plot_stat(DFE1, "faywud")
#
#DFE2$faywuh_avg <- DFE2$thetapi_avg - DFE2$thetah_avg
#DFE2$faywuh_sd <- 0
#DFE2_faywuh <- plot_stat(DFE2, "faywuh")
#DFE2$faywud_avg <- DFE2$thetaw_avg - DFE2$thetapi_avg
#DFE2$faywud_sd <- 0
#DFE2_faywud <- plot_stat(DFE2, "faywud")
#
#DFE3$faywuh_avg <- DFE3$thetapi_avg - DFE3$thetah_avg
#DFE3$faywuh_sd <- 0
#DFE3_faywuh <- plot_stat(DFE3, "faywuh")
#DFE3$faywud_avg <- DFE3$thetaw_avg - DFE3$thetapi_avg
#DFE3$faywud_sd <- 0
#DFE3_faywud <- plot_stat(DFE3, "faywud")

ggarrange(DFE1_inference, DFE1_Dprime, DFE1_tajimasd, DFE1_hapdiv,
                    labels = c("A", "B", "C", "D"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 2,
                    common.legend = TRUE, legend = "right")
DFE2_Dprime <- plot_stat(DFE2, "Dprime") + labs(y = expression(italic("D'"))) + ylim(c(-0.975, -0.850))
DFE2_D <- plot_stat(DFE2, "D")
DFE2_tajimasd <- plot_stat(DFE2, "tajimasd")
DFE2_hapdiv <- plot_stat(DFE2, "hapdiv")
DFE2_inference <- plot_DFEinference(DFE2)
DFE2_brian <- plot_stat(DFE2, "thetadelta") + labs(y = expression(paste(Delta,theta)[w])) + ylim(c(0, 0.75))
DFE2_rsq <- plot_stat(DFE2, "rsq")
DFE2_thetah <- plot_stat(DFE2, "thetah")
DFE2_numSing <- plot_stat(DFE2, "numSing")
DFE2_gamma <- plot_stat(DFE2, "gamma")
DFE2_numSing <- plot_stat(DFE2, "numSing")
DFE2_selfing_B <- plot_stat(DFE2, "selfing_B")
DFE2_pi <- plot_stat(DFE3, "thetapi")

ggarrange(DFE2_inference, DFE2_selfing_B, DFE2_gamma, DFE2_thetah, DFE2_Dprime, DFE2_tajimasd, DFE2_hapdiv, DFE2_brian, DFE2_numSing,
                    labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 3, nrow = 3,
                    common.legend = TRUE, legend = "right")


figure4 <- ggarrange(DFE1_brian, DFE1_Dprime, DFE2_brian, DFE2_Dprime, DFE3_brian, DFE3_Dprime,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 3,
                    common.legend = TRUE, legend = "bottom", hjust=0)

figure4_rsq <- ggarrange(DFE1_brian, DFE1_rsq, DFE2_brian, DFE2_rsq, DFE3_brian, DFE3_rsq,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 3,
                    common.legend = TRUE, legend = "bottom", hjust=0)


df_long_s <- gammabeta_propersd %>%
  pivot_longer(cols = c(s_dfealpha_avg, s_dfealpha_sd, s_grapes_avg, s_grapes_sd), 
    names_to = c(".value", "type"), 
    names_pattern = "(s_dfealpha|s_grapes)_(.*)", 
    values_to = "value"
  ) %>% pivot_longer(cols = c(s_dfealpha, s_grapes), 
               names_to = "measurement", 
               values_to = "s") %>%
  mutate(site = ifelse(measurement == "s_dfealpha", "DFE-alpha", "GRAPES")) %>%
  select(-measurement) %>%
  pivot_wider(names_from = type, values_from = s)


df_long_beta <- gammabeta_propersd %>%
  pivot_longer(cols = c(b_avg, b_sd, GammaZero.negGshape_avg, GammaZero.negGshape_sd), 
    names_to = c(".value", "type"), 
    names_pattern = "(b|GammaZero.negGshape)_(.*)", 
    values_to = "value"
  ) %>% pivot_longer(cols = c(b, GammaZero.negGshape), 
               names_to = "measurement", 
               values_to = "b") %>%
  mutate(site = ifelse(measurement == "b", "DFE-alpha", "GRAPES")) %>%
  select(-measurement) %>%
  pivot_wider(names_from = type, values_from = b)


df_long_beta <- gammabeta_propersd %>%
  pivot_longer(cols = c(b, GammaZero.negGshape), 
               names_to = "measurement", 
               values_to = "beta") %>%
  mutate(site = ifelse(measurement == "b", "DFE-alpha", "GRAPES")) %>%
  select(-measurement) # Remove the 'measurement' column if it's no longer needed

gammabeta_propersd$site<-"DFE-alpha"
DFE1_gammabeta_propersd <- gammabeta_propersd %>% filter(DFE=="DFE1")
DFE2_gammabeta_propersd <- gammabeta_propersd %>% filter(DFE=="DFE2")
DFE3_gammabeta_propersd <- gammabeta_propersd %>% filter(DFE=="DFE3")

DFE1_df_long_s <- df_long_s %>% filter(DFE=="DFE1")
DFE2_df_long_s <- df_long_s %>% filter(DFE=="DFE2")
DFE3_df_long_s <- df_long_s %>% filter(DFE=="DFE3")

DFE1_df_long_beta <- df_long_beta %>% filter(DFE=="DFE1")
DFE2_df_long_beta <- df_long_beta %>% filter(DFE=="DFE2")
DFE3_df_long_beta <- df_long_beta %>% filter(DFE=="DFE3")


plot_stat_dfealphagrapes <- function(df, stat) {
# Plot using ggplot2
ggplot(df, aes(x = selfing, y = avg, fill = site)) +
  geom_errorbar(aes(
    ymin = avg - sd,
    ymax = avg + sd,
    color = site
  ), width = 5, position = position_dodge(1)) +
  geom_line(aes(color = site), position = position_dodge(0.1)) +
geom_point(position = position_dodge(1), aes(fill = site), size = 2.5, shape = 21, stroke=1) +
  labs(
    x = "% Selfing",
    y = paste0(stat, " avg")
  ) +
  scale_color_manual(values = c("#F8766D", "purple")) +
  scale_fill_manual(values = c("#F8766D", "purple")) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 15),
    plot.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  )
}

DFE1_s_dfealpha <-plot_stat_dfealphagrapes(DFE1_df_long_s, "s") +
  geom_hline(aes(yintercept = -5e-04, linetype = "truth"), color = "black", size = 0.5) +
  scale_linetype_manual(name = "Legend", values = "dashed")

DFE2_s_dfealpha <- plot_stat_dfealphagrapes(DFE2_df_long_s, "s") +
  geom_hline(aes(yintercept = -0.005, linetype = "truth"), color = "black", size = 0.5) +
  scale_linetype_manual(name = "Legend", values = "dashed")

DFE3_s_dfealpha <- plot_stat_dfealphagrapes(DFE3_df_long_s, "s") +
  geom_hline(aes(yintercept = -0.1, linetype = "truth"), color = "black", size = 0.5) +
  scale_linetype_manual(name = "Legend", values = "dashed")

DFE1_beta <- plot_stat_dfealphagrapes(DFE1_df_long_beta, "b") +
  geom_hline(aes(yintercept = 0.9, linetype = "truth"), color = "black", size = 0.5) +
  scale_linetype_manual(name = "Legend", values = "dashed")

DFE2_beta <- plot_stat_dfealphagrapes(DFE2_df_long_beta, "b") +
  geom_hline(aes(yintercept = 0.5, linetype = "truth"), color = "black", size = 0.5) +
  scale_linetype_manual(name = "Legend", values = "dashed")

DFE3_beta <- plot_stat_dfealphagrapes(DFE3_df_long_beta, "b") +
  geom_hline(aes(yintercept = 0.3, linetype = "truth"), color = "black", size = 0.5) +
  scale_linetype_manual(name = "Legend", values = "dashed")

dfealpha_s_and_b <- ggarrange(DFE1_s_dfealpha, DFE2_s_dfealpha, DFE3_s_dfealpha, DFE1_beta, DFE2_beta, DFE3_beta,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 3, nrow = 2,
                    common.legend = TRUE, legend = "bottom", hjust=0)


saveRDS(DFE1_s_dfealpha, "/nas/longleaf/home/adaigle/DFESelfing/scripts/summary_figures/DFE1_s_dfealpha.rds")
saveRDS(DFE2_s_dfealpha, "/nas/longleaf/home/adaigle/DFESelfing/scripts/summary_figures/DFE2_s_dfealpha.rds")
saveRDS(DFE3_s_dfealpha, "/nas/longleaf/home/adaigle/DFESelfing/scripts/summary_figures/DFE3_s_dfealpha.rds")
saveRDS(DFE1_beta, "/nas/longleaf/home/adaigle/DFESelfing/scripts/summary_figures/DFE1_beta.rds")
saveRDS(DFE2_beta, "/nas/longleaf/home/adaigle/DFESelfing/scripts/summary_figures/DFE2_beta.rds")
saveRDS(DFE3_beta, "/nas/longleaf/home/adaigle/DFESelfing/scripts/summary_figures/DFE3_beta.rds")

DFE1_s_grapes <- plot_stat(DFE1, "s_grapes") +
  geom_hline(aes(yintercept = -5e-04, linetype = "truth"), color = "black", size = 0.5) +
  scale_linetype_manual(name = "Legend", values = "dashed")

DFE2_s_grapes <- plot_stat(DFE2, "s_grapes") +
  geom_hline(aes(yintercept = -0.005, linetype = "truth"), color = "black", size = 0.5) +
  scale_linetype_manual(name = "Legend", values = "dashed")

DFE3_s_grapes <- plot_stat(DFE3, "s_grapes") +
  geom_hline(aes(yintercept = -0.1, linetype = "truth"), color = "black", size = 0.5) +
  scale_linetype_manual(name = "Legend", values = "dashed")

DFE1_gbeta <- plot_stat(DFE1, "GammaZero.negGshape") +
  geom_hline(aes(yintercept = 0.9, linetype = "truth"), color = "black", size = 0.5) +
  scale_linetype_manual(name = "Legend", values = "dashed")

DFE2_gbeta <- plot_stat(DFE2, "GammaZero.negGshape") +
  geom_hline(aes(yintercept = 0.5, linetype = "truth"), color = "black", size = 0.5) +
  scale_linetype_manual(name = "Legend", values = "dashed")

DFE3_gbeta <- plot_stat(DFE3, "GammaZero.negGshape") +
  geom_hline(aes(yintercept = 0.3, linetype = "truth"), color = "black", size = 0.5) +
  scale_linetype_manual(name = "Legend", values = "dashed")

grapes_s_and_b <- ggarrange(DFE1_s_grapes, DFE2_s_grapes, DFE3_s_grapes, DFE1_gbeta, DFE2_gbeta, DFE3_gbeta,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 3, nrow = 2,
                    common.legend = TRUE, legend = "bottom", hjust=0)
#figure4_D <- ggarrange(DFE1_brian, DFE1_D, DFE2_brian, DFE2_D, DFE3_brian, DFE3_D,
#                    labels = c("A", "B", "C", "D", "E", "F"),
#                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
#                    ncol = 2, nrow = 3,
#                    common.legend = TRUE, legend = "bottom", hjust=0)
#
#DvsDprime <- ggarrange(DFE1_Dprime, DFE1_D, DFE2_Dprime, DFE2_D, DFE3_Dprime, DFE3_D,
#                    labels = c("A", "B", "C", "D", "E", "F"),
#                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
#                    ncol = 2, nrow = 3,
#                    common.legend = TRUE, legend = "bottom", hjust=0)
#annotate_figure(fig, 
#                left = "DFE3      DFE2     DFE1")

#ggsave("/nas/longleaf/home/adaigle/DFESelfing/figures_for_publication/figure9.pdf", plot = figure9, width = 8.5, height = 8.5, dpi = 600)
ggsave("/nas/longleaf/home/adaigle/DFESelfing/figures_for_publication/figure4.svg", plot = figure4, width = 8.5, height = 8.5, dpi = 300)


figure4_poster <- ggarrange(DFE1_brian, DFE1_Dprime, DFE2_brian, DFE2_Dprime, DFE3_brian, DFE3_Dprime,
                    #labels = c("A", "B", "C", "D", "E", "F"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 3,
                    common.legend = TRUE, legend = "bottom", hjust=0)
ggsave("/nas/longleaf/home/adaigle/DFESelfing/figures_for_publication/figure4_poster.png", plot = figure4_poster, width = 6.5, height = 9, dpi = 300, unit="in")


regress_with_B <- function(DFE, stat) {
  ggplot(DFE, aes(x = B_avg, y = get(paste0(stat, "_avg")), fill=site)) +
  geom_errorbar(aes(ymin = get(paste0(stat, "_avg")) - get(paste0(stat, "_sd")), ymax = get(paste0(stat, "_avg")) + get(paste0(stat, "_sd")),color=site), 
    width = .1, position=position_dodge(0)) +
  geom_line(aes(color=site)) +
  geom_point() + 
  labs(x = "B value",
       y = paste0(stat, "_avg")) +
  scale_color_manual(values=c("#619CFF", "#F8766D")) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
  axis.title.x=element_text(size=12),axis.title.y=element_text(size=12), strip.text = element_text(size=15), 
  plot.title= element_text(size=0), legend.title = element_text(size=15), legend.text = element_text(size=15))
}

B_DFE2_Dprime <- regress_with_B(DFE2, "Dprime") + labs(y = expression(italic("D'"))) + ylim(c(-0.975, -0.850))
B_DFE2_D <- regress_with_B(DFE2, "D")
B_DFE2_tajimasd <- regress_with_B(DFE2, "tajimasd")
B_DFE2_hapdiv <- regress_with_B(DFE2, "hapdiv")
B_DFE2_brian <- regress_with_B(DFE2, "thetadelta") + labs(y = expression(paste(Delta,theta)[w])) + ylim(c(0, 0.75))
B_DFE2_rsq <- regress_with_B(DFE2, "rsq")
B_DFE2_thetah <- regress_with_B(DFE2, "thetah")
B_DFE2_numSing <- regress_with_B(DFE2, "numSing")
B_DFE2_selfing_B <- regress_with_B(DFE2, "selfing_B")
B_DFE2_B <- regress_with_B(DFE2, "B")
B_DFE2_gamma <- regress_with_B(DFE2, "gammaerror")


ggarrange(B_DFE2_B, B_DFE2_selfing_B, B_DFE2_gamma, B_DFE2_thetah, B_DFE2_Dprime, B_DFE2_tajimasd, B_DFE2_hapdiv, B_DFE2_brian, B_DFE2_numSing,
                    labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 3, nrow = 3,
                    common.legend = TRUE, legend = "right")


regress_with_gamma <- function(DFE, stat) {
  ggplot(DFE, aes(x = gammaerror_avg, y = get(paste0(stat, "_avg")), fill=site)) +
  #geom_errorbar(aes(ymin = get(paste0(stat, "_avg")) - get(paste0(stat, "_sd")), ymax = get(paste0(stat, "_avg")) + get(paste0(stat, "_sd")),color=site), 
  #  width = .1, position=position_dodge(0)) +
  geom_line(aes(color=site)) +
  geom_point() + 
  labs(x = "Gamma value",
       y = paste0(stat, "_avg")) +
  scale_color_manual(values=c("#619CFF", "#F8766D")) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
  axis.title.x=element_text(size=12),axis.title.y=element_text(size=12), strip.text = element_text(size=15), 
  plot.title= element_text(size=0), legend.title = element_text(size=15), legend.text = element_text(size=15))
}

gamma_DFE2_Dprime <- regress_with_gamma(DFE2, "Dprime") + labs(y = expression(italic("D'"))) + ylim(c(-0.975, -0.850))
gamma_DFE2_D <- regress_with_gamma(DFE2, "D")
gamma_DFE2_tajimasd <- regress_with_gamma(DFE2, "tajimasd")
gamma_DFE2_hapdiv <- regress_with_gamma(DFE2, "hapdiv")
gamma_DFE2_brian <- regress_with_gamma(DFE2, "thetadelta") + labs(y = expression(paste(Delta,theta)[w])) + ylim(c(0, 0.75))
gamma_DFE2_rsq <- regress_with_gamma(DFE2, "rsq")
gamma_DFE2_thetah <- regress_with_gamma(DFE2, "thetah")
gamma_DFE2_numSing <- regress_with_gamma(DFE2, "numSing")
gamma_DFE2_selfing_B <- regress_with_gamma(DFE2, "selfing_B")
gamma_DFE2_B <- regress_with_gamma(DFE2, "B")
gamma_DFE2_gamma <- regress_with_gamma(DFE2, "gamma")

ggarrange(gamma_DFE2_B, gamma_DFE2_selfing_B, gamma_DFE2_gamma, gamma_DFE2_thetah, gamma_DFE2_Dprime, gamma_DFE2_tajimasd, gamma_DFE2_hapdiv, gamma_DFE2_brian, gamma_DFE2_numSing,
                    labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 3, nrow = 3,
                    common.legend = TRUE, legend = "right")


regress_with_beta <- function(DFE, stat) {
  ggplot(DFE, aes(x =betaerror_avg, y = get(paste0(stat, "_avg")), fill=site)) +
  #geom_errorbar(aes(ymin = get(paste0(stat, "_avg")) - get(paste0(stat, "_sd")), ymax = get(paste0(stat, "_avg")) + get(paste0(stat, "_sd")),color=site), 
  #  width = .01, position=position_dodge(0)) +
  geom_line(aes(color=site)) +
  geom_point() + 
  labs(x = "beta error",
       y = paste0(stat, "_avg")) +
  scale_color_manual(values=c("#619CFF", "#F8766D")) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
  axis.title.x=element_text(size=12),axis.title.y=element_text(size=12), strip.text = element_text(size=15), 
  plot.title= element_text(size=0), legend.title = element_text(size=15), legend.text = element_text(size=15))
}

beta_DFE2_Dprime <- regress_with_beta(DFE2, "Dprime") + labs(y = expression(italic("D'"))) + ylim(c(-0.975, -0.850))
beta_DFE2_D <- regress_with_beta(DFE2, "D")
beta_DFE2_tajimasd <- regress_with_beta(DFE2, "tajimasd")
beta_DFE2_hapdiv <- regress_with_beta(DFE2, "hapdiv")
beta_DFE2_brian <- regress_with_beta(DFE2, "thetadelta") + labs(y = expression(paste(Delta,theta)[w])) + ylim(c(0, 0.75))
beta_DFE2_rsq <- regress_with_beta(DFE2, "rsq")
beta_DFE2_thetah <- regress_with_beta(DFE2, "thetah")
beta_DFE2_numSing <- regress_with_beta(DFE2, "numSing")
beta_DFE2_selfing_B <- regress_with_beta(DFE2, "selfing_B")
beta_DFE2_B <- regress_with_beta(DFE2, "B")
beta_DFE2_gamma <- regress_with_beta(DFE2, "gamma")
beta_DFE2_hprime <- regress_with_beta(DFE2, "hprime")

ggarrange(beta_DFE2_B, beta_DFE2_selfing_B, beta_DFE2_gamma, beta_DFE2_thetah, beta_DFE2_Dprime, beta_DFE2_tajimasd, beta_DFE2_hapdiv, beta_DFE2_brian, beta_DFE2_numSing,
                    labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 3, nrow = 3,
                    common.legend = TRUE, legend = "right")

neu_by_sel_DFE2 <- DFE2 %>%
  group_by(selfing) %>%
  summarise(across(starts_with("posn_"):starts_with("selfing_B_"), ~ .[site == "selected"] / .[site == "neutral"], .names = "ratio_{col}"))

neu_by_sel_DFE2 <- left_join(DFE2, neu_by_sel_DFE2)

regress_with_beta(neu_by_sel_DFE2, "ratio_Dprime")
regress_with_beta(neu_by_sel_DFE2, "ratio_thetapi")
regress_with_beta(neu_by_sel_DFE2, "ratio_hprime")
regress_with_beta(neu_by_sel_DFE2, "ratio_thetah")
regress_with_beta(neu_by_sel_DFE2, "ratio_hapdiv")
regress_with_beta(neu_by_sel_DFE2, "ratio_tajimasd")
regress_with_beta(neu_by_sel_DFE2, "ratio_D")
regress_with_beta(neu_by_sel_DFE2, "ratio_rsq")


regress_with_beta(neu_by_sel_DFE2, "hprime")
regress_with_beta(DFE2, "tajimasd")
regress_with_beta(DFE2, "tajimasd")
regress_with_beta(DFE2, "tajimasd")

plot_stat(neu_by_sel_DFE2, "thetapi")
regress_with_B(neu_by_sel_DFE2, "ratio_thetapi")
regress_with_gamma(neu_by_sel_DFE2, "ratio_thetapi")
regress_with_beta(neu_by_sel_DFE2, "ratio_thetapi")


summarize_df$siteDFE <- interaction(summarize_df$site, summarize_df$DFE)

plot_stat_allDFEs <- function(df, stat) {
  ggplot(df, aes(x = selfing, y = get(paste0(stat, "_avg")), fill=DFE)) +
  geom_errorbar(aes(ymin = get(paste0(stat, "_avg")) - get(paste0(stat, "_sd")), ymax = get(paste0(stat, "_avg")) + get(paste0(stat, "_sd")),color=DFE), 
    width = 5, position=position_dodge(0.1)) +
  geom_line(aes(color=DFE)) +
  geom_point() + 
  labs(x = "% Selfing",
       y = paste0(stat, "_avg")) +
  scale_color_manual(values=c("blue", "#F8766D", "#619CFF", "red", "purple", "orange")) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
  axis.title.x=element_text(size=12),axis.title.y=element_text(size=12), strip.text = element_text(size=15), 
  plot.title= element_text(size=0), legend.title = element_text(size=15), legend.text = element_text(size=15))
}

selfing_B <- plot_stat_allDFEs(summarize_df, "selfing_B")
B <- plot_stat_allDFEs(summarize_df, "B")
s_dfealpha <- plot_stat_allDFEs(summarize_df, "s_dfealpha")

dfealpha_beta <- plot_stat_allDFEs(summarize_df, "b")
gammaerror <- plot_stat_allDFEs(summarize_df, "gammaerror")
betaerror <- plot_stat_allDFEs(summarize_df, "betaerror")
betaerror <- plot_stat_allDFEs(summarize_df, "serror")

thetapi <- plot_stat_allDFEs(summarize_df, "thetapi")
thetaw <- plot_stat_allDFEs(summarize_df, "thetaw")
thetah <- plot_stat_allDFEs(summarize_df, "thetah")

tajimasd <- plot_stat_allDFEs(summarize_df, "tajimasd")
thetadelta <- plot_stat_allDFEs(summarize_df, "thetadelta")
faywuh <- plot_stat_allDFEs(summarize_df, "faywuh")

rsq <- plot_stat_allDFEs(summarize_df, "rsq")
D <- plot_stat_allDFEs(summarize_df, "D")
Dprime <- plot_stat_allDFEs(summarize_df, "Dprime")
hapdiv <- plot_stat_allDFEs(summarize_df, "hapdiv")

ggarrange(selfing_B, B, s_dfealpha, dfealpha_beta, gammaerror, betaerror, thetapi, thetaw, thetah, tajimasd, thetadelta, faywuh, rsq, Dprime, hapdiv,
                    #labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 3, nrow=5,
                    common.legend = TRUE, legend = "right")

neu_by_sel <- summarize_df %>%
  group_by(selfing, DFE) %>%
  summarise(across(starts_with("posn_"):starts_with("faywuh"), ~ .[site == "selected"] / .[site == "neutral"], .names = "ratio_{col}"))

gammaerror_vec <- summarize_df %>% filter(site=="selected") %>% pull(gammaerror_avg)
betaerror_vec <- summarize_df %>% filter(site=="selected") %>% pull(betaerror_avg)
serror_vec <- summarize_df %>% filter(site=="selected") %>% pull(serror_avg)

regress_with_gammaerror <- function(df, stat) {
  ggplot(df, aes(x = serror_vec, y = get(paste0("ratio_", stat, "_avg")), fill=DFE)) +
  #geom_errorbar(aes(ymin = get(paste0(stat, "_avg")) - get(paste0(stat, "_sd")), ymax = get(paste0(stat, "_avg")) + get(paste0(stat, "_sd")),color=site), 
  #  width = .1, position=position_dodge(0)) +
  geom_line(aes(color=DFE)) +
  geom_point() + 
  labs(x = "Error",
       y = paste0(stat, "_avg")) +
  scale_color_manual(values=c("#619CFF", "#F8766D", "green")) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
  axis.title.x=element_text(size=12),axis.title.y=element_text(size=12), strip.text = element_text(size=15), 
  plot.title= element_text(size=0), legend.title = element_text(size=15), legend.text = element_text(size=15))
}

gammaerror_selfing_B <- regress_with_gammaerror(neu_by_sel, "selfing_B")
gammaerror_B <- regress_with_gammaerror(neu_by_sel, "B")
gammaerror_s_dfealpha <- regress_with_gammaerror(neu_by_sel, "s_dfealpha")
gammaerror_dfealpha_beta <- regress_with_gammaerror(neu_by_sel, "b")
gammaerror_gammaerror <- regress_with_gammaerror(neu_by_sel, "gammaerror")
gammaerror_betaerror <- regress_with_gammaerror(neu_by_sel, "betaerror")
gammaerror_thetapi <- regress_with_gammaerror(neu_by_sel, "thetapi")
gammaerror_thetaw <- regress_with_gammaerror(neu_by_sel, "thetaw")
gammaerror_thetah <- regress_with_gammaerror(neu_by_sel, "thetah")
gammaerror_tajimasd <- regress_with_gammaerror(neu_by_sel, "tajimasd")
gammaerror_thetadelta <- regress_with_gammaerror(neu_by_sel, "thetadelta")
gammaerror_faywuh <- regress_with_gammaerror(neu_by_sel, "faywuh")
gammaerror_rsq <- regress_with_gammaerror(neu_by_sel, "rsq")
gammaerror_D <- regress_with_gammaerror(neu_by_sel, "D")
gammaerror_Dprime <- regress_with_gammaerror(neu_by_sel, "Dprime")
gammaerror_hapdiv <- regress_with_gammaerror(neu_by_sel, "hapdiv")

ggarrange(gammaerror_thetapi, gammaerror_thetaw, gammaerror_thetah, gammaerror_tajimasd, gammaerror_thetadelta, gammaerror_faywuh, gammaerror_rsq, gammaerror_Dprime, gammaerror_hapdiv,
                    #labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 3, nrow=3,
                    common.legend = TRUE, legend = "right")

plot_stat_allDFEs(neu_by_sel, "ratio_rsq")

regress_ratio_stat_gammaerror_allDFEs <- function(df, stat) {
  ggplot(df, aes(x = -1*gammaerror_avg, y = get(paste0(stat, "_avg")), fill=DFE)) +
  #geom_errorbar(aes(ymin = get(paste0(stat, "_avg")) - get(paste0(stat, "_sd")), ymax = get(paste0(stat, "_avg")) + get(paste0(stat, "_sd")),color=DFE), 
  #  width = 0.01, position=position_dodge(0.001)) +
  geom_line(aes(color=DFE)) +
  geom_point() + 
  labs(x = "Gamma error",
       y = paste0(stat, "_avg")) +
  scale_color_manual(values=c("blue", "#F8766D", "#619CFF", "red", "purple", "orange")) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
  axis.title.x=element_text(size=12),axis.title.y=element_text(size=12), strip.text = element_text(size=15), 
  plot.title= element_text(size=0), legend.title = element_text(size=15), legend.text = element_text(size=15))
}

neu_by_sel <- summarize_df %>%
  group_by(selfing, DFE) %>%
  summarise(across(starts_with("posn_"):starts_with("faywuh"), ~ .[site == "selected"] / .[site == "neutral"], .names = "ratio_{col}"))
neu_by_sel <- left_join(summarize_df, neu_by_sel)

gammaerror_thetapi <- regress_ratio_stat_gammaerror_allDFEs(neu_by_sel, "ratio_thetapi")
gammaerror_thetah <- regress_ratio_stat_gammaerror_allDFEs(neu_by_sel, "ratio_thetah")
gammaerror_thetadelta <- regress_ratio_stat_gammaerror_allDFEs(neu_by_sel, "ratio_thetadelta")
gammaerror_numsing <- regress_ratio_stat_gammaerror_allDFEs(neu_by_sel, "ratio_numSing")
gammaerror_hapdiv <- regress_ratio_stat_gammaerror_allDFEs(neu_by_sel, "ratio_hapdiv")

summarize_df_selected <- summarize_df %>% filter(site=="selected")
gammaerror_faywuh <- regress_ratio_stat_gammaerror_allDFEs(summarize_df_selected, "faywuh")
gammaerror_Dprime <- regress_ratio_stat_gammaerror_allDFEs(summarize_df_selected, "Dprime")
gammaerror_rsq <- regress_ratio_stat_gammaerror_allDFEs(summarize_df_selected, "rsq")
gammaerror_D <- regress_ratio_stat_gammaerror_allDFEs(summarize_df_selected, "D")
gammaerror_betaerror <- regress_ratio_stat_gammaerror_allDFEs(summarize_df_selected, "betaerror")

ggarrange(gammaerror_thetapi, gammaerror_thetah, gammaerror_thetadelta, gammaerror_numsing, gammaerror_hapdiv, gammaerror_faywuh, gammaerror_rsq, gammaerror_D, gammaerror_Dprime,
                    #labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 3, nrow=3,
                    common.legend = TRUE, legend = "right")

ggarrange(gammaerror_thetapi, gammaerror_thetah, gammaerror_numsing, gammaerror_hapdiv, gammaerror_rsq, gammaerror_D,
                    #labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 3, nrow=2,
                    common.legend = TRUE, legend = "right")

regress_ratio_stat_betaerror_allDFEs <- function(df, stat) {
  ggplot(df, aes(x = -1*betaerror_avg, y = get(paste0(stat, "_avg")), fill=DFE)) +
  #geom_errorbar(aes(ymin = get(paste0(stat, "_avg")) - get(paste0(stat, "_sd")), ymax = get(paste0(stat, "_avg")) + get(paste0(stat, "_sd")),color=DFE), 
  #  width = 0.01, position=position_dodge(0.001)) +
  geom_line(aes(color=DFE)) +
  geom_point() + 
  labs(x = "Beta error",
       y = paste0(stat, "_avg")) +
  scale_color_manual(values=c("blue", "#F8766D", "#619CFF", "red", "purple", "orange")) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
  axis.title.x=element_text(size=12),axis.title.y=element_text(size=12), strip.text = element_text(size=15), 
  plot.title= element_text(size=0), legend.title = element_text(size=15), legend.text = element_text(size=15))
}

neu_by_sel <- summarize_df %>%
  group_by(selfing, DFE) %>%
  summarise(across(starts_with("posn_"):starts_with("faywuh"), ~ .[site == "selected"] / .[site == "neutral"], .names = "ratio_{col}"))
neu_by_sel <- left_join(summarize_df, neu_by_sel)

betaerror_thetapi <- regress_ratio_stat_betaerror_allDFEs(neu_by_sel, "ratio_thetapi")
betaerror_thetah <- regress_ratio_stat_betaerror_allDFEs(neu_by_sel, "ratio_thetah")
betaerror_thetadelta <- regress_ratio_stat_betaerror_allDFEs(neu_by_sel, "ratio_thetadelta")
betaerror_numsing <- regress_ratio_stat_betaerror_allDFEs(neu_by_sel, "ratio_numSing")
betaerror_hapdiv <- regress_ratio_stat_betaerror_allDFEs(neu_by_sel, "ratio_hapdiv")

summarize_df_selected <- summarize_df %>% filter(site=="selected")
betaerror_faywuh <- regress_ratio_stat_betaerror_allDFEs(summarize_df_selected, "faywuh")
betaerror_Dprime <- regress_ratio_stat_betaerror_allDFEs(summarize_df_selected, "Dprime")
betaerror_rsq <- regress_ratio_stat_betaerror_allDFEs(summarize_df_selected, "rsq")
betaerror_D <- regress_ratio_stat_betaerror_allDFEs(summarize_df_selected, "D")
betaerror_betaerror <- regress_ratio_stat_betaerror_allDFEs(summarize_df_selected, "betaerror")
betaerror_thetadelta <- regress_ratio_stat_gammaerror_allDFEs(summarize_df_selected, "thetadelta")

ggarrange(betaerror_thetapi, betaerror_thetah, betaerror_numsing, betaerror_hapdiv, betaerror_rsq, betaerror_D,
                    #labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 3, nrow=2,
                    common.legend = TRUE, legend = "right")

# error y axis, U/R on x axis
summarize_df_selected$UoverR <- (2*3.3*1e-9 * 2503000)/((1-(0.01 *summarize_df_selected$selfing/(2-0.01 *summarize_df_selected$selfing))) * 3.12*1e-8 * 2503000)
#summarize_df_selected$UoverR <- (3.3*1e-9 * 2503000)/((1-(0.01 *summarize_df_selected$selfing)) * 3.12 * 1e-8 * 2503000)
regress_stat_U_R_allDFEs <- function(df, stat, y_label) {
  ggplot(df, aes(x = log(UoverR), y = get(paste0(stat, "_avg")), fill=DFE)) +
  #geom_errorbar(aes(ymin = get(paste0(stat, "_avg")) - get(paste0(stat, "_sd")), ymax = get(paste0(stat, "_avg")) + get(paste0(stat, "_sd")),color=DFE), 
  #  width = 0.0001, position=position_dodge(0.000001)) +
  geom_line(aes(color=DFE)) +
  geom_point() + 
  labs(x = expression(frac(U, R)),
       y = y_label) +
  scale_color_manual(values=c("blue", "#F8766D", "#619CFF", "red", "purple", "orange")) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
  axis.title.x=element_text(size=10),axis.title.y=element_text(size=12), strip.text = element_text(size=15), 
  plot.title= element_text(size=0), legend.title = element_blank(), legend.text = element_text(size=15),
  legend.position = "none")
}

U_R_gammaerror <- regress_stat_U_R_allDFEs(summarize_df_selected, "gammaerror", expression(paste("|", italic(gamma), " error", "|"))) + ylim(0, 1)
U_R_betaerror <- regress_stat_U_R_allDFEs(summarize_df_selected, "betaerror", expression(paste("|", italic(beta), " ", "error", "|"))) + ylim(0, 1)
U_R_serror <- regress_stat_U_R_allDFEs(summarize_df_selected, "serror", expression(paste("|", italic(s), " ", "error", "|"))) + ylim(0, 2.2)
error_plot <- ggarrange(
  U_R_gammaerror, U_R_serror, U_R_betaerror,
  ncol = 3, nrow = 1,
  common.legend = F,
  #legend = "right",
  font.label = list(size = 24, color = "black", face = "bold")
)

#ggsave("/nas/longleaf/home/adaigle/DFESelfing/figures_for_publication/revisionfigure1.svg", plot = error_plot, width = 8.5, height = 4.5, dpi = 300)
saveRDS(error_plot, file = "/nas/longleaf/home/adaigle/DFESelfing/scripts/pylibseq/lowrec/selfing_error_plot.RDS")

regress_stat_B_allDFEs <- function(df, stat, y_label) {
  ggplot(df, aes(x = selfing_B_avg, y = get(paste0(stat, "_avg")), fill=DFE)) +
  #geom_errorbar(aes(ymin = get(paste0(stat, "_avg")) - get(paste0(stat, "_sd")), ymax = get(paste0(stat, "_avg")) + get(paste0(stat, "_sd")),color=DFE), 
  #  width = 0.0001, position=position_dodge(0.000001)) +
  geom_line(aes(color=DFE)) +
  geom_point() + 
  labs(x = expression(italic(B)),
       y = y_label) +
  scale_color_manual(values=c("blue", "#F8766D", "#619CFF", "red", "purple", "orange")) + 
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=12),
  axis.title.x=element_text(size=12),axis.title.y=element_text(size=12), strip.text = element_text(size=15), 
  plot.title= element_text(size=0), legend.title = element_blank(), legend.text = element_text(size=15),
  legend.position = "none")
}

B_gammaerror <- regress_stat_B_allDFEs(summarize_df_selected, "gammaerror", expression(paste("|", italic(gamma), " error", "|"))) + ylim(0, 1)
B_betaerror <- regress_stat_B_allDFEs(summarize_df_selected, "betaerror", expression(paste("|", italic(beta), " ", "error", "|"))) + ylim(0, 1)
B_serror <- regress_stat_B_allDFEs(summarize_df_selected, "serror", expression(paste("|", italic(s), " ", "error", "|"))) + ylim(0, 2.2)
error_plot_B <- ggarrange(
  B_gammaerror, B_serror, B_betaerror,
  ncol = 3, nrow = 1,
  common.legend = F,
  #legend = "right",
  font.label = list(size = 24, color = "black", face = "bold")
)
#ggsave("/nas/longleaf/home/adaigle/DFESelfing/figures_for_publication/revisionfigure2.svg", plot = error_plot, width = 8.5, height = 4.5, dpi = 300)
saveRDS(error_plot_B, file = "/nas/longleaf/home/adaigle/DFESelfing/scripts/pylibseq/lowrec/selfing_error_plot_B.RDS")
write.table(summarize_df_selected, file = '/nas/longleaf/home/adaigle/DFESelfing/scripts/pylibseq/summarize_df_selfing.txt', quote=F, row.names=F)
