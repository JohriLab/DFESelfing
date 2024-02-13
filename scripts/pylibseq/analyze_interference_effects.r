library(tidyverse)
library(ggplot2)
library(ggpubr)

#frequency filtered results
df <- read.table("/nas/longleaf/home/adaigle/DFESelfing/pylibseq/allresults_5000_0.stats", header=T)

#non frequency filtered results
#df <- read.table("/nas/longleaf/home/adaigle/DFESelfing/pylibseq/allresults_50000.stats", header=T)

gammabeta <- read.csv("/nas/longleaf/home/adaigle/DFESelfing/pylibseq/gammabeta.csv")
df <- df %>% mutate(
    selfing = as.numeric(str_extract(filename, "(?<=eqm_selfing)\\d+")),
    DFE = str_extract(filename, "(DFE)\\d+")
    ) %>% group_by(selfing, DFE)
gammabeta <- gammabeta %>% select(B, empirical_Ne, DFE, selfing, gamma, b, GammaZero.negGmean, GammaZero.negGshape, selfing_Ne, newNE) %>%
  group_by(selfing, DFE)

join_df <- left_join(df, gammabeta, multiple="all") 

join_df$thetadelta <- 1 - (join_df$thetapi / join_df$thetaw)

summarize_df <- join_df %>% ungroup() %>% group_by(filename,output,site,DFE,selfing) %>% 
    summarize(across(where(is.numeric), list(avg = ~mean(., na.rm = TRUE)), .names = "{col}")) %>% #.names makes sure columns aren't renamed to _avg yet
    ungroup() %>% group_by(DFE, selfing,site) %>%
    summarize(across(where(is.numeric), list(avg = ~mean(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE))))

lowfreqs <- read.csv("/nas/longleaf/home/adaigle/DFESelfing/pylibseq/freq1_5.csv")
summarize_df$Dprime_avg <- lowfreqs$Dprime_avg
summarize_df$Dprime_sd <- lowfreqs$Dprime_sd

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

quick_summary_all <- summarize_df %>% select(DFE,selfing,site,Dprime_avg,thetadelta_avg) %>% arrange(selfing)
write.csv(quick_summary_all, file = "/nas/longleaf/home/adaigle/DFESelfing/pylibseq/quicksumall.csv", quote=F)


plot_stat <- function(DFE, stat) {
  ggplot(DFE, aes(x = selfing, y = get(paste0(stat, "_avg")), fill=site)) +
  geom_errorbar(aes(ymin = get(paste0(stat, "_avg")) - get(paste0(stat, "_sd")), ymax = get(paste0(stat, "_avg")) + get(paste0(stat, "_sd")),color=site), 
    width = 5, position=position_dodge(0.1)) +
  geom_line(aes(color=site)) +
  geom_point() + 
  labs(x = "% Selfing",
       y = paste0(stat, "_avg")) +
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
  axis.title.x=element_text(size=12),axis.title.y=element_text(size=12), strip.text = element_text(size=15), 
  plot.title= element_text(size=0), legend.title = element_text(size=15), legend.text = element_text(size=15))
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
DFE3_Dprime <- plot_stat(DFE3, "Dprime") + labs(y = expression(italic("D'"))) + ylim(c(-0.975, -0.850))
DFE3_D <- plot_stat(DFE3, "D")
DFE3_tajimasd <- plot_stat(DFE3, "tajimasd")
DFE3_hapdiv <- plot_stat(DFE3, "hapdiv")
DFE3_inference <- plot_DFEinference(DFE3)
DFE3_brian <- plot_stat(DFE3, "thetadelta") + labs(y = expression(paste(Delta,theta))) + ylim(c(0, 0.75))
DFE3_rsq <- plot_stat(DFE3, "rsq")


ggarrange(DFE3_inference, DFE3_Dprime, DFE3_tajimasd, DFE3_hapdiv,
                    labels = c("A", "B", "C", "D"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 2,
                    common.legend = TRUE, legend = "right")

 plot_stat(DFE1, "hprime")
DFE1_Dprime <- plot_stat(DFE1, "Dprime") + labs(y = expression(italic("D'"))) + ylim(c(-0.975, -0.850))
DFE1_tajimasd <- plot_stat(DFE1, "tajimasd")
DFE1_hapdiv <- plot_stat(DFE1, "hapdiv")
DFE1_inference <- plot_DFEinference(DFE1)
DFE1_brian <- plot_stat(DFE1, "thetadelta") + labs(y = expression(italic(paste(Delta,theta)))) + ylim(c(0, 0.75))
DFE1_rsq <- plot_stat(DFE1, "rsq")


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
DFE2_brian <- plot_stat(DFE2, "thetadelta") + labs(y = expression(paste(Delta,theta))) + ylim(c(0, 0.75))
DFE2_rsq <- plot_stat(DFE2, "rsq")


ggarrange(DFE2_inference, DFE2_Dprime, DFE2_tajimasd, DFE2_hapdiv,
                    labels = c("A", "B", "C", "D"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 2,
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

figure4_D <- ggarrange(DFE1_brian, DFE1_D, DFE2_brian, DFE2_D, DFE3_brian, DFE3_D,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 3,
                    common.legend = TRUE, legend = "bottom", hjust=0)

DvsDprime <- ggarrange(DFE1_Dprime, DFE1_D, DFE2_Dprime, DFE2_D, DFE3_Dprime, DFE3_D,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 3,
                    common.legend = TRUE, legend = "bottom", hjust=0)
#annotate_figure(fig, 
#                left = "DFE3      DFE2     DFE1")

#ggsave("/nas/longleaf/home/adaigle/DFESelfing/figures_for_publication/figure9.pdf", plot = figure9, width = 8.5, height = 8.5, dpi = 600)
ggsave("/nas/longleaf/home/adaigle/DFESelfing/figures_for_publication/figure4.svg", plot = figure4, width = 8.5, height = 8.5, dpi = 300)
