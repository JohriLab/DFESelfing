rm(list=ls())
library(tidyverse)
library(ggpubr)

base_dir <- "/nas/longleaf/home/adaigle/DFESelfing/"
sim_outputs_dir <- "/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/"
figures_dir <- paste0(base_dir, "figures_for_publication/")
dominance_dir <- paste0(base_dir, "scripts/dominance/")
h01 <- read.csv(paste0(dominance_dir, "h01.csv"))
h05 <- read.csv(paste0(dominance_dir, "h05.csv"))
h025 <- read.csv(paste0(dominance_dir, "h025.csv"))
h075 <- read.csv(paste0(dominance_dir, "h075.csv"))

combined_df <- rbind(h01, h05, h025, h075) %>%
    mutate(generation = recode(generation,
     'f0' = '0<2Nsh<1',
     'f1' = '1<2Nsh<10', 
     'f2' = '10<2Nsh<100', 
     'f3' = '100<2Nsh<inf')) 

plot_99 <- combined_df %>% filter(selfing_class=="99% Selfing") %>%
    mutate(selfing = recode(selfing,
     'truth' = 'Simulated DFE',
     '99' = 'DFEalpha', 
     '99_grapes' = 'GRAPES'))
plt99 <- ggplot(plot_99, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", "Dominance_adjusted_99", "Dominance_adjusted_50","F_adjusted_0", "true0", 0, "0_grapes",
        "F_adjusted_50", "true50", 50, "50_grapes", "F_adjusted_80", "true80", 80, "80_grapes",
        "F_adjusted_90", "true90", 90, "90_grapes", "F_adjusted_95", "true95", 95, "95_grapes",
        "F_adjusted_99", "true99", 'DFEalpha', "GRAPES")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "99% Selfing", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(dominance)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(legend.position="right", axis.text.x=element_text(size=12, angle = 45, vjust=0.5), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), plot.title= element_text(size=25))


plot_50 <- combined_df %>% filter(selfing_class=="50% Selfing") %>%
    mutate(selfing = recode(selfing,
     'Dominance_adjusted_50' = 'Simulated DFE',
     '50' = 'DFEalpha', 
     '50_grapes' = 'GRAPES'))
plt50 <-ggplot(plot_50, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", "F_adjusted_0", "true0", 0, "0_grapes",
        "F_adjusted_50", "true50", 'DFEalpha', "GRAPES", "F_adjusted_80", "true80", 80, "80_grapes",
        "F_adjusted_90", "true90", 90, "90_grapes", "F_adjusted_95", "true95", 95, "95_grapes",
        "F_adjusted_99", "true99", 99, "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "50% Selfing", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(dominance)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(legend.position="right", axis.text.x=element_text(size=12, angle = 45, vjust=0.5), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), plot.title= element_text(size=25))


plot_0 <- combined_df %>% filter(selfing_class=="0% Selfing")%>%
    mutate(selfing = recode(selfing,
     'Dominance_adjusted_0' = 'Simulated DFE', #changed this name so the legend reads correctly after ggarrange
     '0' = 'DFEalpha', 
     '0_grapes' = 'GRAPES'))
plt0 <-ggplot(plot_0, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c('Simulated DFE', 'Dominance_adjusted_0', "F_adjusted_0", "true0", 'DFEalpha', "GRAPES",
        "F_adjusted_50", "true50", 50, "50_grapes", "F_adjusted_80", "true80", 80, "80_grapes",
        "F_adjusted_90", "true90", 90, "90_grapes", "F_adjusted_95", "true95", 95, "95_grapes",
        "F_adjusted_99", "true99", 99, "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "0% Selfing", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(dominance)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(legend.position="right", axis.text.x=element_text(size=12, angle = 45, vjust=0.5), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), strip.text = element_text(size=15), plot.title= element_text(size=25))

ggarrange(plt0, plt99,
                    labels = c("A", "B", "C"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "right", vjust=1)



plot_DFE1 <- combined_df %>% filter(DFE=="DFE1") %>% filter(selfing_class!="50% Selfing") %>%
    mutate(selfing = recode(selfing,
     'Dominance_adjusted_99' = 'Simulated DFE',
     '99' = 'DFEalpha', 
     '99_grapes' = 'GRAPES',
     'Dominance_adjusted_0' = 'Simulated DFE', #changed this name so the legend reads correctly after ggarrange
     '0' = 'DFEalpha', 
     '0_grapes' = 'GRAPES'))
pltDFE1 <- ggplot(plot_DFE1, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", 'DFEalpha', "GRAPES")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(selfing_class), cols = vars(dominance)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040", rep(c("#F8766D", "purple"),6))) + 
    theme(axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),
  axis.title.x=element_text(size=20),axis.title.y=element_text(size=20), strip.text = element_text(size=15),
  plot.title= element_text(size=20), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))+
  scale_x_discrete(labels = c(~f[0], ~f[1], ~f[2], ~f[3]))
pltDFE1


plot_DFE2 <- combined_df %>% filter(DFE=="DFE2") %>% filter(selfing_class!="50% Selfing") %>%
    mutate(selfing = recode(selfing,
     'Dominance_adjusted_99' = 'Simulated DFE',
     '99' = 'DFEalpha', 
     '99_grapes' = 'GRAPES',
     'Dominance_adjusted_0' = 'Simulated DFE', #changed this name so the legend reads correctly after ggarrange
     '0' = 'DFEalpha', 
     '0_grapes' = 'GRAPES'))
pltDFE2 <- ggplot(plot_DFE2, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", 'DFEalpha', "GRAPES")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(selfing_class), cols = vars(dominance)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040", rep(c("#F8766D", "purple"),6))) + 
    theme(axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),
  axis.title.x=element_text(size=20),axis.title.y=element_text(size=20), strip.text = element_text(size=15),
  plot.title= element_text(size=20), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))+
  scale_x_discrete(labels = c(~f[0], ~f[1], ~f[2], ~f[3]))
pltDFE2


plot_DFE3 <- combined_df %>% filter(DFE=="DFE3") %>% filter(selfing_class!="50% Selfing") %>%
    mutate(selfing = recode(selfing,
     'Dominance_adjusted_99' = 'Simulated DFE',
     '99' = 'DFEalpha', 
     '99_grapes' = 'GRAPES',
     'Dominance_adjusted_0' = 'Simulated DFE', #changed this name so the legend reads correctly after ggarrange
     '0' = 'DFEalpha', 
     '0_grapes' = 'GRAPES'))
pltDFE3 <- ggplot(plot_DFE3, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", 'DFEalpha', "GRAPES")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(selfing_class), cols = vars(dominance)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040", rep(c("#F8766D", "purple"),6))) + 
    theme(axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=20),axis.text.x=element_text(size=15), axis.title.y=element_text(size=20), strip.text = element_text(size=15),
  plot.title= element_text(size=20), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))+
  scale_x_discrete(labels = c(~f[0], ~f[1], ~f[2], ~f[3]))
  
pltDFE3


figure6_old <- ggarrange(pltDFE1, pltDFE2, pltDFE3,
                    labels = c("A", "B", "C"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 1, nrow = 3,
                    common.legend = TRUE, legend = "bottom", vjust=1)
#ggsave("/nas/longleaf/home/adaigle/DFESelfing/figures_for_publication/figure6.svg", plot = figure6, width = 8.5, height = 12, dpi = 600)

plot_0 <- combined_df %>% filter(selfing_class=="0% Selfing") %>%
    mutate(selfing = recode(selfing,
     'Dominance_adjusted_99' = 'Simulated DFE',
     '99' = 'DFEalpha', 
     '99_grapes' = 'GRAPES',
     'Dominance_adjusted_0' = 'Simulated DFE', #changed this name so the legend reads correctly after ggarrange
     '0' = 'DFEalpha', 
     '0_grapes' = 'GRAPES')) 
plt0 <- ggplot(plot_0, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", 'DFEalpha', "GRAPES")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(dominance)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), 
    axis.title.x=element_text(size=15),axis.title.y=element_text(size=15), strip.text = element_text(size=15),
    plot.title= element_text(size=0), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(expression(italic(f[0])), expression(italic(f[1])), expression(italic(f[2])), expression(italic(f[3]))))

plot_99 <- combined_df %>% filter(selfing_class=="99% Selfing") %>%
    mutate(selfing = recode(selfing,
     'Dominance_adjusted_99' = 'Simulated DFE',
     '99' = 'DFEalpha', 
     '99_grapes' = 'GRAPES',
     'Dominance_adjusted_0' = 'Simulated DFE', #changed this name so the legend reads correctly after ggarrange
     '0' = 'DFEalpha', 
     '0_grapes' = 'GRAPES')) 
plt99 <- ggplot(plot_99, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", 'DFEalpha', "GRAPES")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(dominance)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), 
    axis.title.x=element_text(size=15),axis.title.y=element_text(size=15), strip.text = element_text(size=15),
    plot.title= element_text(size=0), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(expression(italic(f[0])), expression(italic(f[1])), expression(italic(f[2])), expression(italic(f[3]))))


figure6 <- ggarrange(plt0, plt99,
                    labels = c("A", "B"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "bottom", vjust=1)
ggsave(paste0(figures_dir, "figure6.svg"), plot = figure6, width = 8.5, height = 10, dpi = 300)

plot_50 <- combined_df %>% filter(selfing_class=="50% Selfing") %>%
    mutate(selfing = recode(selfing,
     'Dominance_adjusted_99' = 'Simulated DFE',
     '99' = 'DFEalpha', 
     '99_grapes' = 'GRAPES',
     'Dominance_adjusted_50' = 'Simulated DFE', #changed this name so the legend reads correctly after ggarrange
     '50' = 'DFEalpha', 
     '50_grapes' = 'GRAPES')) 
sfigure07 <- ggplot(plot_50, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", 'DFEalpha', "GRAPES")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(dominance)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), 
    axis.title.x=element_text(size=15),axis.title.y=element_text(size=15), strip.text = element_text(size=15),
    plot.title= element_text(size=0), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(expression(italic(f[0])), expression(italic(f[1])), expression(italic(f[2])), expression(italic(f[3]))))

ggsave(paste0(figures_dir, "sfigure07.svg"), plot = sfigure07, width = 8.5, height = 8.5, dpi = 150)
