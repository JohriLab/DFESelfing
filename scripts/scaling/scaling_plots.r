library(tidyverse)
library(ggpubr)

q10_plot <- read.csv("/nas/longleaf/home/adaigle/DFESelfing/scripts/scaling/q10.csv")
q10_Btable <- read.csv("/nas/longleaf/home/adaigle/DFESelfing/scripts/scaling/q10_B.csv")
q20_plot <- read.csv("/nas/longleaf/home/adaigle/DFESelfing/scripts/scaling/q20.csv")
q20_Btable <- read.csv("/nas/longleaf/home/adaigle/DFESelfing/scripts/scaling/q20_B.csv")
q50_plot <- read.csv("/nas/longleaf/home/adaigle/DFESelfing/scripts/scaling/q50.csv")
q50_Btable <- read.csv("/nas/longleaf/home/adaigle/DFESelfing/scripts/scaling/q50_B.csv")
q100_plot <- read.csv("/nas/longleaf/home/adaigle/DFESelfing/scripts/scaling/q100.csv")
q100_Btable <- read.csv("/nas/longleaf/home/adaigle/DFESelfing/scripts/scaling/q100_B.csv")

plot <- rbind(q10_plot, q20_plot, q50_plot, q100_plot)
Btable <- rbind(q10_Btable, q20_Btable, q50_Btable, q100_Btable)

plot_99 <- plot %>% filter(selfing_class=="99% Selfing")
Btable_99 <- Btable %>% filter(selfing_class=="99% Selfing")

plt99 <- ggplot(plot_99, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", "Adjusted DFE", "DFE-alpha", "GRAPES")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(scale)) +
  scale_fill_manual(values = c("#404040", rep(c("grey", "#F8766D", "purple"),6))) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), 
    axis.title.x=element_text(size=12),axis.title.y=element_text(size=12), strip.text = element_text(size=13),
    plot.title= element_text(size=0), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(expression(italic(f[0])), expression(italic(f[1])), expression(italic(f[2])), expression(italic(f[3])))) +
  geom_text(data = Btable_99, 
    aes(label = paste("B = ", format(round(B_avg, 2), nsmall = 2, digits = 2))), 
    vjust = -0.5, hjust=0.25, size = 4, position = position_dodge(width = 0.9), fontface="italic") 

plot_95 <- plot %>% filter(selfing_class=="95% Selfing") %>%
    filter(!(selfing=="Simulated DFE" & scale=="10" & DFE=="DFE2"))
Btable_95 <- Btable %>% filter(selfing_class=="95% Selfing")

plt95 <- ggplot(plot_95, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", "Adjusted DFE", "DFE-alpha", "GRAPES")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(scale)) +
  scale_fill_manual(values = c("#404040", rep(c("grey", "#F8766D", "purple"),6))) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), 
    axis.title.x=element_text(size=12),axis.title.y=element_text(size=12), strip.text = element_text(size=13),
    plot.title= element_text(size=0), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(expression(italic(f[0])), expression(italic(f[1])), expression(italic(f[2])), expression(italic(f[3])))) +
  geom_text(data = Btable_95, 
    aes(label = paste("B = ", format(round(B_avg, 2), nsmall = 2, digits = 2))), 
    vjust = -0.5, hjust=0.25, size = 4, position = position_dodge(width = 0.9), fontface="italic") 

figure <- ggarrange(plt95, plt99,
                    labels = c("A", "B"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "bottom", vjust=1)

ggsave(paste0(figures_dir, "sfigure7_new.svg"), plot = figure, width = 8.5, height = 10, dpi = 150)
