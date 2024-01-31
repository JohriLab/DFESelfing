library(tidyverse)
base_dir <- "/nas/longleaf/home/adaigle/DFESelfing/"
figures_dir <- paste0(base_dir, "figures_for_publication/")

int500 <- read.csv(paste0(base_dir, "scripts/intergenic_plot/int500.csv"))
int1500 <- read.csv(paste0(base_dir, "scripts/intergenic_plot/int1500.csv"))
int3000 <- read.csv(paste0(base_dir, "scripts/intergenic_plot/int3000.csv"))

plt_data <- rbind(int500, int1500, int3000) %>% 
    mutate(selfing = recode(selfing,
     'truth' = "Simulated DFE",
     '0' = 'DFE-alpha', 
     '0_grapes' = 'GRAPES',
      '50' = 'DFE-alpha', 
     '50_grapes' = 'GRAPES',
     '80' = 'DFE-alpha', 
     '80_grapes' = 'GRAPES',
     '90' = 'DFE-alpha', 
     '90_grapes' = 'GRAPES',
     '95' = 'DFE-alpha', 
     '95_grapes' = 'GRAPES',
     '99' = 'DFE-alpha', 
     '99_grapes' = 'GRAPES'))

plot_99 <- plt_data %>% filter(selfing_class=="99% Selfing")
ggplot_99<- ggplot(plot_99, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", "DFE-alpha", "GRAPES", 0, "0_grapes",
        "F_adjusted_50", "true50", 50, "50_grapes", "F_adjusted_80", "true80", 80, "80_grapes",
        "F_adjusted_90", "true90", 90, "90_grapes", "F_adjusted_95", "true95", 95, "95_grapes",
        "F_adjusted_99", "true99", 99, "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "DFEalpha and Grapes, 99% Selfing, Varying Coding Density", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(intergenic)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(legend.position="none", axis.text.x=element_text(size=15), axis.text.y=element_text(size=12),legend.title = element_blank(),
  axis.title.x=element_text(size=20),axis.title.y=element_text(size=15), strip.text = element_text(size=15), plot.title= element_blank())+
  scale_x_discrete(labels = c(~f[0], ~f[1], ~f[2], ~f[3]))


plot_50 <- plt_data %>% filter(selfing_class=="50% Selfing")
ggplot_50 <- ggplot(plot_50, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", "DFE-alpha", "GRAPES", 0, "0_grapes",
        "F_adjusted_50", "true50", 50, "50_grapes", "F_adjusted_80", "true80", 80, "80_grapes",
        "F_adjusted_90", "true90", 90, "90_grapes", "F_adjusted_95", "true95", 95, "95_grapes",
        "F_adjusted_99", "true99", 99, "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "DFEalpha and Grapes, 50% Selfing, Varying Coding Density", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(intergenic)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(legend.position="none", axis.text.x=element_text(size=15), axis.text.y=element_text(size=12),legend.title = element_blank(),
  axis.title.x=element_text(size=20),axis.title.y=element_text(size=15), strip.text = element_text(size=15), plot.title= element_blank())+
  scale_x_discrete(labels = c(~f[0], ~f[1], ~f[2], ~f[3]))


plot_0 <- plt_data %>% filter(selfing_class=="0% Selfing")
ggplot_0<- ggplot(plot_0, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", "DFE-alpha", "GRAPES", 0, "0_grapes",
        "F_adjusted_50", "true50", 50, "50_grapes", "F_adjusted_80", "true80", 80, "80_grapes",
        "F_adjusted_90", "true90", 90, "90_grapes", "F_adjusted_95", "true95", 95, "95_grapes",
        "F_adjusted_99", "true99", 99, "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "DFEalpha and Grapes, 0% Selfing, Varying Coding Density", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(intergenic)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040", rep(c("#F8766D", "purple"),6))) + 
  theme(legend.position="none", axis.text.x=element_text(size=15), axis.text.y=element_text(size=12), legend.title = element_blank(),
  axis.title.x=element_text(size=20),axis.title.y=element_text(size=15), strip.text = element_text(size=15), plot.title= element_blank())+
  scale_x_discrete(labels = c(~f[0], ~f[1], ~f[2], ~f[3]))

sfigure05 <- ggarrange(ggplot_0, ggplot_50, ggplot_99,
                    labels = c("A", "B", "C"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 1, nrow = 3,
                    common.legend = TRUE, legend = "bottom")
ggsave(paste0(figures_dir, "sfigure05.svg"), plot = sfigure05, width = 8.5, height = 11, dpi = 600)
