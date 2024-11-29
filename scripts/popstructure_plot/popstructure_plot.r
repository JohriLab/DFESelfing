library(ggplot2)
library(ggpubr)

base_dir <- "/nas/longleaf/home/adaigle/DFESelfing/"
figures_dir <- paste0(base_dir, "figures_for_publication/")

Nem2 <-  read.csv(paste0(base_dir, "scripts/popstructure_plot/Nem2.csv"))
Nem1 <-  read.csv(paste0(base_dir, "scripts/popstructure_plot/Nem1.csv"))
Nem01 <- read.csv(paste0(base_dir, "scripts/popstructure_plot/Nem01.csv"))
Nem05 <- read.csv(paste0(base_dir, "scripts/popstructure_plot/Nem05.csv"))

Nem01_even <- Nem01 %>% mutate(sampling="Even sampling") %>% filter(DFE=="DFE1")
Nem01_uneven_1 <- read.csv(paste0(base_dir, "scripts/popstructure_uneven_1_plot/Nem01.csv")) %>%
  mutate(sampling="Primarily one deme") %>% filter(DFE=="DFE1")
Nem01_uneven_2 <- read.csv(paste0(base_dir, "scripts/popstructure_uneven_2_plot/Nem01.csv")) %>% 
  mutate(sampling="Primarily two demes") %>% filter(DFE=="DFE1")

Nem01s <- rbind(Nem01_even,Nem01_uneven_1,Nem01_uneven_2)

Nem01_allsamplings_plt <- ggplot(Nem01s, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", "Adjusted DFE", "DFE-alpha", "GRAPES")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(sampling), cols = vars(selfing_class)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040",rep(c("grey", "#F8766D", "purple"),3))) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), 
    axis.title.x=element_text(size=15),axis.title.y=element_text(size=15), strip.text = element_text(size=15),
    plot.title= element_text(size=0), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(expression(italic(f[0])), expression(italic(f[1])), expression(italic(f[2])), expression(italic(f[3]))))


Nem01_plt <- ggplot(Nem01, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", "Adjusted DFE", "DFE-alpha", "GRAPES")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values = c("#404040",rep(c("grey", "#F8766D", "purple"),3))) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), 
    axis.title.x=element_text(size=0),axis.title.y=element_text(size=15), strip.text = element_text(size=12),
    plot.title= element_text(size=0), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(expression(italic(f[0])), expression(italic(f[1])), expression(italic(f[2])), expression(italic(f[3]))))

Nem1_plt <- ggplot(Nem1, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", "Adjusted DFE", "DFE-alpha", "GRAPES")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  scale_fill_manual(values = c("#404040",rep(c("grey", "#F8766D", "purple"),3))) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), 
    axis.title.x=element_text(size=0),axis.title.y=element_text(size=15), strip.text = element_text(size=12),
    plot.title= element_text(size=0), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(expression(italic(f[0])), expression(italic(f[1])), expression(italic(f[2])), expression(italic(f[3]))))

Nem2_plt <- ggplot(Nem2, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", "Adjusted DFE", "DFE-alpha", "GRAPES")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040",rep(c("grey", "#F8766D", "purple"),3))) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), 
    axis.title.x=element_text(size=0),axis.title.y=element_text(size=15), strip.text = element_text(size=12),
    plot.title= element_text(size=0), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(expression(italic(f[0])), expression(italic(f[1])), expression(italic(f[2])), expression(italic(f[3]))))

Nem05_plt <- ggplot(Nem05, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("Simulated DFE", "Adjusted DFE", "DFE-alpha", "GRAPES")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040",rep(c("grey", "#F8766D", "purple"),3))) + 
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), 
    axis.title.x=element_text(size=0),axis.title.y=element_text(size=15), strip.text = element_text(size=12),
    plot.title= element_text(size=0), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(expression(italic(f[0])), expression(italic(f[1])), expression(italic(f[2])), expression(italic(f[3]))))

sfigure13 <- ggarrange(Nem2_plt, Nem1_plt, Nem05_plt, Nem01_plt,
                    labels = c("A", "B", "C", "D"),
                    font.label = list(size = 22, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 2,
                    common.legend = T, legend = "bottom")


ggsave(paste0(figures_dir, "figure8.svg"), plot = Nem01_allsamplings_plt, width = 8.5, height = 8.5, dpi = 300)
ggsave(paste0(figures_dir, "new_figure7.jpg"), plot = Nem01_allsamplings_plt, width = 8.5, height = 8.5, dpi = 300)
ggsave(paste0(figures_dir, "new_figure7.tiff"), plot = Nem01_allsamplings_plt, width = 8.5, height = 8.5, dpi = 300)

ggsave(paste0(figures_dir, "sfigure13.svg"), plot = sfigure13, width = 8.5, height = 9, dpi = 150)
