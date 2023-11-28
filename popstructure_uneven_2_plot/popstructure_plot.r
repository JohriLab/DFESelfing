library(ggplot2)
library(igraph)
library(ggraph)
library(ggpubr)
# Create an empty graph
island_graph <- make_empty_graph(n = 5)

# Add vertices in a circle
#island_graph <- add_vertices(island_graph, 5)
layout <- layout_in_circle(island_graph)
for (i in 1:5) {
  for (j in 1:5) {
    if (i != j) {
      island_graph <- add_edges(island_graph, c(i, j))
    }
  }
}

# Create a vector of labels
#labels <- rep("m_d=m/5", ecount(island_graph))

# Plot the graph with labels
#ggarrange(islands, Nem2_plt, Nem1_plt, Nem01_plt,
#                        labels = c("A", "B", "C", "D"),
#                        font.label = list(size = 24, color = "black", face = "bold", family = NULL),
#                        ncol = 2, nrow = 2,
#                        common.legend = T, legend = "bottom")
 
    # Plot the graph with labels, evenly spaced nodes, and larger node points
#islands <- ggraph(island_graph, layout = circle_layout) +
#      geom_edge_arc(aes(x = x, y = y), strength = 0, arrow = arrow(type = "closed", length = unit(0.15, "inches")), end_cap = circle(3, 'mm')) +
#      geom_node_point(size = 5) +
#      theme_void()

Nem2 <- read.csv("/nas/longleaf/home/adaigle/DFESelfing/popstructure_uneven_2_plot/Nem2.csv")
Nem1 <- read.csv("/nas/longleaf/home/adaigle/DFESelfing/popstructure_uneven_2_plot/Nem10.csv")
Nem01 <- read.csv("/nas/longleaf/home/adaigle/DFESelfing/popstructure_uneven_2_plot/Nem01.csv")
Nem05 <- read.csv("/nas/longleaf/home/adaigle/DFESelfing/popstructure_uneven_2_plot/Nem05.csv")

Nem01_plt <- ggplot(Nem01, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("truth", "DFE-alpha", "GRAPES","F_adjusted_0", "true0", 0, "0_grapes",
        "F_adjusted_50", "true50", 50, "50_grapes", "F_adjusted_80", "true80", 80, "80_grapes",
        "F_adjusted_90", "true90", 90, "90_grapes", "F_adjusted_95", "true95", 95, "95_grapes",
        "F_adjusted_99", "true99", 99, "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040",rep(c("#F8766D", "purple"),3))) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=15),axis.title.y=element_text(size=20), strip.text=element_text(size=12),
  plot.title= element_text(size=20), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(~f[0], ~f[1], ~f[2], ~f[3]))

Nem1_plt <- ggplot(Nem1, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("truth", "DFE-alpha", "GRAPES","F_adjusted_0", "true0", 0, "0_grapes",
        "F_adjusted_50", "true50", 50, "50_grapes", "F_adjusted_80", "true80", 80, "80_grapes",
        "F_adjusted_90", "true90", 90, "90_grapes", "F_adjusted_95", "true95", 95, "95_grapes",
        "F_adjusted_99", "true99", 99, "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040",rep(c("#F8766D", "purple"),3))) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=15),axis.title.y=element_text(size=20), strip.text=element_text(size=12),
  plot.title= element_text(size=20), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(~f[0], ~f[1], ~f[2], ~f[3]))

Nem2_plt <- ggplot(Nem2, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("truth", "DFE-alpha", "GRAPES","F_adjusted_0", "true0", 0, "0_grapes",
        "F_adjusted_50", "true50", 50, "50_grapes", "F_adjusted_80", "true80", 80, "80_grapes",
        "F_adjusted_90", "true90", 90, "90_grapes", "F_adjusted_95", "true95", 95, "95_grapes",
        "F_adjusted_99", "true99", 99, "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040",rep(c("#F8766D", "purple"),3))) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=15),axis.title.y=element_text(size=20), strip.text=element_text(size=12),
  plot.title= element_text(size=20), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(~f[0], ~f[1], ~f[2], ~f[3]))

Nem05_plt <- ggplot(Nem05, aes(x = generation, y = value, fill = factor(selfing, 
    levels = c("truth", "DFE-alpha", "GRAPES","F_adjusted_0", "true0", 0, "0_grapes",
        "F_adjusted_50", "true50", 50, "50_grapes", "F_adjusted_80", "true80", 80, "80_grapes",
        "F_adjusted_90", "true90", 90, "90_grapes", "F_adjusted_95", "true95", 95, "95_grapes",
        "F_adjusted_99", "true99", 99, "99_grapes")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  expand_limits(y=c(0,1)) +
  facet_grid(rows = vars(DFE), cols = vars(selfing_class)) +
  #scale_fill_manual(values = c("#404040", rep(c("#00BA38", "#619CFF", "#F8766D", "purple"),6))) + 
  scale_fill_manual(values = c("#404040",rep(c("#F8766D", "purple"),3))) + 
  theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15),
  axis.title.x=element_text(size=15),axis.title.y=element_text(size=20), strip.text=element_text(size=12),
  plot.title= element_text(size=20), legend.position = "bottom", legend.text = element_text(size=12)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
  scale_x_discrete(labels = c(~f[0], ~f[1], ~f[2], ~f[3]))

sfigure15 <- ggarrange(Nem2_plt, Nem1_plt, Nem05_plt, Nem01_plt,
                    labels = c("A", "B", "C", "D"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 2,
                    common.legend = T, legend = "bottom")


#ggsave("/nas/longleaf/home/adaigle/DFESelfing/figures_for_publication/figure8.svg", plot = Nem01_plt, width = 8.5, height = 8.5, dpi = 600)
ggsave("/nas/longleaf/home/adaigle/DFESelfing/figures_for_publication/sfigure15.svg", plot = sfigure15, width = 8.5, height = 10, dpi = 600)
