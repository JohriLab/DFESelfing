rm(list=ls())
library(RColorBrewer)
library(tidyverse)
library(reshape2)
library(ggpubr)

base_dir <- "/nas/longleaf/home/adaigle/DFESelfing/"
figures_dir <- paste0(base_dir, "figures_for_publication/")

# Function to create a ggarrange figure
create_figure <- function(DFE_number) {
  
  # Define paths based on DFE_number
  s_dfealpha_path <- paste0("/nas/longleaf/home/adaigle/DFESelfing/scripts/summary_figures/DFE", DFE_number, "_s_dfealpha.rds")
  beta_path <- paste0("/nas/longleaf/home/adaigle/DFESelfing/scripts/summary_figures/DFE", DFE_number, "_beta.rds")
  estimate_path <- paste0("/nas/longleaf/home/adaigle/DFESelfing/scripts/summary_figures/DFE", DFE_number, "_estimate.rds")
  estimate_s_path <- paste0("/nas/longleaf/home/adaigle/DFESelfing/scripts/summary_figures/DFE", DFE_number, "_estimate_s.rds")
  
  # Load the plots
  s_dfealpha <- readRDS(s_dfealpha_path) + theme(legend.position = "none") +
    labs(y = expression(italic(s[d])))
  beta <- readRDS(beta_path) + theme(legend.position = "none") +
    labs(y = expression(beta))
  estimate <- readRDS(estimate_path) + theme(axis.title.y = element_text(size=10), strip.text.y = element_blank())
  estimate_s <- readRDS(estimate_s_path) + theme(axis.title.y = element_text(size=10), axis.text.x=element_text(size=10), legend.position = "none") + 
    theme(strip.text.y = element_blank()) +
    labs(x = expression(plain("selection coefficient ") * " " * "(" * italic(s[d])* ")"))

  legend <- get_legend(estimate)
  # Remove the legend from the estimate plot
  estimate <- estimate + theme(legend.position = "none")
  # Arrange the plots into a single figure
  figure <- ggarrange(
    s_dfealpha, beta, estimate, estimate_s,
    labels = c("A", "B", "C", "D"),
    ncol = 2, nrow = 2, common.legend = FALSE
  )
  # Add the extracted legend below the figure
  figure_with_legend <- ggarrange(figure, legend, ncol = 1, heights = c(1, 0.1))
  
  
  return(figure_with_legend)
}

create_figure_lowrec <- function(DFE_number) {
  
  # Define paths based on DFE_number
  s_dfealpha_path <- paste0("/nas/longleaf/home/adaigle/DFESelfing/scripts/summary_figures/DFE", DFE_number, "_s_dfealpha_lowrec.rds")
  beta_path <- paste0("/nas/longleaf/home/adaigle/DFESelfing/scripts/summary_figures/DFE", DFE_number, "_beta_lowrec.rds")
  estimate_path <- paste0("/nas/longleaf/home/adaigle/DFESelfing/scripts/summary_figures/DFE", DFE_number, "_estimate_lowrec.rds")
  estimate_s_path <- paste0("/nas/longleaf/home/adaigle/DFESelfing/scripts/summary_figures/DFE", DFE_number, "_estimate_s_lowrec.rds")
  
  # Load the plots
  s_dfealpha <- readRDS(s_dfealpha_path) + theme(legend.position = "none") +
    labs(y = expression(italic(s[d])))
  beta <- readRDS(beta_path) + theme(legend.position = "none") +
    labs(y = expression(beta))
  estimate <- readRDS(estimate_path) + theme(axis.title.y = element_text(size=10), strip.text.y = element_blank())
  estimate_s <- readRDS(estimate_s_path) + theme(axis.title.y = element_text(size=10), axis.text.x=element_text(size=10), legend.position = "none") + 
    theme(strip.text.y = element_blank()) +
    labs(x = expression(plain("selection coefficient ") * " " * "(" * italic(s[d])* ")"))

  legend <- get_legend(estimate)
  # Remove the legend from the estimate plot
  estimate <- estimate + theme(legend.position = "none")
  # Arrange the plots into a single figure
  figure <- ggarrange(
    s_dfealpha, beta, estimate, estimate_s,
    labels = c("A", "B", "C", "D"),
    ncol = 2, nrow = 2, common.legend = FALSE
  )
  # Add the extracted legend below the figure
  figure_with_legend <- ggarrange(figure, legend, ncol = 1, heights = c(1, 0.1))
  
  
  return(figure_with_legend)
}

f1 <- create_figure(1)
f2 <- create_figure(2)
f3 <- create_figure(3)

f4 <- create_figure_lowrec(1)
f5 <- create_figure_lowrec(2)
f6 <- create_figure_lowrec(3)

figure <- ggarrange(
  f1, f4,
  #labels = c("A", "B", "C", "D"),
  ncol = 1, nrow = 2, common.legend = T
)

ggsave(paste0(figures_dir, "new_figure_01.svg"), plot = f1, width = 8.5, height = 8.5, dpi = 150)
ggsave(paste0(figures_dir, "new_figure_02.svg"), plot = f4, width = 8.5, height = 8.5, dpi = 150)

#ggsave(paste0(figures_dir, "new_figure_02.png"), plot = f4, width = 10, height = 8.5, dpi = 150)
#ggsave(paste0(figures_dir, "sfigure01.svg"), plot = sfigure01, width = 8.5, height = 9.5, dpi = 150)
#ggsave(paste0(figures_dir, "sfigure01.svg"), plot = sfigure01, width = 8.5, height = 9.5, dpi = 150)
#ggsave(paste0(figures_dir, "sfigure01.svg"), plot = sfigure01, width = 8.5, height = 9.5, dpi = 150)
#ggsave(paste0(figures_dir, "sfigure01.svg"), plot = sfigure01, width = 8.5, height = 9.5, dpi = 150)
