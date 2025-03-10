# Comparing Kconceptus2pu_initial vs Kconceptus2pu_final for all 103 chemicals
# updated: 3/10/25

rm(list=ls())
library(httk)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(dplyr)
library(tidyverse)
library(latex2exp)

# load data from prioritization for testing 
load('./data/invitrodb_v3_5_deiod_filtered_httk.RData', 
     verbose = TRUE)

load_dawson2021() # most data for ToxCast chems
load_sipes2017() # most data for pharma compounds
load_pradeep2020() # best ML method for in silico prediction

parameters <- c("Kconceptus2pu_initial", "Kconceptus2pu_final")

# function to return specific parameter used in 1tri_pbtk
get_param <- function(dtxsid, parameter) {
    params <- parameterize_1tri_pbtk(dtxsid = dtxsid, physchem.exclude = FALSE)
    return(params[[parameter]])
}

df <- matrix(0, ncol = length(parameters), nrow = nrow(ivive.moe.tb))

for(i in 1:nrow(ivive.moe.tb)) {
  df[i, ] = sapply(parameters, get_param, dtxsid = ivive.moe.tb$dtxsid[i])
  cat("Done for", i, "chemicals \n")
}

colnames(df) <- parameters
df <- data.frame(df)
df$dtxsid <- ivive.moe.tb$dtxsid
df$chnm <- ivive.moe.tb$chnm
df <- df[, c("dtxsid", "chnm", parameters)]

# log10 transform the data first 
df[c("Kconceptus2pu_initial", "Kconceptus2pu_final")] <- lapply(df[c("Kconceptus2pu_initial", "Kconceptus2pu_final")], log10)
df$m <- (df$Kconceptus2pu_final - df$Kconceptus2pu_initial)/13
setDT(df)

my_theme <-  theme_bw() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) 

endpoints.gg <- ggplot(df, aes(x = Kconceptus2pu_initial, 
               y = Kconceptus2pu_final)) +
  geom_point(alpha = 0.75) + 
  geom_abline(slope = 1, linewidth = 1, linetype = 'solid') +
  geom_abline(aes(slope = 1, intercept = 1), linetype = 'longdash', color = 'red', linewidth = 1) + # default linesize = 1
  geom_abline(aes(slope = 1, intercept = -1), linetype = 'longdash', color = 'red', linewidth = 1) +
  geom_label_repel(aes(label = chnm),
                   data = df[Kconceptus2pu_initial < -0.1 | Kconceptus2pu_initial > 5],
                   size = 3.5, 
                   force = 500,
                   max.overlaps = 10) +
  labs(x = TeX("$log_{10} K_{conceptus}^{initial}$ = average of maternal tissue partition coeffs"), 
       y = TeX("$log_{10} K_{conceptus}^{final} = log_{10} K_{placenta2pu}$")) +
  my_theme

# plot Kconceptus2pu linear relationships for all chems
ylims <- c(floor(min(df$Kconceptus2pu_initial)), ceiling(max(df$Kconceptus2pu_final)))
linear.plt <- ggplot(df) + 
  geom_abline(aes(intercept = Kconceptus2pu_initial, slope = m, color = factor(chnm)), 
              linewidth = 1.5) +
  scale_x_continuous(breaks = seq(0,13), 
                     labels = as.character(seq(0,13)), 
                     limits = c(0,13)) +
  scale_y_continuous(breaks = seq(ylims[1], ylims[2], 1), 
                     limits = ylims) +
  labs(x = "Gestational Age (weeks)", y = TeX("$log_{10} K_{conceptus2pu}$")) +
  my_theme +
  theme(legend.position = "none")

# combine into one figure 
fig <- plot_grid(endpoints.gg, linear.plt,
          labels = c("A", "B"),
          label_size = 25,
          hjust = 0, align = "h", 
          nrow = 1, ncol = 2)

ggsave(plot = fig, 
       units = "in", 
       dpi = 300, 
       width = 12.5, height = 5.55, 
       device = "jpg", 
       filename = "./figures/supp/Kconceptus2pu-v3.jpg")


