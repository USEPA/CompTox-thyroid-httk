# ==============================================================================
# Reproduces Figure 6 of Truong et al. 2024. 
# 
# Demonstrates ingredients to stitch two HTTK models together. 
# Namely, plots Vconceptus, Wconceptus, and Kconceptus2pu as a function of time. 
# 
# @author: Kimberly Truong
# created: 10/24/23
# updated: 12/10/24
# ==============================================================================

rm(list=ls())
library(httk)
library(data.table)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(scales)
library(latex2exp)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(cowplot)

# Show Mass/Volume Conservation at Week 13 transition --------------------------
# estimate PK half-life of Fipronil
half.life <- calc_half_life(chem.name = "Fipronil")

# distinguish maternal from fetal compartments  
maternal_compts <- c('gut', 'liver', 'kidney', 'lung', 'ven', 'art', 
                     'adipose','thyroid', 'rest')
maternal_states <- paste0('A', maternal_compts)
maternal_states <- c(maternal_states, 
                     'Atubules', 'Ametabolized', 'AUC')

fetal_compts <- c(maternal_compts[!( maternal_compts %in% c('adipose') )], 
                  "brain")
fetal_states <- c(paste0('Af', fetal_compts), 'fAUC')

# dynamical outputs of each model 
firsttri.vols.out <- c("Vven", "Vart", "Vadipose", "Vrest", "Vffmx", "Vallx", "Vconceptus")
fetal.vols.out <- c("fBW", "Vamnf", "Vplacenta", firsttri.vols.out[firsttri.vols.out != "Vconceptus"], 
                    paste0("Vf", fetal_compts))

# we chose Fipronil as input but any chem would give the same output 
# day resolution is sufficient since t1/2 >> 1 day 
out <- solve_1tri_pbtk(chem.name = 'Fipronil', 
                       daily.dose = 1, 
                       doses.per.day = 1, 
                       times = seq(0, 13*7, 1),
                       monitor.vars = firsttri.vols.out)
out <- out[, c("time", firsttri.vols.out)] # remove the "Agutlumen" col 
out <- as.data.frame(out)

out2 <- solve_fetal_pbtk(chem.name = 'Fipronil', 
                         dose = 0, 
                         daily.dose = 1, 
                         doses.per.day = 1, 
                         times = seq(13*7, 40*7, 1), 
                         monitor.vars = fetal.vols.out)
out2 <- out2[, c("time", fetal.vols.out)]
out2 <- as.data.frame(out2)

# add in constant tissue volumes
# these are all maternal organs 
params <- parameterize_1tri_pbtk(chem.name = 'Fipronil')

add_tissuevol_constants <- function(df) {
  df$Vgut = params$Vgutc * params$pre_pregnant_BW / params$gut_density
  df$Vkidney = params$Vkidneyc * params$pre_pregnant_BW / params$kidney_density
  df$Vliver = params$Vliverc * params$pre_pregnant_BW / params$liver_density
  df$Vlung = params$Vlungc * params$pre_pregnant_BW / params$lung_density
  df$Vthyroid = params$Vthyroidc * params$pre_pregnant_BW / params$thyroid_density
  return(df)
}

full.out <- add_tissuevol_constants(out)
full.out2 <- add_tissuevol_constants(out2)
setDT(full.out)
setDT(full.out2)

# Show volume conservation at 13 weeks of gestational age 
ggdata <- rbindlist(list(full.out[, c("time", "Vconceptus")], 
                    full.out2[, c("time", "fBW", "Vplacenta", "Vamnf")]), 
                    fill = TRUE
                    )

# convert to weeks 
ggdata[, time := time / 7]
ggdata[, Vpoc := fBW + Vplacenta + Vamnf]
ggdata2 <- melt.data.table(ggdata[, c("time", "Vconceptus", "Vpoc")], 
                           id.vars = c("time"), 
                           variable.name = "volume")
ggdata2 <- ggdata2[!is.na(value)]

ggdata2[volume == "Vconceptus", model := "1tri_pbtk"]
ggdata2[volume == "Vpoc", model := "fetal_pbtk"]

# for reusing themes
my_theme <-  theme_bw() + 
  theme(title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 14, hjust = 0), 
        legend.background = element_blank()) 

volpp <- ggplot(ggdata2, 
       aes(x = time, y = value, color = model)) +
  geom_line(linewidth = 2) +
  scale_color_manual(labels = c(TeX("$V_{conceptus}$"), TeX("$V_{fetus} + V_{placenta} + V_{amnf}$")), 
                     values = c('1tri_pbtk' = '#7CAE00', 'fetal_pbtk' = '#C77CFF')) +
  geom_vline(xintercept = 13, 
             linetype = 'dashed', color = 'black', linewidth = 1.25) +
  labs(x = 'Gestational Age (weeks)', y = 'Volume (L)', 
       color = "", ) +
  my_theme +
  theme(legend.position = "bottom", 
        legend.direction = "vertical", 
        legend.margin=margin(-15, 0, 0, 0))
  
# check for mass conservation at 13 weeks 
conceptus_density <- 1 # kg/L
full.out2$Wconceptus <- full.out2$fBW + params$placenta_density*full.out2$Vplacenta + params$amnf_density*full.out2$Vamnf
full.out[, Wconceptus := conceptus_density*Vconceptus]

mass_dat <- bind_rows(full.out[, c("time", "Wconceptus")], 
                      full.out2[, c("time", "Wconceptus")])

mass_dat[, time := time/7]
mass_dat.m <- melt.data.table(mass_dat, 
                              id.vars = c("time"), 
                              variable.name = "mass")
mass_dat.m[which(mass_dat.m$time <= 13), model := "1tri_pbtk"]
mass_dat.m[which(mass_dat.m$time > 13), model := "fetal_pbtk"]

masspp <- ggplot(mass_dat.m, 
       aes(x = time, y = value, color = model)) +
  geom_line(linewidth = 2) +
  scale_color_manual(labels = c(TeX("$W_{conceptus}$"), 
                                TeX("$W_{fetus} + W_{placenta} + W_{amnf}$")), 
                     values = c('1tri_pbtk' = '#F8766D', 'fetal_pbtk' = '#00BFC4')) +
  geom_vline(xintercept = 13, 
             linetype = 'dashed', color = 'black', linewidth = 1.25) +
  labs(x = 'Gestational Age (weeks)', y = 'Mass (kg)', 
       color = "") +
  my_theme +
  theme(legend.position = "bottom", 
        legend.direction = "vertical", 
        legend.margin=margin(-15, 0, 0, 0))

# Dynamical Kconceptus2pu parameter --------------------------------------------

# load data from prioritization  
load('./data/invitrodb_v3_5_deiod_filtered_httk_121024.RData', verbose = TRUE)

load_dawson2021() # most data for ToxCast chems
load_sipes2017() # most data for pharma compounds
load_pradeep2020() # best ML method for in silico prediction

parameters <- c("Kconceptus2pu_initial", "Kconceptus2pu_final")

# function to return specific parameter used in 1tri_pbtk
# recall that HTTK functions only work on a single chem at a time 
get_param <- function(dtxsid, parameter) {
  params <- parameterize_1tri_pbtk(dtxsid = dtxsid)
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

# prior to log transforming the data 
df$slopes <- (df$Kconceptus2pu_final-df$Kconceptus2pu_initial)/13

# look at summary stats of slopes
summary(df$slopes)
#> Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> -7715.385   -31.744    -3.553  -268.152    -0.458  6130.769 

# summary rounds up to the nearest int
quantile(df$slopes, c = c(0, 0.25, 0.5, 0.75, 1))
#> 0%           25%           50%           75%          100% 
#> -7715.3846154   -31.7442308    -3.5531923    -0.4576154  6130.7692308 
 
# log10 transform y axis 
df[c("Kconceptus2pu_initial.log10", "Kconceptus2pu_final.log10")] <- lapply(df[c("Kconceptus2pu_initial", "Kconceptus2pu_final")], log10)
df$m <- (df$Kconceptus2pu_final.log10-df$Kconceptus2pu_initial.log10)/13

ylims <- c(floor(min(df$Kconceptus2pu_initial.log10)), ceiling(max(df$Kconceptus2pu_final.log10)))

# show only two example chems with +/- slopes 
Kconc.pp <- ggplot(df[which(df$chnm %in% c("Ergocalciferol", "Terbutylazine")),]) + 
  geom_abline(aes(intercept = Kconceptus2pu_initial.log10, slope = m, color = factor(chnm)), 
              linewidth = 1.5) +
  scale_x_continuous(breaks = seq(0,13), 
                     labels = as.character(seq(0,13)), 
                     limits = c(0,13)) +
  scale_y_continuous(breaks = seq(ylims[1], ylims[2], 1), 
                     limits = ylims) +
  scale_color_manual(values = c('Terbutylazine' = '#FFB900', 'Ergocalciferol' = '#5773CC')) +
  labs(x = "Gestational Age (weeks)", y = TeX("$log_{10} K_{conceptus2pu}$")) +
  my_theme +
  theme(legend.position = "none") +
  annotate("text", x = 11, y = 5.5, label = "Ergocalciferol", 
           size = 5) +
  annotate("text", x = 11, y = 0.6, label = "Terbutylazine", 
           size = 5)

# put it all together
fig <- ggdraw() +
  draw_plot(volpp, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(masspp, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(Kconc.pp, x = 0.10, y = 0, width = 0.75, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), 
                  size = 15, 
                  x = c(0, 0.5, 0), 
                  y = c(1, 1, 0.5))
fig

ggsave(plot = fig, 
       units = "in", 
       dpi = 300, 
       width = 8.1, height = 8.2, 
       device = "tiff", 
       filename = "./figures/300dpi/model_stitching-v2.tiff")

ggsave(plot = fig, 
       units = "in", 
       dpi = 300, 
       width = 8.1, height = 8.2, 
       device = "png", 
       filename = "./figures/model_stitching-v2.png")

