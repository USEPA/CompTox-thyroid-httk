# ==============================================================================
# Script finds critical times at which tissues reach the AC50 related to  
# a specific assay target (DIO enzyme) which I'll call tmax* for short. 
# 
# Tissue concentrations are scaled by the tissue x target-specific AED  
# after assuming dose of 1 mg/kg/day. 
# For BERs < 3, it is the min time (tmax*) at which Cmax*ExpoCast95 > AC50/1000. 
# 
# @author: Kimberly Truong
# created: 2/1/2023
# updated: 11/12/2024
# ==============================================================================

rm(list=ls())
library(httk)
library(dplyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(latex2exp)
library(cowplot)
library(viridis)

# load working function to run the full gestational model
source("./bin/full_pregnancy.R")

load('./data/invitrodb_v3_5_deiod_filtered_httk.RData', verbose = TRUE)

load_dawson2021() # most data for ToxCast chems
load_sipes2017() # most data for pharma compounds 
load_pradeep2020() # best ML method for in silico prediction

# compute half-lives
for (i in 1:nrow(ivive.moe.tb)) {
  ivive.moe.tb$half.life[i] <- calc_half_life(dtxsid = ivive.moe.tb$dtxsid[i]) # in hrs
  ivive.moe.tb$half.life.days <- ivive.moe.tb$half.life/24 # in days
}

ivive.moe.tb$half.life.yrs <- ivive.moe.tb$half.life/24/30/12
ivive.moe.tb$half.life.wks <- ivive.moe.tb$half.life/24/7

# Main Script starts here ------------------------------------------------------

# tissues where DIO enzymes live 
impacted_tissues <- c('Cliver', 'Cthyroid', 'Cfliver', 'Cfbrain', 'Cfthyroid', 
                      'Cplacenta', 'Cconceptus')

# focus on tissues specific to each target 
tissue_list <- list(DIO1 = c('Cliver', 'Cfliver', 'Cthyroid', 'Cfthyroid', 'Cconceptus'), 
                    DIO2 = c('Cthyroid', 'Cfthyroid', 'Cfbrain', 'Cconceptus'), 
                    DIO3 = c('Cthyroid', 'Cfthyroid', 'Cfbrain', 'Cplacenta', 'Cconceptus'), 
                    IYD = c('Cthyroid', 'Cfthyroid', 'Cconceptus'))

# preallocate list of dataframes containing times at which relevant tissues
# reach the AC50 for each DIO target 
crit.times <- list()

for(i in names(tissue_list)) {
  empty.mat <- matrix(NA, nrow = 0, ncol = length(tissue_list[[i]]) + 1)
  empty.df <- data.frame(empty.mat)
  colnames(empty.df) <- c('dtxsid', tissue_list[[i]])

  crit.times[[i]] <- empty.df
}

for(i in 1:nrow(ivive.moe.tb)) {
  
  cat('Calculating for chemical: ', ivive.moe.tb$chnm[i], '\n')
  

  # output solution every hr
  sol.out <- full_pregnancy(dtxsid = ivive.moe.tb$dtxsid[i],
                            daily.dose = 1, 
                            doses.per.day = 1,
                            time.course = seq(0, 40*7, 1/24),
                            track.vars = impacted_tissues)

  
  # function to find time to reach Cmax from scaled output dataframe   
  # acts on each column of the dataframe but also need the entire dataframe 
  find_threshold <- function(v, threshold, conc_df) {
    e = 1e-2 
    ind <- which(abs(v-threshold) < e, arr.ind = TRUE)[1]
    return(conc_df$time[ind])
  }
  
  # for ith chem, find critical times for each target it was found to be 
  # a selective and potent hit for from invitrodb
  for (k in names(tissue_list)) {
    enz.df <- tissue.aeds[[k]]
    chem.rvec <- subset(enz.df, dtxsid == ivive.moe.tb$dtxsid[i]) # this is still a dataframe obj of 1 row 

    if (nrow(chem.rvec) != 0) {
      seem3.dose <- ivive.moe.tb[[i,'exposure95']]
      enz.tissues <- tissue_list[[k]]
      # scaled.res <- sapply(sol.out[, enz.tissues], function(x) x*seem3.dose) # right scaling when looking at priority chems
      
      # scaled by tissue-specific AEDs
      # this covers all chems 
      scaled.res <- mapply(function(x,y) x*y, x = sol.out[, enz.tissues], y = chem.rvec[, enz.tissues]) 
      scaled.res <- as.data.frame(scaled.res) # outputs as a list of columns 
      scaled.res$time <- sol.out$time
      
      # find the tmax* for each relevant tissue (column) associated with a target
      ans <- sapply(scaled.res[, enz.tissues], find_threshold, 
                    threshold = ivive.moe.tb[[i, k]], conc_df = scaled.res)
      
      crit.times[[k]][nrow(crit.times[[k]]) + 1, ] <- c(ivive.moe.tb$dtxsid[i], as.numeric(ans))
  }

  }
  
}


add_chnm_col <- function(df) {
  df$chnm <- ivive.moe.tb$chnm[match(df$dtxsid, ivive.moe.tb$dtxsid)]
  return(df)
}

# add chemical names to each dataframe (DIO target)
crit.times <- lapply(crit.times, add_chnm_col)
chem.cols <- c('dtxsid', 'chnm')

# Wrangle data for eCDF --------------------------------------------------------

melt_times <- function(name, df) {
  df.m <- df %>% 
    reshape2::melt(id.vars = c("dtxsid", "chnm"), 
                   variable.name = "compt", 
                   value.name = "tmax", 
                   na.rm = T)
  df.m$tmax <- as.numeric(df.m$tmax)
  df.m$target <- name 
  return(df.m)
}

new.times.list <- Map(melt_times, names(tissue_list), crit.times)
times.m <- bind_rows(new.times.list)
setDT(times.m)

times.m[, body := "maternal"]
times.m[substr(compt,1,2) == "Cf" | compt %in% c("Cconceptus", "Cplacenta"), body := "fetal"]
View(times.m[, .SD, by = .(dtxsid, body)])

# get the minimum day to reach Cmax for each chem x (maternal, fetal) tissue 
ecdf.data <- times.m[times.m[, .I[tmax == min(tmax)], by = .(dtxsid, body)]$V1]
setnames(ecdf.data, old = "tmax", new = "min_tstar")

ecdf.data$compt <- as.character(ecdf.data$compt)
ecdf.data[, tissues := paste(unique(gsub("C[f]*", "", compt)), collapse = "; "), by = .(dtxsid, body)]
ecdf.data <- unique(ecdf.data[, `:=` (compt = NULL, target = NULL)])

# convert to weeks 
ecdf.data[, min_tstar := (min_tstar/7)]

# reuse theme aesthetics
my_theme <-  theme_bw() + 
  theme(title = element_text(size = 14),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.background = element_blank()) 

my_viridis_colors <- viridis(n = 3)

# plot ECDF for times to reach Cmax in maternal vs fetal tissues
cdfpp <- ggplot(ecdf.data, aes(min_tstar, color = body)) + 
  stat_ecdf(linewidth = 1) +
  scale_x_continuous(limits = c(0, 40), 
                     breaks = c(seq(0,10,2), 13, seq(15,40,5))) +
  scale_color_manual(labels = c(TeX("$min(T^{*}_{f})$ Min week to reach Cmax in the conceptus, placenta, \n or fetal compartment (thyroid, liver, or brain)"), 
                                TeX("$min(T^{*}_{m})$: Min week to reach Cmax in the liver or thyroid of the mother")), 
                     values = c(my_viridis_colors[2], my_viridis_colors[1])) +
  geom_vline(xintercept = 13, 
             linetype = 'dashed', color = my_viridis_colors[3], linewidth = 1) +
  my_theme + 
  theme(legend.position = "bottom", 
        legend.direction = "vertical", 
        legend.title = element_blank()) +
  labs(x = "Minimum Time to reach Cmax by Chemical (weeks of gestation)", 
       y = "Cumulative Frequency") 

# Toxicokinetic Properties for Chems achieving Cmax in 1st vs 2nd Trimester-----

# check the rank order of each CDF curve 
fetal.order <- ecdf.data[body == "fetal"][order(min_tstar)]
mat.order <- ecdf.data[body == "maternal"][order(min_tstar)]

# let's see if chems are distinct toxicokinetically (1st vs 2nd trimester)
set1 <- ecdf.data[body == "maternal" & min_tstar <= 13, dtxsid]
set2 <- ecdf.data[body == "maternal" & min_tstar > 13, dtxsid]

parameters <- c("Clint", "Funbound.plasma", "Fraction_unbound_plasma_fetus", "Pow")

# function to return specific parameter used in fetal_pbtk
# recall that HTTK functions only work on a single chem at a time 
get_param <- function(dtxsid, parameter) {
  params <- parameterize_fetal_pbtk(dtxsid = dtxsid)
  return(params[[parameter]])
}

# this was done once and saved in the RData file 
physchem.tb <- matrix(0, ncol = length(parameters), nrow = nrow(ivive.moe.tb))

for(i in 1:nrow(ivive.moe.tb)) {
  physchem.tb[i, ] = sapply(parameters, get_param, dtxsid = ivive.moe.tb$dtxsid[i])
  cat("Done for", i, "chemicals \n")
}

colnames(physchem.tb) <- parameters
physchem.tb <- data.frame(physchem.tb)
physchem.tb$dtxsid <- ivive.moe.tb$dtxsid
physchem.tb$chnm <- ivive.moe.tb$chnm
physchem.tb <- physchem.tb[, c("dtxsid", "chnm", parameters)]

# make everything on a log10 scale
physchem.tb[parameters] <- lapply(physchem.tb[parameters], log10)

# melt data for violin plots
ecdf.tk <- reshape2::melt(physchem.tb[, c("dtxsid", parameters[parameters != "Fraction_unbound_plasma_fetus"])], 
                          id.vars = c("dtxsid"), 
                          variable.name = "TK_prop", 
                          value.name = "value")
ecdf.tk <- ecdf.tk %>% 
  rowwise() %>%
  mutate(flag = case_when(
    dtxsid %in% set1 ~ "set1", 
    dtxsid %in% set2 ~ "set2"
  ))

setDT(ecdf.tk)

# Clint = 0 => log10(Clint) = -Inf <= -5 
ecdf.tk[TK_prop == "Clint" & value == "-Inf", value := -5]

# make the violin plots
vpp <- ggplot(ecdf.tk, aes(x = TK_prop, y = value, fill = flag)) +
  geom_violin(alpha = 0.5, 
              draw_quantiles = c(0.5)) +
  geom_point(position = position_jitterdodge(dodge.width = 1), 
             alpha = 0.8, 
             show.legend = F) +
  scale_fill_manual(labels = c(TeX("$min(T^{*}_{m}) \\leq $ 13 weeks"), 
                               TeX("$min(T^{*}_{m}) >$ 13 weeks")), 
                      values = c("#E1BEE7", "#4A148C") ) +
  labs(x = "", 
       y = "Value (log10 units)", 
       fill = "Chemical Sets") +
  my_theme +
  guides(fill = guide_legend(position = "inside")) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.82, hjust = 0.80), 
        legend.position.inside = c(0.17, 0.85))
vpp 

# put it all together
fig <- ggdraw() + 
  draw_plot(cdfpp, x = 0, y = 0.2, width = 0.5, height = 0.7) +
  draw_plot(vpp, x = 0.5, y = 0.10, width = 0.5, height = 0.8) +
  draw_plot_label(label = c("A", "B"), size = 15, 
                  x = c(0, 0.5), y = c(0.90, 0.90))

fig 

ggsave(plot = fig, 
       units = "in", 
       dpi = 300, 
       width = 13.9, height = 8.1, 
       device = "tiff", 
       filename = "./figures/300dpi/ecdf-v2.tiff")

ggsave(plot = fig, 
       units = "in", 
       dpi = 300, 
       width = 13.9, height = 8.1, 
       device = "png", 
       filename = "./figures/ecdf-v2.png")

# update RData file with min times to reach Cmax in fetal vs. maternal tissues
# as well as physicochemical property values
e <- new.env(parent = emptyenv())
load('./data/invitrodb_v3_5_deiod_filtered_httk.RData', envir = e)
e$ecdf.data <- ecdf.data
e$physchem.tb <- physchem.tb
do.call("save", c(ls(envir = e), list(envir = e, file ='./data/invitrodb_v3_5_deiod_filtered_httk.RData')))
