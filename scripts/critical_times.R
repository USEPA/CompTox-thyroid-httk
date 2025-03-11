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
# updated: 3/7/2025
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
library(openxlsx)

load('./data/invitrodb_v3_5_deiod_filtered_httk.RData', verbose = TRUE)

load_dawson2021() # most data for ToxCast chems
load_sipes2017() # most data for pharma compounds 
load_pradeep2020() # best ML method for in silico prediction

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

  # output solution every hr
  sol.out <- solve_full_pregnancy(dtxsid = ivive.moe.tb$dtxsid[i],
                            daily.dose = 1, 
                            doses.per.day = 1,
                            time.course = seq(0, 40*7, 1/24),
                            track.vars = impacted_tissues, 
                            physchem.exclude = FALSE)

  
  # function to find time to reach Cmax from scaled output dataframe   
  # acts on each column of the dataframe but also need the entire dataframe 
  find_threshold <- function(v, threshold, conc_df) {
    e = 1e-1
    ind <- which(abs(v-threshold) < e, arr.ind = TRUE)[1]
    return(conc_df$time[ind])
  }
  
  # for ith chem, find critical times for each target it was found to be 
  # a selective and potent hit for from invitrodb
  for (k in names(tissue_list)) {
    enz.df <- tissue.aeds[[k]]
    chem.rvec <- subset(enz.df, dtxsid == ivive.moe.tb$dtxsid[i]) # this is still a dataframe obj of 1 row 

    if (nrow(chem.rvec) != 0) {
      seem3.dose <- ivive.moe.tb[[i,'seem3.u95']]
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
Cmax.times.m <- bind_rows(new.times.list)
setDT(Cmax.times.m)

Cmax.times.m[, lifestage := "maternal"]
Cmax.times.m[substr(compt,1,2) == "Cf" | compt %in% c("Cconceptus", "Cplacenta"), lifestage := "fetal"]
View(Cmax.times.m[, .SD, by = .(dtxsid, lifestage)])

# add this to Table 7
targets <- c('DIO1', 'DIO2', 'DIO3', 'IYD')
plasma.tissue.bers[, target := mapply(function(x) targets[x], 
                                      apply(plasma.tissue.bers[, paste0(targets, "_BER"), with = FALSE], 1, which.min))]
setnames(plasma.tissue.bers, "most_sensitive_tissue", "compt")

table7 <- merge.data.table(plasma.tissue.bers, Cmax.times.m[, -c("lifestage")], 
                           by = c("dtxsid", "chnm", "target", "compt"), 
                           all.x = TRUE) 

# collapse multiple sensitive tissues and their tmaxes (time to reach Cmax in each tissue)
table7[, most_sensitive_tissues := paste(compt, collapse = ", ") , by = .(dtxsid, target)]
table7[, tmax := round(tmax, digits = 1)]
table7[, day_at_Cmax := paste(tmax, collapse = ", "), by = .(dtxsid, target)]

keep_cols <- c("dtxsid", "chnm", "plasma_BER50", 
               paste0(targets, "_BER"), "min_BER", 
               "most_sensitive_tissues", "day_at_Cmax")
table7 <- unique(table7[, ..keep_cols])
table7 <- table7[order(min_BER)]
write.xlsx(table7, "./tables/preg_vs_nonpregBERs_v3.xlsx", colnames = T)

# get the minimum day to reach Cmax for each chem x (maternal, fetal) tissue 
ecdf.data <- Cmax.times.m[Cmax.times.m[, .I[tmax == min(tmax)], by = .(dtxsid, lifestage)]$V1]
setnames(ecdf.data, old = "tmax", new = "min_tstar")

ecdf.data$compt <- as.character(ecdf.data$compt)
ecdf.data[, tissues := paste(unique(gsub("C[f]*", "", compt)), collapse = "; "), by = .(dtxsid, lifestage)]
ecdf.data <- unique(ecdf.data[, `:=` (compt = NULL, target = NULL)])

# convert to weeks 
ecdf.data[, min_tstar := (min_tstar/7)]

# PLOTTING ECDF CURVES ---------------------------------------------------------
# reuse theme aesthetics
my_theme <-  theme_bw() + 
  theme(title = element_text(size = 14),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.background = element_blank()) 

my_viridis_colors <- viridis(n = 3)

# plot ECDF for times to reach Cmax in maternal vs fetal tissues
cdfpp <- ggplot(ecdf.data, aes(min_tstar, color = lifestage)) + 
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

# How many chems reached Cmax in the mother by 1st tri?  (88/103 = ~85%)
ecdf.data[lifestage == "maternal" & min_tstar < 13, length(unique(dtxsid))]
#> [1] 88

# 15 chems reached Cmax in the mother in the 2nd-3rd tri 
# but 13 of these chems reached Cmax in 37th-39th weeks
View(ecdf.data[min_tstar > 13 & lifestage == "maternal"])

# Toxicokinetic Properties for Chems achieving Cmax in 1st vs 2nd Trimester-----

# let's see if chems are distinct toxicokinetically (1st vs 2nd trimester)
set1 <- ecdf.data[lifestage == "maternal" & min_tstar <= 13, dtxsid]
set2 <- ecdf.data[lifestage == "maternal" & min_tstar > 13, dtxsid]

parameters <- c("Clint", "Funbound.plasma", "Fraction_unbound_plasma_fetus", "Pow")

# function to return specific parameter used in fetal_pbtk
# recall that HTTK functions only work on a single chem at a time
get_param <- function(dtxsid, parameter) {
  params <- parameterize_fetal_pbtk(dtxsid = dtxsid, physchem.exclude = FALSE)
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

# capture medians before log10 transformation
tk.medians <- ecdf.tk[, .(median = median(value)), by = .(TK_prop, flag)]
tk.medians$median.log <- log10(tk.medians$median)
tk.medians[median.log == -Inf, median.log := -5]
tk.medians$prop_int <- as.numeric(tk.medians$TK_prop)
tk.medians[flag == "set1", prop_int := prop_int - 0.5]
tk.medians[flag == "set2", prop_int := prop_int + 0.03]
setnames(tk.medians, "prop_int", "x1")

# convert to log10 scale
ecdf.tk[, value := log10(value)]

# Clint = 0 => log10(Clint) = -Inf <= -5 
ecdf.tk[TK_prop == "Clint" & value == "-Inf", value := -5]

# make the violin plots
vpp <- ggplot(ecdf.tk, aes(x = TK_prop, y = value, fill = flag)) +
  geom_violin(alpha = 0.5) + 
              # draw_quantiles = c(0.5)) +
  geom_point(position = position_jitterdodge(dodge.width = 1),
             alpha = 0.8,
             show.legend = F) +
  # geom_boxplot(position = position_dodge(width = 0.9), 
  #   width=0.1, color="black", alpha=0.5) +
  scale_fill_manual(labels = c(TeX("$min(T^{*}_{m}) \\leq $ 13 weeks"), 
                               TeX("$min(T^{*}_{m}) >$ 13 weeks")), 
                      values = c("#E1BEE7", "#4A148C") ) +
  labs(x = "", 
       y = "Value (log10 units)", 
       fill = "Chemical Sets") +
  my_theme +
  guides(fill = guide_legend(position = "inside")) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.82, hjust = 0.80), 
        legend.position.inside = c(0.17, 0.85)) +
  geom_segment(
               aes(x = x1, xend = x1 + 0.55, y = median.log, yend = median.log),
               color = "red", linetype = "dashed", linewidth = 1, 
               data = tk.medians)
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
       filename = "./figures/300dpi/ecdf-v3.tiff")

ggsave(plot = fig, 
       units = "in", 
       dpi = 300, 
       width = 13.9, height = 8.1, 
       device = "png", 
       filename = "./figures/ecdf-v3.png")

# update RData file with min times to reach Cmax in fetal vs. maternal tissues
# as well as physicochemical property values
e <- new.env(parent = emptyenv())
load('./data/invitrodb_v3_5_deiod_filtered_httk.RData', envir = e)
e$table7 <- table7
e$Cmax.times.m <- Cmax.times.m
e$ecdf.data <- ecdf.data
e$physchem.tb <- physchem.tb
do.call("save", c(ls(envir = e), list(envir = e, file ='./data/invitrodb_v3_5_deiod_filtered_httk.RData')))

# Investigating Chems that reach Cmax in the Mother End of Term ----------------

Cmax.term.chems <- ecdf.data[min_tstar > 13 & lifestage == "maternal", dtxsid]

# 13 out of 15 chems achieved Cmax in the 37th-39th weeks 
# exceptions: 2-Chloro-N-phenylacetamide and Quinoxyfen
View(ecdf.data[lifestage == "maternal" & dtxsid %in% Cmax.term.chems])

# 10 of these chems have Clint == 0 and all chems except 
# 2-Chloro-N-phenylacetamide have Fup ~ 0
View(physchem.tb[physchem.tb$dtxsid %in% Cmax.term.chems,])

# examine some simulations
# Kepone (DTXSID1020770): Clint = 0; Fup = 0
# 2,2',6,6'-Tetrachlorobisphenol A (DTXSID3021770):  Clint = 511.1; fup = 0
test.chem <- c("DTXSID3021770")

# output solution every hr
sol.out <- solve_full_pregnancy(dtxsid = test.chem,
                                daily.dose = 1, 
                                doses.per.day = 1,
                                time.course = seq(0, 40*7, 1/24),
                                track.vars = impacted_tissues, 
                                physchem.exclude = FALSE)


# function to find time to reach Cmax from scaled output dataframe   
# acts on each column of the dataframe but also need the entire dataframe 
find_threshold <- function(v, threshold, conc_df) {
  e = 1e-2 
  ind <- which(abs(v-threshold) < e, arr.ind = TRUE)[1]
  return(conc_df$time[ind])
}

# for ith chem, find critical times for each target it was found to be 
# a selective and potent hit for from invitrodb
for (k in c("DIO3")) {
  browser()
  enz.df <- tissue.aeds[[k]]
  chem.rvec <- subset(enz.df, dtxsid == test.chem) # this is still a dataframe obj of 1 row 
  
  if (nrow(chem.rvec) != 0) {
    seem3.dose <- ivive.moe.tb[[which(ivive.moe.tb$dtxsid == test.chem),'seem3.u95']]
    enz.tissues <- tissue_list[[k]]
    # scaled.res <- sapply(sol.out[, enz.tissues], function(x) x*seem3.dose) # right scaling when looking at priority chems
    
    # scaled by tissue-specific AEDs
    # this covers all chems 
    scaled.res <- mapply(function(x,y) x*y, x = sol.out[, enz.tissues], y = chem.rvec[, enz.tissues]) 
    scaled.res <- as.data.frame(scaled.res) # outputs as a list of columns 
    scaled.res$time <- sol.out$time
    
    # find the tmax* for each relevant tissue (column) associated with a target
    ans <- sapply(scaled.res[, enz.tissues], find_threshold, 
                  threshold = ivive.moe.tb[[which(ivive.moe.tb$dtxsid == test.chem), k]], conc_df = scaled.res)
  }
  
}
