# ==============================================================================
# Script checks if conceptus/fetal concentrations (Cfx, Cconceptus, Cfplasma)
# are ever higher than maternal plasma (Cplasma). 
# 
# @author: Kimberly Truong
# created: 5/30/24
# updated: 3/12/25
# ==============================================================================

rm(list=ls())
library(httk)
library(dplyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(latex2exp)
library(cowplot)

load('./data/invitrodb_v3_5_deiod_filtered_httk.RData', verbose = TRUE)

load_dawson2021() # most data for ToxCast chems
load_sipes2017() # most data for pharma compounds 
load_pradeep2020() # best ML method for in silico prediction

# Main Script starts here ------------------------------------------------------

# tissues where DIO enzymes live 
tissues <- c('Cfliver', 'Cfbrain', 'Cfthyroid', 'Cconceptus', 'Cplasma', 'Cfplasma')

# initialize empty dataframe containing the max difference from the Cplasma for 
# every given chemical x tissue 
max.diff <-  data.frame()

for (i in 1:nrow(ivive.moe.tb)) {
  
  # output solution every hr
  sol.out <- solve_full_pregnancy(dtxsid = ivive.moe.tb$dtxsid[i],
                            daily.dose = 1, 
                            doses.per.day = 1,
                            time.course = seq(0, 40*7, 1/24),
                            track.vars = tissues, 
                            physchem.exclude = FALSE)

  # this ignores Inf values after division 
  my.max <- function(x) {
    
    if (all(is.finite(x))) {
      return(max(x, na.rm = T))
    }
    else {
      return(max(x[is.finite(x)], na.rm = T))
    }
  } 
    
  diff.out <- mapply(function(x) x / sol.out[, 'Cplasma', drop = FALSE], 
                     x = sol.out[, tissues[tissues != "Cplasma"] ])
  
  diff.df <- data.frame(diff.out)
  
  max.diff[i, paste0(names(apply(diff.df, 2, my.max)), ".ratio") ] <- apply(diff.df, 2, my.max)
  max.diff[i, 'dtxsid'] <- ivive.moe.tb$dtxsid[i]
  max.diff[i, 'chnm'] <- ivive.moe.tb$chnm[i]

}

setDT(max.diff)

# melt for ggplot 
max.diff.m <- melt.data.table(max.diff, 
                                 id.vars = c("dtxsid", "chnm"), 
                                 variable.name = "ftissue", 
                                 value.name = "max_diff")
max.diff.m$ftissue <- as.character(max.diff.m$ftissue)
max.diff.m[, ftissue := str_split_i(ftissue, "\\.", 1) ]
View(max.diff.m[, .SD, by = .(dtxsid, chnm)])

max.diff.m[, max_diff := log10(max_diff)]
max.diff.m[, max.max_diff := max(max_diff), by = .(dtxsid, chnm)] 
max.diff.m <- max.diff.m[order(max.max_diff)]

# update RData file with max differences
e <- new.env(parent = emptyenv())
load('./data/invitrodb_v3_5_deiod_filtered_httk.RData', envir = e)
e$max.diff.m <- max.diff.m
do.call("save", c(ls(envir = e), list(envir = e, file ='./data/invitrodb_v3_5_deiod_filtered_httk.RData')))

# PLOTTING MAX DIFF FROM CPLASMA -----------------------------------------------
boundary.chem <- tail(max.diff.m[max.max_diff < 0.5])[6, chnm] #pirimicarb 

min.maxdiff <- max.diff.m[, min(max_diff)]
max.maxdiff <- max.diff.m[, max(max_diff)]
ylims <- c(floor(min.maxdiff), ceiling(max.maxdiff)-0.5)

my_theme <- theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12), 
        legend.background = element_blank()) 

Cdiff <- ggplot(max.diff.m, 
       mapping = aes(x = reorder(chnm, -max.max_diff), 
                     y = max_diff)) +
  geom_point(aes(shape = ftissue, color = ftissue), size = 3, 
             alpha = 0.7, stroke = 1.5) +
  geom_hline(yintercept = 0.5, color = 'red', linetype = "longdash", linewidth = 1.2) +
  scale_y_continuous(breaks = seq(ylims[1], ylims[2], 1), 
                     limits = ylims) +
  scale_shape_manual(name = "x: Fetal Tissue Concentration", 
                     labels = c("Cfplasma", "Cfbrain", "Cconceptus", "Cfthyroid", "Cfliver"), 
                     values = c('Cfbrain' = 16, 
                                'Cfliver' = 17,
                                'Cfthyroid' = 8, 
                                'Cconceptus' = 7, 
                                'Cfplasma' = 18)) +
  scale_color_manual(name = "x: Fetal Tissue Concentration", 
                     labels = c("Cfplasma", "Cfbrain", "Cconceptus", "Cfthyroid", "Cfliver"), 
                     values = c('Cfbrain' = "#E69F00", 
                                'Cfliver' = "#CC79A7",
                                'Cfthyroid' = "#0072B2", 
                                'Cconceptus' = "#009E73", 
                                'Cfplasma' = "#D55E00")) +
  annotate("rect", xmin = Inf, xmax =  boundary.chem, 
           ymin = -Inf, ymax = Inf, alpha = 0.2) +
  annotate("text", x = 89, y = 2.3, label = "Protected by Cplasma", 
           fontface = "bold", size = 14/.pt) +
  annotate("text", x = 41, y = 2.3, label = "Cftissue > Cplasma", 
           fontface = "bold", size = 14/.pt) +
  my_theme + 
  theme(axis.title = element_text(face = "bold"), 
        legend.position = 'bottom', 
        legend.direction = 'horizontal') + 
  labs(x = 'Chemical', y = "log10(x) - log10 Cplasma") +
  coord_flip()

# Exploring TK properties ------------------------------------------------------

# let's see if higher concs are driven by TK properties
# pare down data to the max difference by chemical 
TK.data <- max.diff.m[max_diff == max.max_diff]
TK.data[, max.max_diff := NULL]

# relative max is achieved in both Cfplasma and Cfliver for Sodium 2-mercaptobenzothiolate (DTXSID1026035)
TK.data[, if(.N > 1) .SD, by = dtxsid]

# Pirimicarb has diff of 0.5
TK.data[, max_diff := signif(max_diff,digits = 2)]
protective.chems <- TK.data[max_diff <= 0.5, unique(dtxsid)]
nonprotective.chems <- TK.data[max_diff > 0.5, unique(dtxsid)]
parameters <- c("Clint", "Funbound.plasma", "Fraction_unbound_plasma_fetus", "Pow")

# melt data for violin plots
tk.datm <- reshape2::melt(physchem.tb[, c("dtxsid", parameters[parameters != "Fraction_unbound_plasma_fetus"])], 
                id.vars = c("dtxsid"), 
                variable.name = "TK_prop", 
                value.name = "value")

tk.datm <- tk.datm %>% 
  rowwise() %>%
  mutate(flag = case_when(
    dtxsid %in% protective.chems ~ "protective", 
    dtxsid %in% nonprotective.chems ~ "nonprotective"
  ))
tk.datm$flag <- factor(tk.datm$flag, levels = c("protective", "nonprotective"))
setDT(tk.datm)

# extract medians
tk.medians <- tk.datm[, .(median = median(value)), by = .(TK_prop, flag)]
tk.medians$med.log10 <- log10(tk.medians$median)
tk.medians$prop_int <- as.numeric(tk.medians$TK_prop)
tk.medians[flag == "protective", prop_int := prop_int - 0.5]
tk.medians[flag == "nonprotective", prop_int := prop_int + 0.03]
setnames(tk.medians, "prop_int", "x1")

tk.datm[, value := log10(value)]
# Clint = 0 => log10(Clint) = -Inf <= -5 
tk.datm[TK_prop == "Clint" & value == "-Inf", value := -5]

# make the violin plots
vpp <- ggplot(tk.datm, aes(x = TK_prop, y = value, fill = flag)) +
  geom_violin(alpha = 0.5) + 
              # draw_quantiles = c(0.5)) +
  geom_point(aes(x = TK_prop, y = value),
             position = position_jitterdodge(dodge.width = 1), 
             alpha = 0.8, 
             show.legend = F) +
  scale_fill_manual(labels = c("Always Protected by Maternal Plasma Concentration", 
                                 "Fetal Tissue Concentration > Maternal Plasma Concentration"), 
                    values = c("#828585", "#DC322F")) +
  guides(fill = guide_legend(position = "inside")) +
  labs(x = "", y = "Value (log10 units)", fill = "") +
  my_theme + 
  theme(legend.position.inside = c(0.35, 0.88), 
        legend.direction = "vertical", 
        axis.text.x = element_text(angle = 60, vjust = 0.82, hjust = 0.80)) +
  geom_segment(
    aes(x = x1, xend = x1 + 0.55, y = med.log10, yend = med.log10),
    color = "#00017D", linetype = "dashed", linewidth = 1, 
    data = tk.medians)

# put it all together
# try plot_grid
# plot_grid(Cdiff, vpp, 
#           rel_widths = c(1,1), 
#           rel_heights = c(1, 0.33), 
#           labels = c("A", "B"), 
#           ncol = 2, nrow = 1)

fig <- ggdraw() + 
  draw_plot(Cdiff, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(vpp, x = 0.55, y = 0.46, width = 0.45, height = 0.5) +
  draw_plot_label(label = c("A", "B"), size = 30, 
                  x = c(0, 0.55), y = c(1, 1))

fig

ggsave(plot = fig, 
       units = "in", 
       dpi = 300, 
       width = 18, height = 16.5, 
       device = "tiff", 
       filename = "./figures/300dpi/is_Cplasma_protective-v3.tiff")

ggsave(plot = fig, 
       units = "in", 
       dpi = 300, 
       width = 18, height = 16.5, 
       device = "png", 
       filename = "./figures/is_Cplasma_protective-v3.png")

