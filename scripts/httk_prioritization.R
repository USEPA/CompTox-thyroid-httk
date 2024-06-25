# ==============================================================================
# This script plots the results of
# applying the full gestational model in httk for IVIVE. 
# 
# This script generates Figures 8-10 of Truong et al. 2024.  
# @author: Kimberly Truong 
# created: 12/7/23
# updated: 6/25/24
# ==============================================================================

rm(list=ls())
library(dplyr)
library(tidyverse)
library(data.table)
library(openxlsx)
library(ggpubr)
library(reshape2) # pkg to melt data frame
library(ggrepel)
library(latex2exp)
library(cowplot)

load('./data/Deiodinases/deiod_invitrodbv3_5_filtered_ivive_fullterm_preg_branch.RData', 
     verbose = TRUE)

# tissues where DIO enzymes live 
impacted_tissues <- c('Cliver', 'Cthyroid', 'Cfliver', 'Cfbrain', 'Cfthyroid', 
                      'Cplacenta', 'Cconceptus')

targets <- c('DIO1', 'DIO2', 'DIO3', 'IYD')

# Wrangle Data for ggplot ------------------------------------------------------

get_ber_data <- function(name, df) {
  
  df$exposure <- ivive.moe.tb$exposure95[match(df$dtxsid, ivive.moe.tb$dtxsid)]
  df.m <- df %>% select(-exposure) %>%
    reshape2::melt(id.vars = c('dtxsid'), 
                   variable.name = 'compt', value.name = 'aed')
  df.m$exposure <- ivive.moe.tb$exposure95[match(df.m$dtxsid, ivive.moe.tb$dtxsid)]
  df.m$ber = df.m$aed/df.m$exposure
  df.m2 <- df.m %>% 
    select(-exposure) %>%
    reshape2::melt(id.vars = c('dtxsid', 'compt', 'ber'))
  
  # add melted ExpoCast data as rows for each unique chem 
  expocast.dat <- df %>% select(dtxsid, exposure) %>% 
    reshape2::melt(id.vars = c('dtxsid'))
  
  df.m2 <- bind_rows(df.m2, expocast.dat)
  df.m2$value <- log10(df.m2$value)
  df.m2$ber <- log10(df.m2$ber)
  df.m2$chnm <- ivive.moe.tb$chnm[match(df.m2$dtxsid, ivive.moe.tb$dtxsid)]
  df.m2$enzyme <- name
  return(df.m2)
}

new.list <- Map(get_ber_data, names(tissue.aeds), tissue.aeds)
new_dat= bind_rows(new.list) 
setDT(new_dat)

# order chemicals by the most sensitive tissue endpoint 
new_dat[, min_ber := min(ber, na.rm = TRUE), keyby = .(enzyme, dtxsid)]

# exposure values are dependent on dtxsid not httk compt 
new_dat[variable == 'exposure', compt := 'exposure']
new_dat[variable == 'aed' & min_ber == ber, most_sensitive_tissue := compt, by = .(enzyme, dtxsid)]

new_dat[compt %in% impacted_tissues[grep('f', impacted_tissues)], body := 'fetal']
new_dat[compt %in% c('Cliver', 'Cthyroid'), body := 'maternal']
new_dat[compt %in% c('Cplacenta', 'Cconceptus'), body := 'conceptus']
new_dat[compt == "exposure", body := "exposure"]

new_dat[compt %in% impacted_tissues[grep('liver', impacted_tissues)], 
        tissue := 'liver']
new_dat[compt %in% impacted_tissues[grep('thyroid', impacted_tissues)], 
        tissue := 'thyroid']
new_dat[compt == 'Cfbrain', tissue := 'brain']
new_dat[compt == 'Cplacenta', tissue := 'placenta']
new_dat[compt == "Cconceptus", tissue := 'conceptus']
new_dat[compt == 'exposure', tissue := 'exposure'] 

# convert tissue and body variables to factors
new_cols <- c('body', 'tissue')
new_dat[, (new_cols) := lapply(.SD, factor), .SDcols = new_cols]
new_dat <- new_dat[order(min_ber), .SD, keyby = .(enzyme)]

min_val <- new_dat[min_ber < 3, min(value), by = enzyme][, min(V1)]
max_val <- new_dat[value != "Inf" & min_ber < 3, max(value), by = enzyme][, max(V1)]
ylims <- c(floor(min_val), ceiling(max_val))
  
listgg <- list()
for(k in names(tissue.aeds)) {
  
  aeds <- new_dat[enzyme == k & ber < 3]
  chems <- new_dat[enzyme == k & ber < 3, dtxsid]
  exposures <- new_dat[enzyme == k & variable == "exposure" & dtxsid %in% chems]
  
  listgg[[k]] <- ggplot(data = rbind(aeds, exposures), 
                        mapping = aes(x = reorder(chnm, -min_ber), 
                                      y = value)) +
    geom_point(aes(shape = tissue, color = body), 
               alpha = 0.65, size = 3, stroke = 1.5, 
               show.legend = T) +
    scale_shape_manual(values = c('brain' = 19,
                                  'exposure' = 18, 
                                  'liver' = 17,
                                  'placenta' = 15, 
                                  'thyroid' = 8, 
                                  'conceptus' = 7), 
                       drop = F) + # this makes sure all levels of a factor are displayed even if unused 
    scale_color_manual(values = c('maternal' = 'darkviolet', 'fetal' = 'green', 
                                  'conceptus' = 'pink', 'exposure' = "dimgrey")) +
    labs(x = 'Chemical', y = 'Oral AED or Exposure (log10 mg/kg-bw/day)', title = k) +
    theme_bw() +
    theme(legend.position = 'bottom',
          legend.direction = 'horizontal', 
          plot.title = element_text(size = 18, face = "bold"), 
          axis.title = element_text(size = 16), 
          legend.title = element_blank(), 
          axis.text = element_text(size = 14), 
          legend.text = element_text(size = 14)) +
    scale_y_continuous(breaks = seq(ylims[1], ylims[2], 1), 
                       limits = ylims) +
    coord_flip()
}

# remove extra axis labels 
listgg$IYD <- listgg$IYD + labs(x = '', y = '')
listgg$DIO3 <- listgg$DIO3 + 
  labs(x = '') + 
  theme(axis.text.y = element_text(size = 13))

listgg$DIO1 <- listgg$DIO1 + labs(y = '')
bers.lessthan3 <- ggarrange(listgg$DIO1, listgg$IYD, listgg$DIO2, listgg$DIO3, 
          nrow = 2, ncol = 2, 
          heights = c(1, 1, 0.1, 1), 
          common.legend = TRUE, legend = 'bottom', 
          align = "hv")

ggsave(plot = bers.lessthan3, 
       units = "in", 
       dpi = 300, 
       width = 21.33, height = 11, 
       device = "tiff", 
       filename = "./doc/comptox-thyroid-httk/figures/preg_prioritization_by_BER.tiff")

ggsave(plot = bers.lessthan3, 
       units = "in", 
       dpi = 300, 
       width = 21.33, height = 11, 
       device = "png", 
       filename = "./doc/comptox-thyroid-httk/figures/preg_prioritization_by_BER.png")

# update RData file with melted data for ggplot
# i.e. chemicals ranked by tissue-based BERs with exposure estimates from ExpoCast
# e <- new.env(parent = emptyenv())
# load('./data/Deiodinases/deiod_invitrodbv3_5_filtered_ivive_fullterm_preg_branch.RData', envir = e)
# e$tissue.bers <- new_dat
# do.call("save", c(ls(envir = e), list(envir = e, file = './data/Deiodinases/deiod_invitrodbv3_5_filtered_ivive_fullterm_preg_branch.RData')))

# Interpreting BERs ------------------------------------------------------------

# table matching scatterplot subplots with all calculated BERs
# every enzyme x chem x min BER constitutes a row 
enzxchem.tb <- new_dat[!is.na(most_sensitive_tissue), 
                       .(chnm, min_ber, most_sensitive_tissue, body, tissue), by = .(enzyme, dtxsid)]
# rows ordered by enzyme and within each enzyme group, ordered by min_ber 
enzxchem.tb.sorted <- enzxchem.tb[order(min_ber), .SD, keyby = .(enzyme)] 

# collapse multiple most-sensitive tissues 
# this really only applies to Ergocalciferol 
melt.foreach.dio <- enzxchem.tb.sorted[, most_sensitive_compts := paste(most_sensitive_tissue, collapse = ", "), 
                                       by = .(enzyme, dtxsid)][, .(chnm, min_ber, most_sensitive_compts), by = .(enzyme, dtxsid)]
melt.foreach.dio2 <- unique(melt.foreach.dio)

# min BER matrix
min.ber.mat <- dcast(melt.foreach.dio2[, 1:4], 
                     dtxsid + chnm ~ enzyme, value.var = c('min_ber'))
setDT(min.ber.mat)
min.ber.mat[, min_across_enzyme := pmin(DIO1, DIO2, DIO3, IYD, na.rm = T)]
min.ber.mat <- min.ber.mat[order(min_across_enzyme)]

# find the target corresponding to the min_across_enzyme
min.ber.mat[, enzyme := mapply(function(x) targets[x], apply(min.ber.mat[, targets, with = FALSE], 1, which.min))]

# merge data on most sensitive tissues reflecting min BER 
# this ensures an extra row for Ergocalciferol to represent the multiple sensitive tissue endpoints 
new.ranking.tissues <- merge.data.table(min.ber.mat, 
                                 enzxchem.tb.sorted, 
                                 by = c('dtxsid', 'chnm', 'enzyme'), 
                                 all.x = TRUE)
setnames(new.ranking.tissues, old = 'enzyme', new = 'target_min')

# round all BERs to 2 decimal places for easy table viewing 
ranked.by.plasma50 <- ivive.moe.tb[, lmoe50 := round(log_pBER50, digits = 2)][order(lmoe50)]

final.new.ranking <- new.ranking.tissues[, -c("target_min", "min_ber", "most_sensitive_compts")]
final.new.ranking[, c(targets, "min_across_enzyme") := lapply(.SD, round, digits = 2), 
                  by = seq_len(nrow(final.new.ranking)), .SDcols = c(targets, "min_across_enzyme")]
final.new.ranking <- final.new.ranking[order(min_across_enzyme)]

# comparing tissue BER with plasma BER from 3compss model 
plasma.tissue.bers <- merge.data.table(ranked.by.plasma50[, c('dtxsid', 'chnm', 'lmoe50')], 
                           final.new.ranking[, -c('chnm')], 
                           by = c('dtxsid'))

setnames(plasma.tissue.bers, old = colnames(plasma.tissue.bers)[3:8], 
         new = c('plasma_BER50', 'DIO1_BER', 'DIO2_BER', 'DIO3_BER', 'IYD_BER', 'min_BER'))

plasma.tissue.bers <- plasma.tissue.bers[order(min_BER)]
write.xlsx(plasma.tissue.bers, "./doc/comptox-thyroid-httk/tables/plasma_and_tissueBERs.xlsx", colnames = T)

# hide double label of Ergocalciferol 
# this occurred because min AED for this chemical was achieved in both maternal/fetal thyr
plasma.tissue.bers[which(plasma.tissue.bers$chnm == "Ergocalciferol")[2], chnm := NA_character_]

max_ber <- max(plasma.tissue.bers[, c('plasma_BER50', 'DIO1_BER', 'DIO2_BER', 'DIO3_BER', 'IYD_BER')], na.rm = T)
min_ber <- min(plasma.tissue.bers[, c('plasma_BER50', 'DIO1_BER', 'DIO2_BER', 'DIO3_BER', 'IYD_BER')], na.rm = T)
ber.lims <- c(floor(min_ber), ceiling(max_ber))

# theme for the next 2 figures 
my_theme <- theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12))
  
preg.vs.nonpreg.pp <- ggplot(data = plasma.tissue.bers, aes(x = plasma_BER50, y = min_BER)) +
  geom_point(alpha = 0.75, size = 3, 
             show.legend = T) +
  scale_shape_manual(values = c('brain' = 19,
                                'liver' = 17,
                                'placenta' = 15, 
                                'thyroid' = 8, 
                                'conceptus' = 7), 
                     drop = F) + # this makes sure all levels of a factor are displayed even if unused 
  scale_color_manual(values = c('maternal' = 'darkviolet', 
                                'fetal' = 'forestgreen', 
                                conceptus = 'pink')) +
  geom_abline(slope = 1, linewidth = 1, linetype = 'solid', color = "#666666") +
  geom_abline(aes(slope = 1, intercept = 1), linetype = 'longdash', color = 'grey', linewidth = 1) + # default linesize = 1
  geom_abline(aes(slope = 1, intercept = -1), linetype = 'longdash', color = 'grey', linewidth = 1) +
  geom_label_repel(aes(label = chnm), 
                   size = 5, 
                   force = 28) +
  my_theme +
  theme(legend.position = 'bottom', 
        legend.direction = 'horizontal',
        legend.title = element_blank()) + 
  scale_x_continuous(breaks = seq(ber.lims[1], ber.lims[2], 1), 
                     limits = ber.lims) +
  scale_y_continuous(breaks = seq(ber.lims[1], ber.lims[2], 1), 
                     limits = ber.lims) +
  labs(x = 'plasma BER using 3compss model', y = 'min BER across DIO/IYD targets and tissues')

# plot SEEM3 uncertainty
ivive.moe.tb$uncertainty = log10(ivive.moe.tb$exposure95/ivive.moe.tb$exposure_med)
plasma.tissue.bers$uncertainty <- ivive.moe.tb$uncertainty[match(plasma.tissue.bers$dtxsid, ivive.moe.tb$dtxsid)]

min.min.ber <- min(plasma.tissue.bers$min_BER)
max.min.ber <- max(plasma.tissue.bers$min_BER)
min.ber.lims <- c(floor(min.min.ber), ceiling(max.min.ber))
max.uncertainty <- max(plasma.tissue.bers$uncertainty)

seem3pp <- ggplot(data = plasma.tissue.bers[-3, ], 
       mapping = aes(x = min_BER, 
                     y = uncertainty)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_label_repel(aes(label = chnm),
                   data = subset(plasma.tissue.bers, min_BER < 0 | uncertainty > 4),
                   force = 25,
                   size = 5) +
  my_theme +
  theme(legend.position = 'bottom', 
        legend.direction = 'horizontal',
        legend.title = element_blank()) + 
  scale_x_continuous(breaks = seq(min.ber.lims[1], min.ber.lims[2], 1), 
                     limits = min.ber.lims) +
  scale_y_continuous(breaks = seq(0, ceiling(max.uncertainty), 1), 
                     limits = c(0, ceiling(max.uncertainty))) +
  labs(x = 'min BER by Chemical', y = TeX(r'($log_{10}ExpoCast_{u95} - log_{10}ExpoCast_{med}$)'))

median(plasma.tissue.bers$min_BER)
#> [1] 3.375

bers2 <- plot_grid(preg.vs.nonpreg.pp, seem3pp, 
          labels = "AUTO", label_size = 20)

ggsave(plot = bers2, 
       units = "in", 
       dpi = 300, 
       width = 18.53, height = 8.32, 
       device = "tiff", 
       filename = "./doc/comptox-thyroid-httk/figures/BER_uncertainty.tiff")

ggsave(plot = bers2, 
       units = "in", 
       dpi = 300, 
       width = 18.53, height = 8.32, 
       device = "png", 
       filename = "./doc/comptox-thyroid-httk/figures/BER_uncertainty.png")

