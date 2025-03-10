# ==============================================================================
# This script compares the DIO-bioactivity-based PODs with ToxValDB and 
# the 5th percentile PODs (AEDs) obtained from all of ToxCast. 
# 
# @author: Kimberly Truong 
# created: 12/7/23
# updated: 3/10/2025
# ==============================================================================

rm(list=ls())
library(dplyr)
library(tidyverse)
library(data.table)
library(openxlsx)
library(ggpubr)
library(reshape2)
library(ggrepel)
library(latex2exp)
library(cowplot)

load('./data/invitrodb_v3_5_deiod_filtered_httk.RData', 
     verbose = TRUE)

targets <- c("DIO1", "DIO2", "DIO3", "IYD")

# Plot ToxValDB PODs vs. min AEDs by Chemical ----------------------------------
pdata <- tissue.bers[variable == 'aed', .(dtxsid, chnm, variable, value)]
pdata[, min_aed := min(value), by = .(dtxsid)]
pdata <- pdata[order(dtxsid)]

# take the min AED for each of the 103 chems 
pdata2 <- unique(pdata[, .(dtxsid, chnm, min_aed)])

# add in PODs from ToxValDB, if available 
toxval <- read.xlsx('./tables/toxval pods chemical level oral mgkgday.xlsx')
toxval <- toxval %>% filter(dtxsid %in% pdata2$dtxsid)
pdata2$pod <- log10(toxval$pod[match(pdata2$dtxsid, toxval$dtxsid)])
pdata3 <- pdata2[which(!is.na(pdata2$pod)), ] 

# 66 chems with complementary ToxValDB POD
pdata3[, length(unique(dtxsid))]
#> [1] 66

pod.min <- min(pdata3$min_aed, pdata3$pod)
pod.max <- max(pdata3$min_aed, pdata3$pod)
pod.lims <- c(floor(pod.min), ceiling(pod.max))

# to standardize themes for both plots 
my_theme <- theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12))
  
toxval.pp <- ggplot(data = pdata3, aes(x = min_aed, y = pod)) +
  geom_point(size = 3, alpha = 0.75) + # default size = 1.5 
  geom_abline(slope = 1, linewidth = 1, linetype = 'solid') +
  geom_abline(aes(slope = 1, intercept = 1), 
              linetype = 'longdash', color = 'grey', linewidth = 1) + # default linesize = 1
  geom_abline(aes(slope = 1, intercept = -1), 
              linetype = 'longdash', color = 'grey', linewidth = 1) +
  geom_label_repel(aes(label = chnm), 
                   data = pdata3[abs(pod - min_aed) > 1],
                   size = 4, 
                   force = 65, max.overlaps = 9) +
  my_theme +
  scale_x_continuous(breaks = seq(pod.lims[1], pod.lims[2], 1), 
                     limits = pod.lims) +
  scale_y_continuous(breaks = seq(pod.lims[1], pod.lims[2], 1), 
                     limits = pod.lims) +
  labs(x = 'min DIO/IYD-bioactivity-based AED (log10 mg/kg/day)', 
       y = 'ToxVal POD (log10 mg/kg/day)')

# interesting chems above grey line 
pdata3 %>% 
  filter(min_aed + 1 < pod) %>% 
  select(dtxsid, chnm) %>% print()

chems.to.investigate <- pdata3[min_aed + 1 < pod, dtxsid]
# print source(s) underlying ToxValDB pod 
toxval.studies <- toxval %>% filter(dtxsid %in% chems.to.investigate) %>% 
  select(name, studies, pod_10, pod_hra, pod_hra_source, pod) %>% print()

# DIO-based AEDs vs. 5th %-tile ToxCast AC50 => AEDs using 3compss -------------

# focus on the prioritized chems 
dio.aeds <- tissue.bers[variable == 'aed' & min_ber <= 3, 
                    .(dtxsid, chnm, enzyme, variable, compt, value, min_ber)]
dio.aeds[, min_aed_bytarget := min(value), by = .(dtxsid, enzyme)]
dio.aeds <- dio.aeds[order(dtxsid)]

# filter down to the min AED for each target for each chemical
dio.aeds <- dio.aeds[variable == "aed" & value == min_aed_bytarget]
dio.aeds[, variable := paste0(enzyme, ".minAED")]
dio.aeds[, enzyme := NULL]

# how many priority chems do we get from our workflow?
dio.aeds[, length(unique(dtxsid))]
#> [1] 39

priority.chems <- dio.aeds[, unique(dtxsid)]

# load invitrodb v3.5 data 
load("./data/invitrodb_v3_5_snapshot.RData", verbose = TRUE)

# filter down to prioritized chems with BER < 3 
mc5.dio.priority <- mc5[dsstox_substance_id %in% priority.chems][use.me == 1]

# take the 5%-tile AC50 for each chem 
mc5.dio.priority[hitc == 1, modl_ga.5percentile := quantile(modl_ga, c(0.05), type = 2), 
                 by = dsstox_substance_id]
toxcast.pods <- mc5.dio.priority[hitc == 1, .(dsstox_substance_id, chnm, aeid, aenm,
                                              modl,modl_ga, modl_ga.5percentile)]
toxcast.pods <- toxcast.pods[order(dsstox_substance_id)]
toxcast.pods.httk <- unique(toxcast.pods[, .(dsstox_substance_id, chnm, 
                                             ac50_5percentile = modl_ga.5percentile)])
toxcast.pods.httk[, lowerbd.ac50_uM := 10^(ac50_5percentile)]

# convert to AEDs using httk's 3compss model 
library(httk)
load_dawson2021() # most data for ToxCast chems
load_sipes2017() # most data for pharma compounds 
load_pradeep2020() # best ML method for in silico prediction 

for (i in 1:nrow(toxcast.pods.httk)) {
  cat('Calculating for chemical: ', toxcast.pods.httk$chnm[i], '\n')
  set.seed(123)
  out <- calc_mc_css(dtxsid = toxcast.pods.httk$dsstox_substance_id[i], 
                     species = 'Human', 
                     daily.dose = 1, 
                     which.quantile = c(0.5, 0.95), 
                     model = '3compartmentss', 
                     output.units = 'uM', 
                     parameterize.arg.list=list(physchem.exclude=FALSE))
  
  toxcast.pods.httk$mc.css50[i] <- out["50%"]
  toxcast.pods.httk$mc.css95[i] <- out["95%"]
  
}

toxcast.pods.httk[, `:=` (aed50 = lowerbd.ac50_uM/mc.css50, 
                          aed95 = lowerbd.ac50_uM/mc.css95)]
setnames(toxcast.pods.httk, old = "dsstox_substance_id", new = "dtxsid")

# melt the data for ggplot
toxcast.pods.m <- melt.data.table(toxcast.pods.httk, 
                                  id.vars = c("dtxsid", "chnm"),
                                  measure.vars = c("aed50", "aed95"), 
                                  variable.name = "variable", 
                                  value.name = "value")

# convert back to log10 scale
toxcast.pods.m[, value := log10(value)]
toxcast.pods.m[variable == "aed50", variable := "ToxCast.aed50"]
toxcast.pods.m[variable == "aed95", variable := "ToxCast.aed95"]

pods.data <- rbind(dio.aeds[, -c("min_aed_bytarget")], 
                   toxcast.pods.m, 
                   fill = T) # fill = T will add "NA" to min_ber for ToxCast rows 

# order by ToxCast AED95 i.e. the smallest bioactivity-related POD
chnm.order <- toxcast.pods.httk[order(aed95), chnm]
pods.data$chnm <- factor(pods.data$chnm, levels = chnm.order)

pod.breaks <- c(paste0(targets, ".minAED"), paste0("ToxCast.aed", c("50", "95")))
pod.labels <- c("DIO1 min AED", 
                "DIO2 min AED", 
                "DIO3 min AED",
                "IYD min AED", 
                TeX("ToxCast 5%-ile $POD_{50}$"), 
                TeX("ToxCast 5%-ile $POD_{95}$"))

dxp <- ggplot(data = pods.data[variable != "ToxCast.aed50"], 
       mapping = aes(x = chnm, 
                     y = value)) +
  geom_point(aes(shape = variable, 
                 color = variable, 
                 alpha = variable), 
             size = 3, 
             show.legend = T) +
  scale_alpha_manual(breaks = pod.breaks, 
                     labels = pod.labels, 
                     values = c("DIO1.minAED" = 0.8, 
                                "DIO2.minAED" = 0.8,
                                "DIO3.minAED" = 0.8, 
                                "IYD.minAED" = 0.8, 
                                "ToxCast.aed50" = 1, 
                                "ToxCast.aed95" = 1)) +
  scale_shape_manual(breaks = pod.breaks, 
                     labels = pod.labels,  
                     values = c("DIO1.minAED" = 17, 
                                "DIO2.minAED" = 17, 
                                "DIO3.minAED" = 17, 
                                "IYD.minAED" = 17, 
                                "ToxCast.aed50" = 19, 
                                "ToxCast.aed95" = 19)) +
  scale_color_manual(breaks = pod.breaks, 
                       labels = pod.labels, 
                     values = c("DIO1.minAED" = "#56B4E9", 
                                "DIO2.minAED" = "#009E73", 
                                "DIO3.minAED" = "#E69F00", 
                                "IYD.minAED" = "#F0E442", 
                                "ToxCast.aed50" = "#9E9E9E", 
                                "ToxCast.aed95" = "#9E9E9E")) +
  my_theme +
  theme(legend.position = 'bottom', 
        legend.title = element_blank()) + 
  guides(color = guide_legend(nrow = 3, ncol = 2, byrow = F), 
         shape = guide_legend(nrow = 3, ncol = 2, byrow = F), 
         alpha = guide_legend(nrow = 3, ncol = 2, byrow = F)) +
  labs(x = 'Chemical', y = 'Oral AED (log10 mg/kg-bw/day)') +
  coord_flip()

toxcast.pods.httk[, c("aed50", "aed95") := lapply(.SD, log10), 
                  .SDcols = c("aed50", "aed95")]

bioactive.pods <- dxp + geom_segment(data = toxcast.pods.httk, 
                              aes(x = toxcast.pods.httk$chnm, 
                                  y = toxcast.pods.httk$aed95, 
                                  xend = toxcast.pods.httk$chnm, 
                                  yend = toxcast.pods.httk$aed50), 
                              color = "#9E9E9E", linewidth = 1)

# put it altogether in one figure
pods.fig <- plot_grid(toxval.pp, bioactive.pods, 
          nrow = 1, ncol = 2, 
          align = "h", 
          labels = c("A", "B"), 
          label_size = 20)

ggsave(plot = pods.fig, 
       units = "in", 
       dpi = 300, 
       width = 18.53, height = 8.32, 
       device = "tiff", 
       filename = "./figures/300dpi/bioactivity_pods-v3.tiff")

ggsave(plot = pods.fig, 
       units = "in", 
       dpi = 300, 
       width = 18.53, height = 8.32, 
       device = "png", 
       filename = "./figures/bioactivity_pods-v3.png")

