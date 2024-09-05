# ==============================================================================
# Script creates final table of chemicals tested in DIO/IYD assays with relevant
# fields for better interpretation i.e. sc/mc hitcalls, efficacy/potencies, 
# selectivity, class identification by Richard, and flags for assay interference 
# (namely, if is positive for chain_alkaneLinear{8,10,12,14} ToxPrint or highly 
# promiscuous across ToxCast for colorful substances). 
# @author: Kimberly Truong
# created: 11/16/2023
# updated: 8/22/2024
# ==============================================================================

# Set up Env and Load Data -----------------------------------------------------
rm(list=ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(readxl) # to read in Excel files with multiple sheets
library(openxlsx)

# create a `not in` operator
`%notin%` <- Negate(`%in%`)

# load all data files 
setwd('C:/Users/ktruong/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/thyroid_prioritization/doc/comptox-thyroid-httk/')
load('./data/invitrodb_v3_5_thyroid_data.RData', verbose = TRUE)
load('./data/invitrodb_v3_5_deiod_Enrichment.RData', verbose = TRUE)

# load in Richard's latest classifications 
new.classified <- read_excel(file.path("data", "chem_classifications.xlsx")) %>% 
  as.data.table()

# combine QACs from "Expert" to "ClassyFire" group
new.classified[class_type == "Expert" & chosen_class == "Quaternary ammonium salts", class_type := "ClassyFire"]

# subset down to DIO/IYD test chems 
classified_dt <- new.classified[dtxsid %in% rownames(M.hitc)]
setnames(classified_dt, "dtxsid", "dsstox_substance_id")

interference.tbl <- read_excel('./data/olker19_assay_interference_suppTable5.xlsx')
setnames(interference.tbl, 
         old = c('CASRN', 'Reason Excluded'), 
         new = c('casn', 'flag'))

surfactants <- read_excel('./data/Chemical_List_ALLSURFACTANTS-2024-09-05.xlsx')
solvents <- read_excel('./data/Chemical_List_PARISIII-2024-09-05.xlsx')

# DIO/IYD assays
deiod.meta <- ace.thyroid[thyroid_related_target == 'Deiodinases', 
                          .(aeid, aenm = assay_component_endpoint_name)]

aenm.abbrevs <- c('DIO1', 'DIO2', 'DIO3', 'IYD')

# ignore amphibian assays 
h.deiod.aeids <- deiod.meta[grepl("\\w+h(DIO|IYD)\\w+", deiod.meta$aenm), aeid]
h.deiod.aenm <- deiod.meta[grepl("\\w+h(DIO|IYD)\\w+", deiod.meta$aenm), aenm]
names(h.deiod.aeids) <- aenm.abbrevs
names(h.deiod.aenm) <- aenm.abbrevs

# get efficacy/potency fields for each
sc.deiod <- sc2[aeid %in% h.deiod.aeids & !is.na(dsstox_substance_id), 
                .(casn, chnm, hitc, max_med), keyby = .(aeid, dsstox_substance_id)]
mc.deiod <- mc5[aeid %in% h.deiod.aeids & !is.na(dsstox_substance_id), 
                      .(casn, chnm, hitc, modl_acc, modl_ga), keyby = .(aeid, dsstox_substance_id)]

# calculate selectivity for mc data 
mc.acc <- mc.deiod[, .(hitc, modl_acc), by = .(dsstox_substance_id, aeid)]
mc.acc[is.na(modl_acc), modl_acc := 5]
mc.acc$cytomed <- cytotox$cytotox_median_log[match(mc.acc$dsstox_substance_id, cytotox$dsstox_substance_id)]

# check if there are rows where cytotox point is missing 
# these are precisely the pos controls in DIO2/3/IYD assays 
mc.acc[is.na(cytomed)]
mc.acc[is.na(cytomed), cytomed := 3] # fill in with the default cytotox point 

mc.acc[, selectivity := cytomed - pmin(modl_acc, 5)]

# add selectivity as a column for mc table 
mc.deiod2 <- merge.data.table(mc.deiod, 
                              mc.acc[, .(dsstox_substance_id, aeid, selectivity)], 
                              by = c('aeid', 'dsstox_substance_id'), 
                              all.x = TRUE)

# merge sc/mc data 
scmc.dat <- merge.data.table(sc.deiod, mc.deiod2, 
                                by = c('aeid', 'dsstox_substance_id', 'casn', 'chnm'), 
                                all.x = TRUE)
setnames(scmc.dat, 
         old = c('hitc.x', 'hitc.y'), 
         new = c('sc.hitc', 'mc.hitc'))

# add class information for qSAR
scmc.sar <- merge.data.table(scmc.dat, 
                             classified_dt[, .(dsstox_substance_id, class_type, chosen_class)], 
                             by = c('dsstox_substance_id'), all.x = TRUE)

# flag chems that were experimentally verified to have interfered with assay 
all_dat <- merge.data.table(scmc.sar,
                            interference.tbl[, c('casn', 'flag')], 
                            by = c('casn'), 
                            all.x = TRUE)

setnames(all_dat, "flag", "flags")

# flag solvents
all_dat[dsstox_substance_id %in% solvents$DTXSID,
        flags := ifelse(is.na(flags), 'solvent', paste(flags, 'solvent', sep = '; '))]

# flag surfactants
all_dat[dsstox_substance_id %in% surfactants$DTXSID,
        flags := ifelse(is.na(flags), 'surfactant', paste(flags, 'surfactant', sep = '; '))]

# flag quaternary ammonium salts 
all_dat[chosen_class == "Quaternary ammonium salts",
        flags := "QAC"]

# how many QACs are there?
all_dat[flags == "QAC", length(unique(dsstox_substance_id))]
#>[1] 10

# label all long-chain Cs containing the enriched CTs
all_dat[dsstox_substance_id %in% longchainc.chems, flags := ifelse(is.na(flags), 
                                                             "contain enriched CT", 
                                                             paste(flags, "contain enriched CT", sep = "; "))]

# add extra annotation for long-chain Cs positive for all DIO targets 
all_dat[dsstox_substance_id %in% dio.panactive, flags := ifelse(is.na(flags), 
                                                                "DIO panactive", 
                                                                paste(flags, "DIO panactive", sep = "; "))]

# label colorful substances that are promiscuous across ToxCast (hit rate > 30%)
all_dat[dsstox_substance_id %in% colors.panactive, flags := ifelse(is.na(flags), 
                                                                   "promiscuous across db", 
                                                                   paste(flags, "promiscuous across db", sep = "; "))]

# Add flag for chems that are ill-defined in terms of structure (i.e. have no ToxPrint)
toxprints <- fread("./data/DIO_IYD_testlib_invitrodb_v3_5_input_toxprints.csv")
setnames(toxprints, "DTXSID", "dsstox_substance_id")
toxprints[!is.na(cid), has_smile := 1]
toxprints[cid == "", has_smile := 0]

all_dat$has_smile <- toxprints$has_smile[match(all_dat$dsstox_substance_id, toxprints$dsstox_substance_id)]

# map aeids to enzyme for easy reference 
for (i in 1:4) {
  all_dat[aeid == h.deiod.aeids[i], deiodinase_type := aenm.abbrevs[i]]
}

# reorganize columns in a sensical order
all_dat <- all_dat[, c("aeid", "deiodinase_type", "casn", "dsstox_substance_id", "chnm", "has_smile",  
                       "sc.hitc", "max_med", 
                       "mc.hitc", "modl_acc", "modl_ga", "selectivity", 
                       "class_type", "chosen_class", 
                       "flags")]

# filter down to selective chems for each target (with score > 0.3)
sel <- all_dat[mc.hitc == 1 & selectivity > 0.3 & is.na(flags)][
  order(-selectivity), .SD, keyby = deiodinase_type] 

sel[, .N, by = deiodinase_type] 
#>   deiodinase_type  N
#>1:            DIO1 52
#>2:            DIO2 59
#>3:            DIO3 88
#>4:             IYD 15

length(unique(sel$dsstox_substance_id))
#> [1] 131

# add flag for cytotoxicity using selectivity 
all_dat[selectivity < 0.3, flags := ifelse(is.na(flags), "nonselective", 
                                           paste(flags, "nonselective", sep = "; "))]

# compare above to pre-filtering out likely false positives
# flags (alone) clean up 58%-87% of the mc positives! 
all_dat[mc.hitc == 1, .N, keyby = deiodinase_type]
#>  deiodinase_type   N
#>1:            DIO1 154
#>2:            DIO2 209
#>3:            DIO3 209
#>4:             IYD 119

1-(52/154)
#> [1] 0.6623377
1-(59/209)
#> [1] 0.7177033
1-(88/209)
#> [1] 0.5789474
1-(15/119)
#> [1] 0.8739496

# 3-6% were mixtures/UVCBs
sel2 <- sel[has_smile == 1][
  order(-selectivity), .SD, keyby = deiodinase_type] 

sel2[, .N, by = deiodinase_type]
#>   deiodinase_type  N
#>1:            DIO1 48
#>2:            DIO2 52
#>3:            DIO3 76
#>4:             IYD 12

sel2[, length(unique(dsstox_substance_id))]
#> [1] 117

# AC50 table by chem x aeid  
sel2[, ac50_uM := 10^modl_ga]
ac50.dt <- dcast(sel2, dsstox_substance_id + chnm ~ deiodinase_type, value.var = 'ac50_uM')
setDT(ac50.dt)

# determine minimum bioactive threshold for each chem 
ac50.dt[, min_ac50 := pmin(DIO1, DIO2, DIO3, IYD, na.rm = TRUE)]

# create RData for IVIVE
save(all_dat,
     sel,
     ac50.dt,
     file = './data/invitrodb_v3_5_deiod_filtered_httk.RData')

# save table for Supplemental File 4 
write.xlsx(all_dat, "./tables/Supp4_DIO-IYD_invitrodbv3_5_processed.xlsx")

