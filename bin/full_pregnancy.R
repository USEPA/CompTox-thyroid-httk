
library(httk)
library(data.table)
library(dplyr)
library(tidyverse)
library(reshape2)
library(RColorBrewer)

# typical compartments for individuals 
maternal_compts <- c('gutlumen', 'gut', 'liver', 'kidney', 'lung', 'ven', 'art', 
                     'adipose','thyroid', 'rest')
maternal_states <- paste0('A', maternal_compts)
maternal_states <- c(maternal_states, 
                     'Atubules', 'Ametabolized', 'AUC')
maternal_concs <- paste0("C", maternal_compts[! (maternal_compts %in% c("gutlumen"))])

fetal_compts <- c(maternal_compts[!( maternal_compts %in% c('adipose', 'gutlumen') )], 
                  "brain")
fetal_states <- c(paste0('Af', fetal_compts), 'fAUC')
fetal_concs <- paste0("Cf", fetal_compts)

# all compts by the 2nd/3rd trimesters
mf.states <- c(maternal_states, fetal_states, "Aplacenta")
mf.outputs <- c(mf.states, # all possible outputs of fetal_pbtk
                maternal_concs, 
                fetal_concs, 
                "Cplacenta", 
                "Cplasma", "Aplasma", "Rblood2plasma", 
                "Cfplasma", "Afplasma", "Rfblood2plasma")

# spell out all possible outputs of 1tri_pbtk
firsttri.states <- c(maternal_states, "Aconceptus")
firsttri.outputs <- c(firsttri.states, 
                      maternal_concs, "Cconceptus", 
                      "Cplasma", "Aplasma", "Rblood2plasma")

full_pregnancy <- function(dtxsid, track.vars = NULL, plt = FALSE,
                           return.units = "amt",
                           time.course = seq(0,40*7,1/24),
                           ...) {
  cat("Solving for chemical: ", dtxsid, "\n")
  
  # split time.course to the appropriate domain for each model
  t1 = time.course[time.course <= 13*7]
  t2 = time.course[time.course > 13*7]
  
  # we add day 91 to both time series to stitch the models together
  t1 = sort(unique(c(0,t1, 13*7))) # ode solver needs the first value of "times" to be the initial time (T0=0)
  t2 = c(13*7,t2) # initial time T0 = 91
  
  # track chemical amounts (i.e. state vars of each model) 
  firsttri.out <- solve_1tri_pbtk(dtxsid = dtxsid, 
                                  times = t1, 
                                  dose = 0,
                                  # rtol = 1e-7, # default: rtol = 1e-8, these specs are for Ergocalciferol
                                  # atol = 1e-11, # default: atol = 1e-12
                                  monitor.vars = firsttri.outputs, 
                                  ...)
  
  # initialize vector for "initial.values" input to fetal_pbtk 
  initial.dat <- setNames(rep(0, length(mf.states)), mf.states)
  
  # populate end values for maternal compts 
  ind <- which(firsttri.out[, 'time'] == 13*7)
  initial.dat[maternal_states] = firsttri.out[ind, maternal_states]
  
  # compute initial amts for fetal compts (and placenta)
  missing.amts <- c(paste0("Af", fetal_compts),
                    "Aplacenta")
  
  # partition out chemical into fetal compts + placenta assuming a well-mixed conceptus
  missing.vols <- str_replace(missing.amts, "^A", "V")
  
  # get volumes from C model implementation
  vols.out <- solve_fetal_pbtk(dtxsid = dtxsid,
                               dose = 0, 
                               times = c(13*7, 13*7+1), # times needs to contain at least 2 values
                               monitor.vars = c(missing.vols, "fhematocrit")
  )
  
  # # check that vols.out first row is for day 91 
  # Af = firsttri.out[ind, "Cconceptus"]*vols.out[1, missing.vols]
  # names(Af) <- missing.amts 
  # initial.dat[names(Af)] = Af 
  
  fetal.parms <- parameterize_fetal_pbtk(dtxsid = dtxsid)
  
  # get fetal tissue partition coefficients 
  fetal.pcs <- c(fetal.parms[substr(names(fetal.parms),1,2) == 'Kf'], 
                 fetal.parms["Kplacenta2pu"])
  
  fetal.tissues <- fetal_compts[! (fetal_compts %in% c("ven", "art"))]
  
  # reorder fetal volumes to match order of fetal.pcs 
  Vf = vols.out[1, c(paste0("Vf", fetal.tissues), "Vplacenta")]
  Kf = unlist(fetal.pcs[c(paste0("Kf", fetal.tissues,"2pu"), "Kplacenta2pu")])
  
  KVtotal = sum(Kf*Vf)
  
  # add art,ven blood components
  bV = vols.out[1, c("Vfart", "Vfven")]
  KVtotal = KVtotal + vols.out[[1, "fhematocrit"]]*fetal.pcs[["Kfrbc2pu"]]*sum(bV) #RBCs
  KVtotal = KVtotal + + (1-vols.out[[1, "fhematocrit"]])/fetal.parms$Fraction_unbound_plasma_fetus*sum(bV) #plasma 
  
  # check that vols.out first row is for day 91 
  Af = firsttri.out[ind, "Aconceptus"]*Vf*Kf/KVtotal
  names(Af) <- sub("^V", "A", names(Af))
  initial.dat[names(Af)] = Af 
  
  # compute amounts for Afart, Afven based on avg partition coefficient of RBCs and plasma 
  initial.dat["Afart"] = firsttri.out[ind, "Aconceptus"]*vols.out[1, "Vfart"]*vols.out[[1, "fhematocrit"]]*fetal.pcs[["Kfrbc2pu"]] 
  initial.dat["Afart"] = initial.dat["Afart"] + firsttri.out[ind, "Aconceptus"]*vols.out[1, "Vfart"]*(1-vols.out[[1, "fhematocrit"]])/fetal.parms$Fraction_unbound_plasma_fetus 
  initial.dat["Afart"] = initial.dat["Afart"]/KVtotal
  
  initial.dat["Afven"] = firsttri.out[ind, "Aconceptus"]*vols.out[1, "Vfven"]*vols.out[[1, "fhematocrit"]]*fetal.pcs[["Kfrbc2pu"]] 
  initial.dat["Afven"] = initial.dat["Afven"] + firsttri.out[ind, "Aconceptus"]*vols.out[1, "Vfven"]*(1-vols.out[[1, "fhematocrit"]])/fetal.parms$Fraction_unbound_plasma_fetus 
  initial.dat["Afven"] = initial.dat["Afven"]/KVtotal
  
  # assuming no placental barrier 
  initial.dat["fAUC"] = initial.dat["AUC"]
  
  # modified input to fetal_pbtk
  mod.fetal.out <- solve_fetal_pbtk(dtxsid = dtxsid, 
                                    times = t2, 
                                    dose = 0, 
                                    monitor.vars = mf.outputs, 
                                    initial.values = initial.dat, 
                                    ...)
  
  # get full solution by concatenating 2 outputs
  full_sol <- bind_rows(data.frame(firsttri.out), data.frame(mod.fetal.out))
  
  # add initial amts approximated at day 91 from Cconceptus 
  # always return the solution with calculations and fetal_pbtk solution for day 91
  full_sol[ind, c(missing.amts, "fAUC")] = initial.dat[c(missing.amts, "fAUC")]
  
  # plot all available states (the amts)
  if (plt == TRUE) {
    
    if (return.units == "amt") {
      # cols <- c(maternal_states[maternal_states != "AUC"], 
      #           fetal_states[fetal_states != "fAUC"], 
      #           "Aconceptus", "Aplacenta")
      
      cols <- c(maternal_states, fetal_states, 
                "Aconceptus", "Aplacenta")
      
      out <- full_sol[, c("time", cols)]
      
      # subset down to requested times 
      out <- out[which(out$time %in% time.course), ]
      
      setDT(out)
      out[, Mtotal := Agutlumen + Agut + Aliver + Akidney + Alung + Aven + Aart + Aadipose + Athyroid + Arest]
      out[, ftotal := Afgut + Afliver + Afkidney + Aflung + Afven + Afart + Afthyroid + Afrest + Afbrain]
      
      # melt the data to ggplot 
      out.m <- reshape2::melt(out, id.vars = c("time"), 
                              variable.name = 'tissue', 
                              value.name = return.units)
      setDT(out.m)
      out.m[, model := '1st trimester']
      out.m[time > 91, model := '2nd-3rd trimester']
      # out.m[, body := "maternal"]
      # out.m[tissue %in% colnames(out)[grep("^Af", colnames(out))], body := "fetal"]
      # out.m[tissue %in% c("Aplacenta", "Aconceptus"), body := "fetal"]
      # out.m[, tissue := sub("^A[f]*", "", tissue)]
      # out.m[tissue == "ftotal", body := "fetal"]
      # out.m[tissue == "ftotal", tissue := "total"]
      # out.m[tissue == "Mtotal", tissue := "total"]
      
      # plot all compartments on a graph faceted by mother/fetus 
      # all.tissues <- c(union(maternal_compts, fetal_compts), "conceptus", "placenta") 
      # set3.colors <- brewer.pal(12, 'Set3') # max num of colors for Set 3 is 12 
      # tissue.colors <- setNames(set3.colors, all.tissues[all.tissues != "gutlumen"])
      # tissue.colors <- c(tissue.colors, "#7FC97F", "#E41A1C", "#FF7F00", "#000000")
      # names(tissue.colors)[13:length(tissue.colors)] = c("gutlumen", "metabolized", "tubules", "total")
      
      p <- ggplot(out.m[!is.na(get(return.units))], aes(x = time, y = log10(get(return.units)))) + 
        geom_point(aes(color = model)) +
        scale_color_manual(values = c('1st trimester' = 'red', '2nd-3rd trimester' = 'black')) +
        facet_wrap(~tissue) + 
        theme_bw() + 
        # theme(title = element_text(size = 16),
        #       axis.title.x = element_text(size = 16),
        #       axis.title.y = element_text(size = 16), 
        #       axis.text.y = element_text(size = 13), 
        #       legend.text = element_text(size = 13), 
        #       legend.title = element_text(size = 13), 
        #       strip.text = element_text(size = 13)) + 
        labs(x = 'time (days)', y = 'Amount (log10 umol)', 
             title = 'Chemical Amounts in Compartments of full gestational model') +
        # xlim(0,300) +
        guides(colour = guide_legend(override.aes = list(size = 5)))
      
      print(p)
      
    } else if (return.units == "conc") {
      cols <- c(maternal_concs, fetal_concs, "Cconceptus", "Cplacenta")
      
      out <- full_sol[, c("time", cols)]
      
      # subset down to requested times 
      out <- out[which(out$time %in% time.course), ]
      
      setDT(out)
      
      # melt the data to ggplot 
      out.m <- reshape2::melt(out, id.vars = c("time"), 
                              variable.name = 'tissue', 
                              value.name = return.units)
      setDT(out.m)
      # out.m[, body := "maternal"]
      # out.m[tissue %in% colnames(out)[grep("^Cf", colnames(out))], body := "fetal"]
      # out.m[tissue %in% c("Cconceptus", "Cplacenta"), body := "conceptus"]
      # out.m[, tissue := sub("^C[f]*", "", tissue)]
      # 
      # # plot all compartments on a graph faceted by mother/fetus 
      # all.tissues <- c(union(maternal_compts, fetal_compts), "conceptus", "placenta") 
      # set3.colors <- brewer.pal(12, 'Set3') # max num of colors for Set 3 is 12 
      # tissue.colors <- setNames(set3.colors, all.tissues[all.tissues != "gutlumen"])
      # tissue.colors <- c(tissue.colors, "#7FC97F")
      # names(tissue.colors)[13] = c("gutlumen")
      
      p <- ggplot(out.m[!is.na(get(return.units))], 
                  aes(x = time, y = log10(get(return.units)))) + 
        geom_point(aes(color = tissue)) +
        # scale_color_manual(values = tissue.colors) +
        facet_wrap(~body) + 
        theme_bw() + 
        # theme(title = element_text(size = 16),
        #       axis.title.x = element_text(size = 16),
        #       axis.title.y = element_text(size = 16), 
        #       axis.text.y = element_text(size = 13), 
        #       legend.text = element_text(size = 13), 
        #       legend.title = element_text(size = 13), 
        #       strip.text = element_text(size = 13)) + 
        labs(x = 'time (days)', y = 'Concentration (log10 uM)', 
             title = 'Chemical Concentration in Compartments of full gestational model') +
        # xlim(0,300) +
        guides(colour = guide_legend(override.aes = list(size = 5)))
      
      print(p)
      
    } else{
      stop("Acceptable values for return.units to plot is 'conc' or 'amt.'")
    }
  }  
  
  # The monitored variables can be altered by the user 
  if (is.null(track.vars)) {
    
    # have the default output columns be selected concs 
    default.track.vars <- c("Agutlumen",maternal_concs,
                              "Aconceptus", "Cconceptus",
                              "Cplasma",
                              "Atubules","Ametabolized","Rblood2plasma",
                              "AUC","fAUC", 
                              "Aplacenta", "Cplacenta",
                              fetal_concs, 
                              "Cfplasma","Rfblood2plasma")
    return(full_sol[, c("time", default.track.vars)])
  }
  else {
    
    # however, always include the compartment that receives the dose 
    return(full_sol[, unique(c("time", "Agutlumen", track.vars))])
  }
  
}