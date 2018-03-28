library(incgraph)
library(tidyverse)
library(cowplot)
library(PRISM)

out.folder <- "data/synthetic"
scratch.folder <- "data_tmp/synthetic"
plots.folder <- "plots/synthetic"

load(paste0(out.folder, "/synth_networks.RData"))
source("scripts/synthetic_h_delta_functions.R")

timeout <- 60

qsub.out <- qsub_lapply(
  X = seq_along(networks),
  qsub_config = override_qsub_config(num_cores = 8, wait = F, remove_tmp_folder = F, stop_on_error = F),
  qsub_environment = c("dplyr", "incgraph"),
  FUN = function(i) {
    net <- networks[[i]]
    prm <- params[i,]
    list2env(prm, globalenv())

    out <- dplyr::bind_rows(lapply(delta.functions, function(dfun) {
      state <- dfun$init(num.nodes, as.matrix(net$network))

      start <- Sys.time()
      time.diff <- as.numeric(Sys.time() - start, units = "secs")
      j <- 0
      while(j < nrow(net$operations) && time.diff <= timeout) {
        j <- j + 1
        state <- dfun$calc.delta(state, net$operations[j,-1], T)$state
        time.diff <- as.numeric(Sys.time() - start, units = "secs")
      }

      data.frame(method = dfun$name, prm, num.executed = j, time.diff)
    }))
    out
  })

save.image(file = paste0(scratch.folder, "/synth_networks_timings.RData"))

# retrieve output
load(file = paste0(scratch.folder, "/synth_networks_timings.RData"))
outs <- qsub_retrieve(qsub.out)
save.image(file = paste0(out.folder, "/synth_networks_timings.RData"))


outs <- pbapply::pblapply(
  X = seq_along(networks),
  # cl = 8,
  FUN = function(i) {
    net <- networks[[i]]
    prm <- params[i,]
    list2env(prm, globalenv())
    
    out <- dplyr::bind_rows(lapply(delta.functions, function(dfun) {
      state <- dfun$init(num.nodes, as.matrix(net$network))
      
      start <- Sys.time()
      time.diff <- as.numeric(Sys.time() - start, units = "secs")
      j <- 0
      while(j < nrow(net$operations) && time.diff <= timeout) {
        j <- j + 1
        state <- dfun$calc.delta(state, net$operations[j,-1], T)$state
        time.diff <- as.numeric(Sys.time() - start, units = "secs")
      }
      
      data.frame(method = dfun$name, prm, num.executed = j, time.diff)
    }))
    out
  })
save.image(file = paste0(out.folder, "/synth_networks_timings_laptop.RData"))