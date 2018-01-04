library(incgraph)
library(tidyverse)
library(cowplot)
library(PRISM)

out.folder <- "data/synthetic"
scratch.folder <- "data_tmp/synthetic"
plots.folder <- "plots/synthetic"
dir.create(out.folder, recursive = T)
dir.create(scratch.folder, recursive = T)
dir.create(plots.folder, recursive = T)

params <- expand.grid(
  model = c("BA", "ER", "GEO"), 
  num.nodes = round(2^seq(0, 4, by = .5) * 1000), 
  avg.degree = 2^seq(1, 6),
  num.operations = 10000,
  repeats = 1:10,
  stringsAsFactors = F
) %>% mutate(
  num.edges = num.nodes * avg.degree / 2
)

qsub.out <- qsub.lapply(
  X = seq_len(nrow(params)), 
  qsub.config = qsub.configuration(wait = F),
  FUN = function(i) {
  model <- params$model[[i]]
  num.nodes <- params$num.nodes[[i]]
  avg.degree <- params$avg.degree[[i]]
  num.edges <- params$num.edges[[i]]
  num.operations <- params$num.operations[[i]]
  
  incgraph::generate.dynamic.network(model = model, num.nodes, num.edges, num.operations, trace = F)
})

save(params, qsub.out, file = paste0(scratch.folder, "/synth_networks.RData"))
load(file = paste0(scratch.folder, "/synth_networks.RData"))

qsub.out$remove.tmpdirs <- F
networks <- qsub.retrieve(qsub.out)

save(params, networks, file = paste0(out.folder, "/synth_networks.RData"))
