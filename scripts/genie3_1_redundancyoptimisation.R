library(igraph)
library(GENIE3)
library(tidyverse)
library(incgraph)
library(PRISM)
library(jsonlite)
library(pbapply)
library(parallel)
library(ggplot2)
library(cowplot)

out.folder <- "data/genie3"
scratch.folder <- "data_tmp/genie3"

source("scripts/genie3_h_methods.R")

# define redundancy vectors
redundant.edges <- c(0, 0, 1, 0, 0, 1, 1, 2, 3, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 6)
names(redundant.edges) <- paste0("G", seq_along(redundant.edges)-1)
redundant.edges2 <- c(0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0,
                      1, 1, 1, 0, 2, 2, 2, 2, 2, 0, 1, 1, 2, 2, 2, 1, 1, 2, 3, 3, 0, 3, 2, 2, 3, 3, 1, 3, 3, 1, 3, 3, 4, 3, 4, 5, 6)
names(redundant.edges2) <- paste0("O", seq_along(redundant.edges2)-1)


load(file = paste0(out.folder, "/datasets.RData"))

for (i in seq_along(datasets)) {
  datasets[[i]]$genes <- colnames(datasets[[i]]$expression)
  datasets[[i]]$undirected.links <- datasets[[i]]$links %>%
    mutate(i2 = pmin(i, j), j2 = pmax(i, j)) %>%
    group_by(i2, j2) %>%
    arrange(desc(value)) %>%
    summarise(value = value[[1]], gold = max(gold)) %>%
    select(i = i2, j = j2, value, gold) %>%
    ungroup %>%
    arrange(desc(value))
}

#const.param <- list(passes.max = 5, window.size = 100, quan.pct = .7, quan.mult = 1)
const.param <- list(passes.max = 5, window.size = 100, quan.pct = .7, quan.mult = 1, max.time = 24*3600, mc.cores = 24)

params <- do.call(expand.grid, c(const.param, list(dataset = seq_along(datasets), method = seq_along(methods))))

qsub <- qsub_lapply(
  X = seq_len(nrow(params)),
  qsub_config = override_qsub_config(name = "IncGraph", num_cores = 24, wait = F, remove_tmp_folder = F, stop_on_error = F, max_wall_time = NULL),
  qsub_packages = c("GENIE3", "dplyr", "tidyr", "incgraph"),
  FUN = function(i) {
    param <- params[i,]
    method.index <- param$method
    dataset.index <- param$dataset
    passes.max <- param$passes.max
    window.size <- param$window.size
    quan.pct <- param$quan.pct
    quan.mult <- param$quan.mult
    max.time <- param$max.time
    mc.cores <- param$mc.cores

    method <- methods[[method.index]]
    dataset <- datasets[[dataset.index]]
    undirected.links <- dataset$undirected.links
    genes <- dataset$genes
    dataset.name <- dataset$name

    #' ## GROWING NETWORK, USE FIRST EDGE WITH LOW ENOUGH REDUNDANCY
    undirected.links2 <- undirected.links
    comparisons <- list()

    eval.first <- GENIE3::evaluate_ranking_direct(undirected.links$value, undirected.links$gold, sum(undirected.links$gold), nrow(undirected.links))
    comparisons[[1]] <- list(ranking = undirected.links, evaluation = eval.first, total.time = 0, individual.times = 0)

    time.start1 <- Sys.time()
    time.diff.tmp <- 0

    state <- incgraph::new.incgraph.network(length(genes), links = NULL)

    number.i <- 1
    todos <- seq(1, window.size)

    output <- list()

    while (time.diff.tmp < max.time) {
      if (number.i %% 100 == 0) cat(method$name, ", dataset: ", dataset.index, ", edge: ", number.i, ", time elapsed: ", round(time.diff.tmp), " / ", max.time, "\n", sep="")

      time.start2 <- Sys.time()

      if (number.i != 1) {
        start.out <- method$start(state)
      } else {
        start.out <- matrix(0, nrow = length(genes), ncol = 73)
        colnames(start.out) <- paste0("O", seq_len(ncol(start.out))-1)
      }

      redundancies <- unlist(parallel::mclapply(todos, mc.cores = mc.cores, function(number.j) {
        ei <- undirected.links2$i[[number.j]]
        ej <- undirected.links2$j[[number.j]]
        delta.gr <- method$delta(state, ei, ej, start.out)
        sum(delta.gr[ei,] * redundant.edges2 + delta.gr[ej,] * redundant.edges2)
      }))
      if (!is.finite(quan.mult)) {
        quan.redundancy <- -Inf
      } else {
        quan.redundancy <- quantile(redundancies, quan.pct) * quan.mult
      }

      do.includes <- redundancies <= quan.redundancy
      if (any(do.includes)) {
        first.to.include <- first(which(do.includes))
      } else {
        first.to.include <- which.min(redundancies)
      }
      this.redundancy <- redundancies[[first.to.include]]

      number.k <- todos[[first.to.include]]
      todos <- c(todos[-first.to.include], number.i + window.size)

      ei <- undirected.links2$i[[number.k]]
      ej <- undirected.links2$j[[number.k]]

      gold <- undirected.links2$gold[[number.k]]
      incgraph::flip(state, ei, ej)

      time.stop2 <- Sys.time()
      time.diff2 <- as.numeric(time.stop2 - time.start2, units = "secs")
      time.diff.tmp <- as.numeric(time.stop2 - time.start1, units = "secs")

      output[[number.i]] <- data.frame(i = number.i, method = method$name, link = number.k, gold, this.redundancy, quan.redundancy, time = time.diff2)

      number.i <- number.i + 1
    }

    output.df <- bind_rows(output) %>% as_data_frame %>% mutate(gold.str = ifelse(gold, "gold", "notgold"))


    ix <- c(output.df$link, setdiff(seq_len(nrow(undirected.links2)), output.df$link))
    undirected.links2 <- undirected.links2[ix,]

    eval <- GENIE3::evaluate_ranking_direct(-seq_len(nrow(undirected.links2)), undirected.links2$gold, sum(undirected.links2$gold), nrow(undirected.links2))

    time.stop1 <- Sys.time()
    time.diff1 <- as.numeric(time.stop1 - time.start1, units = "secs")

    comparisons[[2]] <- list(ranking = undirected.links2, evaluation = eval, total.time = time.diff1, individual.times = output.df$time)

    eval <- bind_rows(lapply(seq_along(comparisons), function(j) {
      data.frame(i, param, dataset.name, iteration = j - 1, comparisons[[j]]$evaluation$au.score, total.time = comparisons[[j]]$total.time, mean.ind.time = mean(comparisons[[j]]$individual.times))
    }))
    list(eval = eval, comparisons = comparisons)
  })
save(qsub, file=paste0(scratch.folder, "/qsub_parallel.RData"))

# retrieve output
load(file=paste0(scratch.folder, "/qsub_parallel.RData"))
outs <- qsub_retrieve(qsub)
save(outs, file=paste0(out.folder, "/outs_parallel.RData"))
