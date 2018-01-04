library(igraph)
library(tidyverse)
library(jsonlite)
library(pbapply)
library(ggplot2)
library(cowplot)

out.folder <- "data/genie3"
scratch.folder <- "data_tmp/genie3"
plots.folder <- "plots/genie3"

load(file = paste0(out.folder, "/datasets.RData"))
load(file = paste0(out.folder, "/outs_parallel.RData"))
source("genie3_h_methods.R")

organism.map <- c("ecoli"="E. coli", "yeast"="S. cerevisiae")

df <- bind_rows(lapply(outs, function(out) {
  co <- out$comparisons[[2]]
  ra <- co$ranking
  
  data_frame(
    dataset = names(datasets)[out$eval$dataset[[1]]],
    method = names(methods)[out$eval$method[[1]]],
    time = co$individual.times,
    datasetf = organism.map[dataset],
    num.edges = seq_along(time),
    num.nodes = sapply(num.edges, function(i) {
      rai <- ra[seq_len(i),]
      length(unique(c(rai$i, rai$j)))
    }),
    avg.degree = num.edges / num.nodes
  )
}))
dfl <- df %>% group_by(dataset) %>% slice(which(num.edges %in% range(num.edges) | time %in% range(time))) %>% ungroup()
df2 <- df %>% spread(method, time) %>% mutate(ratio = orca / incgraph) %>% na.omit
dfl2 <- df2 %>% group_by(dataset) %>% slice(which(num.edges %in% range(num.edges) | ratio %in% range(ratio))) %>% ungroup()

g0 <- ggplot(df) + 
  geom_point(aes(num.edges/1000, time), size = .1) + 
  facet_wrap(~datasetf, scale = "free", ncol = 2) +
  labs(x = "Number of edges", y = "Execution time (s)")
g1 <- ggplot(df %>% filter(method == "incgraph")) + 
  geom_line(aes(num.edges/1000, num.nodes)) + 
  facet_wrap(~datasetf, scale = "free", ncol = 2) +
  labs(x = "Number of edges", y = "Number of nodes")
g2 <- ggplot(df %>% filter(method == "incgraph")) + 
  geom_line(aes(num.edges/1000, avg.degree)) + 
  facet_wrap(~datasetf, scale = "free", ncol = 2) +
  labs(x = "Number of edges", y = "Average degree")
g3 <- ggplot(df2) + 
  geom_point(aes(num.edges/1000, ratio), size = .1) + 
  facet_wrap(~datasetf, scale = "free", ncol = 2) +
  labs(x = "Number of edges", y = "Speedup") +
  geom_point(aes(num.edges/1000, 1), alpha = 0, dfl)

pdf(paste0(plots.folder, "/timings.pdf"), 5, 10)
cowplot::plot_grid(g0, g1, g2, g3, ncol = 1)
dev.off()


g0 <- ggplot(dfl) + 
  geom_point(aes(num.edges/1000, time), size = .1) + 
  facet_wrap(~datasetf, scale = "free", ncol = 2) +
  labs(x = "Number of edges", y = "Execution time (s)")
g1 <- ggplot(dfl %>% filter(method == "incgraph")) + 
  geom_line(aes(num.edges/1000, num.nodes)) + 
  facet_wrap(~datasetf, scale = "free", ncol = 2) +
  labs(x = "Number of edges", y = "Number of nodes")
g2 <- ggplot(dfl %>% filter(method == "incgraph")) + 
  geom_line(aes(num.edges/1000, avg.degree)) + 
  facet_wrap(~datasetf, scale = "free", ncol = 2) +
  labs(x = "Number of edges", y = "Average degree")
g3 <- ggplot(dfl2) + 
  geom_point(aes(num.edges/1000, ratio), size = .1) + 
  facet_wrap(~datasetf, scale = "free", ncol = 2) +
  labs(x = "Number of edges", y = "Speedup") +
  geom_point(aes(num.edges/1000, 1), alpha = 0, dfl)

pdf(paste0(plots.folder, "/timings_lite1.pdf"), 5, 10)
cowplot::plot_grid(g0, g1, g2, g3, ncol = 1)
dev.off()

system(paste0("convert -density 600 ", plots.folder, "/timings.pdf -quality 100 ", plots.folder, "/timings_lite2.png"))


eval <- bind_rows(lapply(outs, function(o) o$eval))
eval


eval.long <- 
  bind_rows(lapply(outs, function(out) {
    meth <- names(methods)[out$eval$method[[1]]]
    df <- data.frame(
      dataset = names(datasets)[out$eval$dataset[[1]]],
      method = meth,
      out$eval[2,]
    )
    if (meth == "incgraph") {
      df.before <- data.frame(
        dataset = names(datasets)[out$eval$dataset[[1]]],
        method = "original",
        out$eval[1,]
      )
      bind_rows(df, df.before)
    } else {
      df
    }
  })) %>% as_data_frame %>% arrange(dataset, method)

length.df <- bind_rows(lapply(outs, function(out) {
  method <- names(methods)[out$eval$method[[1]]]
  dataset <- names(datasets)[out$eval$dataset[[1]]]
  co <- out$comparisons[[2]]
  num.edges <- co$individual.times %>% length
  ra <- co$ranking[seq_len(num.edges),]
  num.nodes <- length(unique(c(ra$i, ra$j)))
  avg.degree <- num.edges / num.nodes
  data_frame(method, dataset, num.edges, num.nodes, avg.degree)
}))
max.length.df <- length.df %>% group_by(dataset) %>% summarise(length = max(num.edges))
max.length.df <- max.length.df %>% mutate(length = 1000)
combined.rankings <- 
  bind_rows(lapply(outs, function(out) {
    datas <- names(datasets)[out$eval$dataset[[1]]]
    meth <- names(methods)[out$eval$method[[1]]]
    
    df <- data.frame(dataset = datas, method = meth, out$comparisons[[2]]$ranking) %>% mutate(rank = seq_len(n()))
    if (meth == "incgraph") {
      df.before <- data.frame(dataset = datas, method = "original", out$comparisons[[1]]$ranking) %>% mutate(rank = seq_len(n()))
      bind_rows(df, df.before)
    } else {
      df
    }
  })) %>% as_data_frame
rankings.limited <- combined.rankings %>% left_join(max.length.df, by = "dataset") %>% filter(rank <= length)
new.evals <- lapply(seq_len(nrow(eval.long)), function(i) {
  ds <- eval.long$dataset[[i]]
  ms <- eval.long$method[[i]]
  ranking.part <- rankings.limited %>% filter(dataset == ds, method == ms)
  GENIE3::evaluate.ranking.direct(-ranking.part$rank, ranking.part$gold)
})
new.eval.long <- bind_rows(lapply(seq_len(nrow(eval.long)), function(i) {
  ds <- eval.long$dataset[[i]]
  ms <- eval.long$method[[i]]
  data.frame(dataset = ds, method = ms, new.evals[[i]]$au.score)
})) %>% as_data_frame()
new.eval.metrics <- bind_rows(lapply(seq_len(nrow(eval.long)), function(i) {
  ds <- eval.long$dataset[[i]]
  ms <- eval.long$method[[i]]
  data.frame(dataset = ds, method = ms, new.evals[[i]]$metrics)
})) %>% as_data_frame()

new.eval.long

ggplot(new.eval.long) + geom_point(aes(auroc, aupr, colour = method)) + facet_wrap(~dataset, scales = "free")

make.plot.df <- new.eval.long %>%
  filter(method %in% c("original", "incgraph")) %>%
  mutate(
    methodf = c("original"="Original network", "incgraph"="IncGraph")[method], 
    datasetf = c("ecoli"="E. coli", "yeast"="S. cerevisiae")[dataset],
    y = c("original"=1, "incgraph"=0)[method], y1 = y - .45, y2 = y + .45)

g <- ggplot(make.plot.df) +
  geom_rect(aes(xmin = 0, xmax = F1, ymin = y1, ymax = y2), stat = "identity") + 
  geom_text(aes(-.01, y, label = methodf), hjust = 1) +
  geom_text(aes(.4, y, label = round(F1, 3)), hjust = 0) +
  xlim(-.1, .43) +
  theme_nothing() +
  facet_wrap(~datasetf, scales = "free", ncol = 1)

pdf(paste0(plots.folder, "/f1score.pdf"), 6, 2.5)
g
dev.off()

em2 <- new.eval.metrics %>% group_by(dataset, method) %>% arrange(desc(cutoff)) %>% slice(seq_len(1000)) %>% ungroup()
g1 <- ggplot(em2, aes(fpr, rec, colour=method)) + geom_step() + facet_wrap(~dataset, scales = "free", ncol = 1)
g2 <- ggplot(em2, aes(rec, prec, colour=method)) + geom_step() + facet_wrap(~dataset, scales = "free", ncol = 1)
plot_grid(g1, g2, ncol = 2)



num.links <- 100
net1 <- outs[[2]]$comparisons[[1]]$ranking %>% head(num.links)
net2 <- outs[[2]]$comparisons[[2]]$ranking %>% head(num.links)
netc <- bind_rows(net1, net2)
grc <- igraph::graph.data.frame(netc, directed = F)
verts <- igraph::as_data_frame(grc, "vertices")
gr1 <- igraph::graph.data.frame(net1, directed = F, vertices = verts)
gr2 <- igraph::graph.data.frame(net2, directed = F, vertices = verts)

lay <- layout_nicely(grc)
isolates1 <- degree(gr1) == 0
isolates2 <- degree(gr2) == 0

pdf(paste0(plots.folder, "/net1.pdf"), 4, 4)
plot(delete.vertices(gr1, which(isolates1)), layout = lay[!isolates1,], vertex.size = 2, vertex.label = NA, vertex.color = "#444444", vertex.frame.color = NA)
dev.off()
pdf(paste0(plots.folder, "/net2.pdf"), 4, 4)
plot(delete.vertices(gr2, which(isolates2)), layout = lay[!isolates2,], vertex.size = 2, vertex.label = NA, vertex.color = "#444444", vertex.frame.color = NA)
dev.off()