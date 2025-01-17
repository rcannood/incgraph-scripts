library(incgraph)
library(tidyverse)
library(cowplot)
library(PRISM)

out.folder <- "data/synthetic"
scratch.folder <- "data_tmp/synthetic"
plots.folder <- "plots/synthetic"

load(paste0(out.folder, "/synth_networks.RData"))

load(file = paste0(out.folder, "/synth_networks_timings.RData"))

# investigate h-index
h_indices <- sapply(seq_along(networks), function(i) {
  init.net <- networks[[i]]$network
  
  gr <- igraph::graph_from_data_frame(init.net, directed = FALSE, vertices = seq_len(params$num.nodes[[i]]))
  
  sorted_degree <- sort(igraph::degree(gr), decreasing = T)
  position <- seq_along(sorted_degree)
  
  sum(sorted_degree >= position)
})

outs2 <- seq_along(outs) %>% map_df(function(i) {
  outs[[i]] %>% mutate(network_index = i, h_index = h_indices[[i]], rate = time.diff / num.executed)
})

speedup.df <- outs2 %>% 
  select(-time.diff, -num.executed) %>% 
  spread(method, rate) %>% 
  mutate(
    speedup = orca / incgraph,
    num.nodes.f = factor(num.nodes, levels = sort(unique(num.nodes))),
    num.nodes.l = round(log2(num.nodes/1000), 1),
    num.nodes.fl = factor(paste0("1000 × 2^", num.nodes.l), levels = paste0("1000 × 2^", sort(unique(num.nodes.l)))),
    avg.degree.f = factor(avg.degree, levels = sort(unique(avg.degree))),
    num.operations.f = factor(num.operations, levels = sort(unique(num.operations))),
    num.edges.f = factor(num.edges, levels = sort(unique(num.edges)))
  )

summ.df <- speedup.df %>% 
  group_by(model, num.nodes, avg.degree, num.operations, num.nodes.l, avg.degree.f) %>%
  summarise(speedup = mean(speedup), incgraph = mean(incgraph), orca = mean(orca), h.index = mean(h_index)) %>% 
  ungroup()

ggplot(speedup.df) + geom_density(aes(speedup))

ggplot(speedup.df) +
  geom_tile(aes(avg.degree.f, num.nodes.fl, fill = speedup)) +
  facet_wrap(~model) + 
  scale_fill_distiller(palette = "RdBu", trans = "log", breaks = 10^seq(1,6)) +
  labs(x = "Average Degree", y = "Number of nodes", fill = "Speedup") +
  theme(legend.position = "bottom")


ggplot(summ.df) +
  geom_line(aes(num.nodes.l, incgraph, colour = avg.degree.f), alpha = .5) +
  geom_point(aes(num.nodes.l, incgraph, colour = avg.degree.f)) +
  facet_wrap(~model) + 
  scale_y_log10(breaks = 10^seq(-10,10)) + 
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Number of nodes", y = "IncGraph Time", colour = "Placeholder")


ggplot(summ.df) +
  geom_line(aes(num.nodes.l, orca, colour = avg.degree.f), alpha = .5) +
  geom_point(aes(num.nodes.l, orca, colour = avg.degree.f)) +
  facet_wrap(~model) + 
  scale_y_log10(breaks = 10^seq(-10,10)) + 
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Number of nodes", y = "Orca Time", colour = "Placeholder")

ggplot(summ.df) +
  geom_path(aes(orca, incgraph, colour = avg.degree.f), alpha = .5) +
  geom_point(aes(orca, incgraph, colour = avg.degree.f)) +
  facet_wrap(~model) + 
  scale_x_log10(breaks = 10^seq(-10,10)) +
  scale_y_log10(breaks = 10^seq(-10,10)) + 
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Number of nodes", y = "Orca Time", colour = "Placeholder")



g1 <- ggplot(summ.df) +
  geom_line(aes(num.nodes.l, speedup, colour = avg.degree.f), alpha = .5) +
  geom_point(aes(num.nodes.l, speedup, colour = avg.degree.f)) +
  facet_wrap(~model) + 
  scale_y_log10(breaks = c(round(min(summ.df$speedup)), 10^seq(0,10))) + 
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Number of nodes", y = "Orca Time / IncGraph Time", colour = "Placeholder")
g2 <- ggplot(summ.df) +
  geom_line(aes(avg.degree, speedup, colour = factor(num.nodes.l)), alpha = .5) +
  geom_point(aes(avg.degree, speedup, colour = factor(num.nodes.l))) +
  facet_wrap(~model) + 
  scale_x_continuous(trans = "log2", breaks = 2^seq(1, 6)) +
  scale_y_log10(breaks = c(round(min(summ.df$speedup)), 10^seq(0,10))) + 
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Degree", y = "Orca Time / IncGraph Time", colour = "Placeholder")

pdf(paste0(plots.folder, "/timingsplot.pdf"), 10, 7)
print(plot_grid(g1, g2, ncol = 1))
dev.off()


g1 <- ggplot(speedup.df) +
  geom_smooth(aes(num.nodes.l, speedup, colour = avg.degree.f, fill = avg.degree.f)) +
  facet_wrap(~model) + 
  scale_y_log10(breaks = c(round(min(summ.df$speedup)), 10^seq(0,10))) + 
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Number of nodes", y = "Orca Time / IncGraph Time", colour = "Placeholder", fill = "Placeholder")
g2 <- ggplot(speedup.df) +
  geom_smooth(aes(avg.degree, speedup, colour = factor(num.nodes.l), fill = factor(num.nodes.l))) +
  facet_wrap(~model) + 
  scale_x_continuous(trans = "log2", breaks = 2^seq(1, 6)) +
  scale_y_log10(breaks = c(round(min(summ.df$speedup)), 10^seq(0,10))) + 
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Degree", y = "Orca Time / IncGraph Time", colour = "Placeholder", fill = "Placeholder")

pdf(paste0(plots.folder, "/timingsplot_smooth.pdf"), 10, 7)
print(plot_grid(g1, g2, ncol = 1))
dev.off()

g <- ggplot(summ.df) + 
  geom_point(aes(h.index^3, incgraph, colour = model)) + 
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_colour_brewer(palette = "Set2") +
  coord_equal() +
  labs(x = "h-index³", y = "IncGraph time (s)")

pdf(paste0(plots.folder, "/hindex.pdf"), 8, 5)
print(g)
dev.off()

ggplot(summ.df %>% mutate(num.edges = avg.degree * num.nodes / 2, pct_network_change = speedup / num.edges)) +
  geom_line(aes(num.nodes.l, pct_network_change*100, colour = avg.degree.f), alpha = .5) +
  geom_point(aes(num.nodes.l, pct_network_change*100, colour = avg.degree.f)) +
  facet_wrap(~model) + 
  # scale_y_log10(breaks = c(round(min(summ.df$speedup)), 10^seq(0,10))) + 
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Number of nodes", y = "Max. percentage of network to change", colour = "Placeholder")

ggplot(summ.df) +
  geom_line(aes(num.nodes.l, incgraph, colour = avg.degree.f), alpha = .5) +
  geom_point(aes(num.nodes.l, incgraph, colour = avg.degree.f)) +
  facet_wrap(~model) + 
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Number of nodes", y = "IncGraph Time", colour = "Placeholder")


ggplot(summ.df) +
  geom_line(aes(num.nodes.l, orca, colour = avg.degree.f), alpha = .5) +
  geom_point(aes(num.nodes.l, orca, colour = avg.degree.f)) +
  facet_wrap(~model) + 
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Number of nodes", y = "IncGraph Time", colour = "Placeholder")
