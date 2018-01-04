library(GENIE3)
library(tidyverse)
library(incgraph)
library(jsonlite)

out.folder <- "data/genie3"
scratch.folder <- "data_tmp/genie3"
plots.folder <- "plots/genie3"
dir.create(out.folder, recursive = T)
dir.create(scratch.folder, recursive = T)
dir.create(plots.folder, recursive = T)

#' READ ECOLI DATA
expression <- read.table(paste0(out.folder, "/expression/ecoli/E.tsv"), header = T, sep = "\t")
gmap <- read.table(paste0(out.folder, "/expression/ecoli/Gmap.tsv"), header = F, sep = "\t", quote = "")
colnames(gmap) <- c("symbol", "anon")
tfs <- fromJSON(paste0(out.folder, "/expression/ecoli/regulators.json"))
adjacency <- read.table(paste0(out.folder, "/expression/ecoli/adjacency.csv"), header = T, sep = "\t", row.names = 1) %>% t

adjacency <- adjacency[rownames(adjacency) %in% gmap$symbol, colnames(adjacency) %in% gmap$symbol]
rownames(adjacency) <- gmap$anon[match(rownames(adjacency), gmap$symbol)]
colnames(adjacency) <- gmap$anon[match(colnames(adjacency), gmap$symbol)]

gene.sel <- intersect(colnames(expression), union(rownames(adjacency), colnames(adjacency)))
expression <- expression[,gene.sel]
tfs <- tfs[tfs %in% gene.sel]
adjacency2 <- matrix(0, nrow = length(tfs), ncol = length(gene.sel), dimnames = list(tfs, gene.sel))
adjacency <- adjacency[rownames(adjacency) %in% tfs, colnames(adjacency) %in% gene.sel]
adjacency2[rownames(adjacency), colnames(adjacency)] <- adjacency
adjacency <- adjacency2
rm(adjacency2)

ecoli <- list(
  name = "E. coli",
  expression = expression,
  adjacency = adjacency,
  tfs = tfs
)

rm(expression, gmap, tfs, adjacency, gene.sel)

#' READ YEAST DATA
expression <- read.table(paste0(out.folder, "/expression/yeast/E.tsv"), header = T, sep = "\t", row.names = 1, check.names = F)
gmap <- read.table(paste0(out.folder, "/expression/yeast/Gmap.tsv"), header = F, sep = "\t", quote = "", stringsAsFactors = F)
colnames(gmap) <- c("symbol", "anon")
tfs <- fromJSON(paste0(out.folder, "/expression/yeast/regulators.json"))
adjacency <- read.table(paste0(out.folder, "/expression/yeast/adjacency.csv"), header = T, sep = "\t", row.names = 1, check.names = F) %>% t

adjacency <- adjacency[rownames(adjacency) %in% gmap$symbol, colnames(adjacency) %in% gmap$symbol]
rownames(adjacency) <- gmap$anon[match(rownames(adjacency), gmap$symbol)]
colnames(adjacency) <- gmap$anon[match(colnames(adjacency), gmap$symbol)]

gene.sel <- intersect(colnames(expression), union(rownames(adjacency), colnames(adjacency)))
expression <- expression[,gene.sel]
tfs <- tfs[tfs %in% gene.sel]
adjacency2 <- matrix(0, nrow = length(tfs), ncol = length(gene.sel), dimnames = list(tfs, gene.sel))
adjacency <- adjacency[rownames(adjacency) %in% tfs, colnames(adjacency) %in% gene.sel]
adjacency2[rownames(adjacency), colnames(adjacency)] <- adjacency
adjacency <- adjacency2
rm(adjacency2)

yeast <- list(
  name = "S. cerevisiae",
  expression = expression,
  adjacency = adjacency,
  tfs = tfs
)

rm(expression, gmap, tfs, adjacency, gene.sel)

#' COMBINE DATA
datasets <- list(
  yeast = yeast,
  ecoli = ecoli
)
rm(yeast, ecoli)

#' RUN GENIE3
for (n in names(datasets)) {
  cat("Processing ", n, "\n", sep="")
  d <- datasets[[n]]
  datasets[[n]]$genie3 <- GENIE3::genie3(d$expression, regulators = d$tfs, mc.cores = "qsub")
}

save(datasets, file = paste0(out.folder, "/datasets.RData"))
