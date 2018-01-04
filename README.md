incgraph-scripts
================

This repository contains the code and the data required to reproduce the results shown in the IncGraph manuscript (under review).

Data
----

Finally, you should uncompress the zip files at [data/genie3/expression/\*.7z](https://github.com/rcannood/incgraph-scripts/tree/master/data/genie3/expression).

Install
-------

A few packages need to be installed.

``` r
install.packages("devtools")
install.packages("incgraph")
devtools::install_github("rcannood/GENIE3")
```

Cluster
-------

The PRISM package allows for easy communication between Rstudio and a gridengine-based cluster. For this, you will need to set up an openssh connection with your own gridengine cluster. If you do not wish to make use of PRISM's qsub\_lapply, you can change all qsub\_lapply's in the code with regular lapply's or parallelMap's.

``` r
devtools::install_github("rcannood/PRISM") 

conf <- PRISM::create_qsub_config(
  remote = "yourremote",
  local_tmp_path = "/home/yourusername/.r2gridengine",
  remote_tmp_path = "/home/yourusername/.r2gridengine"
)
PRISM::set_default_qsub_config(conf)
```
