install.packages("devtools")

install.packages("incgraph")

devtools::install_github("rcannood/GENIE3")

# install PRISM for easy communication with gridengine-based clusters
# or change all qsub_lapply's with lapply's if you do not wish to use PRISM
devtools::install_github("rcannood/PRISM") 

conf <- PRISM::create_qsub_config(
  remote = "yourremote",
  local_tmp_path = "/home/yourusername/.r2gridengine",
  remote_tmp_path = "/home/yourusername/.r2gridengine"
)
PRISM::set_default_qsub_config(conf)