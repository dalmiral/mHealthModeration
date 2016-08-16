library("geepack")
library("zoo")

source("xgeepack.R")
source("xzoo.R")

system("R CMD SHLIB rsnmm.c")
dyn.load(ifelse(Sys.info()["sysname"] == "Windows", "rsnmm.dll", "rsnmm.so"))
source("rsnmm.R")
