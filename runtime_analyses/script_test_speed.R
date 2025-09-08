library(tidyverse)
library(microseq)
library(Rsearch)
library(dada2)


sample.tbl <- read_delim("metadata.tsv")
n.sim <- 5
time.tbl <- data.frame(dada2 = rep(0, n.sim),
                       rsearch_unoise = rep(0, n.sim),
                       rsearch_id99 = rep(0, n.sim),
                       rsearch_id95 = rep(0, n.sim))


#################
### dada2
###
for(r in 1:n.sim){
  if(dir.exists("tmp")){
    ff <- list.files("tmp", full.names = T)
    file.remove(ff)
  } else {
    dir.create("tmp")
  }
  
  t0 <- Sys.time()
  source("dada2.R")
  t1 <- Sys.time()
  time.tbl$dada2[r] <- difftime(t1, t0, units = "secs")
}


###################
### Rsearch UNOISE
###
for(r in 1:n.sim){
  if(dir.exists("tmp")){
    ff <- list.files("tmp", full.names = T)
    OK <- file.remove(ff)
  } else {
    dir.create("tmp")
  }
  
  t0 <- Sys.time()
  source("rsearch_unoise.R")
  t1 <- Sys.time()
  time.tbl$rsearch_unoise[r] <- difftime(t1, t0, units = "secs")
}


#####################################
### Rsearch identity 99% clustering
###
identity <- 0.99
for(r in 1:n.sim){
  if(dir.exists("tmp")){
    ff <- list.files("tmp", full.names = T)
    file.remove(ff)
  } else {
    dir.create("tmp")
  }
  
  t0 <- Sys.time()
  source("rsearch_identity.R")
  t1 <- Sys.time()
  print(difftime(t1, t0, units = "secs"))
  time.tbl$rsearch_id99[r] <- difftime(t1, t0, units = "secs")
}


#####################################
### Rsearch identity 95% clustering
###
identity <- 0.95
for(r in 1:n.sim){
  if(dir.exists("tmp")){
    ff <- list.files("tmp", full.names = T)
    file.remove(ff)
  } else {
    dir.create("tmp")
  }
  
  t0 <- Sys.time()
  source("rsearch_identity.R")
  t1 <- Sys.time()
  time.tbl$rsearch_id95[r] <- difftime(t1, t0, units = "secs")
}

save(time.tbl, file = "time_mac_Cassandra_v2.tbl.RData")

#####################################
### Rsearch plot_base_quality
###
for(r in 1:n.sim){
  t0 <- Sys.time()
  source("plot_base_quality.R")
  t1 <- Sys.time()
  print(difftime(t1, t0, units = "secs"))
}

#####################################
### Rsearch vs_merging_lengths
###
for(r in 1:n.sim){
  t0 <- Sys.time()
  source("vs_merging_lengths.R")
  t1 <- Sys.time()
  print(difftime(t1, t0, units = "secs"))
}

#####################################
### Rsearch vs_optimize_truncqual
###
for(r in 1:n.sim){
  t0 <- Sys.time()
  source("vs_optimize_truncqual.R")
  t1 <- Sys.time()
  print(difftime(t1, t0, units = "secs"))
}

#####################################
### Rsearch vs_optimize_truncee_rate
###
for(r in 1:n.sim){
  t0 <- Sys.time()
  source("vs_optimize_truncee_rate.R")
  t1 <- Sys.time()
  print(difftime(t1, t0, units = "secs"))
}
