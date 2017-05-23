#!/usr/bin/env Rscript

# 20170202_figuresofmerit_batch.R
# wrapper that feeds peptide-subsetted dataframe to nonlinear_quantlim()
# last updated: 2017-04-12

# NOTE!! BEFORE YOU RUN THIS R SCRIPT, INSTALL THIS PACKAGE:
#install.packages("/net/maccoss/vol6/home/lpino/proj/MSstats_3.7.3.tar.gz", lib = "/net/maccoss/vol6/home/lpino/R/3.3.0/libs", type="source")
# test MSstats install
#?nonlinear_quantlim()

rm(list = ls())
library(readr)
library(tidyr)
library(dplyr)
library(gplots)
library(lme4)
library(ggplot2)
library(ggrepel)
library(reshape)
library(reshape2)
library(data.table)
library(Rcpp)
library(survival)
library(limma)
library(marray)
library(preprocessCore)
library(MSnbase)
library(MSstats, lib.loc= "/net/maccoss/vol6/home/lpino/R/3.3.0/libs")
library(minpack.lm)

# test MSstats version
packageVersion("MSstats")

args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

curve.df <- read.csv(args[1], header=TRUE, stringsAsFactors = FALSE)
fo_dir <- args[2]
precursor.start <- args[3]
precursor.end <- args[4]
fo <- paste(fo_dir, "/", precursor.start, ".txt", sep="")

# 
# CALCULATION OF LOD/LOQ WITH A NONLINEAR FIT
#

# vector of peptides to iterate through
peptides <- unique(curve.df$NAME)
peptide.batch <- peptides[precursor.start:precursor.end]

# calculate LOD/LOQ for each peptide, storing plots in a designated directory
counter <- 1 # inititalize counter for "for" loop below
for (peptide in peptide.batch){

  time <- Sys.time()
  print(paste("counter=",counter)) # sanity check
  print(time)
  print(peptide)

  counter <- counter + 1 # increment counter to ensure update
  #df_in contains data for peptide i  
  df_in <- curve.df %>% filter(NAME == peptide)
  
  #Count the number of blank samples for the peptide (with a non NA intensity) [commented out 'and with different values']
  df_blank = df_in %>% filter(CONCENTRATION == 0 & !is.na(INTENSITY)) #%>%  distinct(INTENSITY) 
  
  #n_blank = number of "acceptable" blank samples:
  n_blank = nrow(df_blank)
  print(paste("n_blank=",n_blank))
  if(n_blank <= 1) {next}
  
  df_out <- nonlinear_quantlim(df_in)
  try(plot_quantlim(spikeindata = df_in, quantlim_out = df_out, dir_output=fo_dir))

  # write the nonlinear_quantlim() results to an outfile for more downstream processing
  if(counter == 1){
    write.table(df_out, file=fo, append=FALSE, col.names=TRUE)
  } else {
    try(write.table(df_out, file=fo, append=TRUE, col.names=FALSE))
  }

  #print(df_out)

}

close(fo)