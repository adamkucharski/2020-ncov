# nCoV stochastic model
# Main script
# Author: AJ Kucharski (2020)

# Set up libraries and paths ----------------------------------------------

library(foreach)
library(doMC)
library(mvtnorm)
library(lubridate)
library(magrittr)
library(coda)
library(tidyverse)
  
registerDoMC(4)  #change the 2 to your number of CPU cores

rm(list=ls(all=TRUE))

# Set user-specific directory path and load datasets
if(Sys.info()["user"]=="adamkuchars" | Sys.info()["user"]=="adamkucharski") {
  setwd("~/Documents/GitHub/2020-nCov/stoch_model/")
  dropbox_path <- "~/Dropbox/Shared_nCov2019/"
}

# Load datasets, functions and parameters ----------------------------------------------

# Load datasets
travel_data <- read_csv(paste0(dropbox_path,"data_sources/mobility_data/mobs_connectivity_data.csv"))
case_data <- read_csv(paste0(dropbox_path,"data_sources/case_data/international_case_data.csv"))

# Load model and plotting functions
source("R/model_functions.R")
source("R/plotting_functions.R")

# Define time period and global transport paramters
passengers_daily <- 3300
wuhan_area <- 19e6

# Load timeseries
source("R/load_timeseries_data.R",local=TRUE)

# Load model parameters
thetaR_IC <- read_csv("inputs/theta_initial_conditions.csv")
theta <- c( r0=as.numeric(thetaR_IC[thetaR_IC$param=="r0","value"]),
            beta=NA,
            betavol=as.numeric(thetaR_IC[thetaR_IC$param=="betavol","value"]),
            gentime=as.numeric(thetaR_IC[thetaR_IC$param=="gentime","value"]), # not used currently
            report=1/as.numeric(thetaR_IC[thetaR_IC$param=="report","value"]),
            recover=1/as.numeric(thetaR_IC[thetaR_IC$param=="recover","value"]),
            init_cases=as.numeric(thetaR_IC[thetaR_IC$param=="init_cases","value"]))

#theta[["r0"]] <- 2.5
#theta[["betavol"]] <- 0.3

theta[["beta"]] <- theta[["r0"]]*theta[["recover"]]

theta_initNames <- c("inf","cases","reports") # also defines groups to use in model


# Run models --------------------------------------------------------------

# Run SMC and output likelihood
output_smc <- smc_model(theta,
                        nn=1e2 # number of particles
                        )
output_smc$lik

# Run multiple SMC and output plots
plot_outputs(rep_plot=10, # number of repeats
             nn=1e2 # number of particles
             )




