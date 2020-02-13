# nCoV stochastic model
# Main script
# Author: AJ Kucharski (2020)

# Set up libraries and paths ----------------------------------------------

library(foreach)
library(doMC)
library(lubridate)
library(magrittr)
library(coda)
library(tidyverse)
library(rootSolve)
  
registerDoMC(4)  #change the 2 to your number of CPU cores

rm(list=ls(all=TRUE))

# - - -
# Set user-specific directory path and load datasets
if(Sys.info()["user"]=="adamkuchars" | Sys.info()["user"]=="adamkucharski") {
  setwd("~/Documents/GitHub/2020-nCov/stoch_model/")
  dropbox_path <- ""
}

# Load datasets, functions and parameters ----------------------------------------------

# - - -
# Load datasets
travel_data_mobs <- read_csv(paste0(dropbox_path,"data/connectivity_data_mobs.csv"))
#travel_data_mobs <- read_csv(paste0(dropbox_path,"data/connectivity_data_worldpop.csv"))
international_conf_data_in <- read_csv(paste0(dropbox_path,"data/international_case_data.csv"))
international_onset_data_in <- read_csv(paste0(dropbox_path,"data/time_series_WHO_report.csv"))
china_onset_data_in <- read_csv(paste0(dropbox_path,"data/time_series_data_bioRvix_Liu_et_al.csv"))
wuhan_onset_data_in <- read_csv(paste0(dropbox_path,"data/time_series_data_lancet_huang_et_al.csv"))
wuhan_onset_2020_01_30 <- read_csv(paste0(dropbox_path,"data/time_series_data_qui_li_nejm_wuhan.csv"))
wuhan_conf_data_in <- read_csv(paste0(dropbox_path,"data/time_series_HKU_Wuhan.csv"))
case_data_in <- international_conf_data_in
travel_data <- travel_data_mobs

# - - -
# Load model and plotting functions
source("R/model_functions.R")
source("R/plotting_functions.R")

# - - -
# Load model parameters
thetaR_IC <- read_csv("inputs/theta_initial_conditions.csv")
theta <- c( r0=as.numeric(thetaR_IC[thetaR_IC$param=="r0","value"]), # note this is only IC - SMC estimates this
            beta=NA,
            betavol=as.numeric(thetaR_IC[thetaR_IC$param=="betavol","value"]),
            gentime=as.numeric(thetaR_IC[thetaR_IC$param=="gentime","value"]), # not used currently
            incubation = 1/as.numeric(thetaR_IC[thetaR_IC$param=="incubation","value"]),
            report = 1/as.numeric(thetaR_IC[thetaR_IC$param=="report","value"]),
            report_local = 1/as.numeric(thetaR_IC[thetaR_IC$param=="report_local","value"]),
            recover = 1/as.numeric(thetaR_IC[thetaR_IC$param=="recover","value"]),
            init_cases=as.numeric(thetaR_IC[thetaR_IC$param=="init_cases","value"]),
            passengers=as.numeric(thetaR_IC[thetaR_IC$param=="outbound_travel","value"]),
            pop_travel=as.numeric(thetaR_IC[thetaR_IC$param=="population_travel","value"]),
            local_rep_prop=as.numeric(thetaR_IC[thetaR_IC$param=="local_rep_prop","value"]), # local propn reported
            onset_prop=as.numeric(thetaR_IC[thetaR_IC$param=="onset_prop","value"]), # propn onsets known
            travel_frac=NA
)

theta[["travel_frac"]] <- theta[["passengers"]]/theta[["pop_travel"]] # Estimate fraction that travel

theta[["beta"]] <- theta[["r0"]]*(theta[["recover"]]) # Scale initial value of R0

theta_initNames <- c("sus","tr_exp1","tr_exp2","exp1","exp2","inf1","inf2","tr_waiting","cases","reports","waiting_local","cases_local","reports_local") # also defines groups to use in model


# - - -
# Load timeseries -  specify travel data being used
# NOTE: USES REPORTING DELAY AS INPUT
source("R/load_timeseries_data.R",local=TRUE)



# Run set up check --------------------------------------------------------------

# - - -
# Run SMC and check likelihood
output_smc <- smc_model(theta,
                        nn=1e3, # number of particles
                        dt=0.25
)
output_smc$lik

# Run main outputs --------------------------------------------------------------

source("R/outputs_main.R")



