# Process model for simulation --------------------------------------------

process_model <- function(t_start,t_end,dt,theta,simTab,simzetaA,travelF){

  # simTab <- storeL[,tt-1,]; t_start = 1; t_end = 2; dt = 0.1; simzetaA <- simzeta[1,]; travelF=theta[["travel_frac"]]

  susceptible_t <- simTab[,"sus"] # input function
  exposed_t1 <- simTab[,"exp1"] # input function
  exposed_t2 <- simTab[,"exp2"] # input function
  tr_exposed_t1 <- simTab[,"tr_exp1"] # input function
  tr_exposed_t2 <- simTab[,"tr_exp2"] # input function

  infectious_t1 <- simTab[,"inf1"] # input function
  infectious_t2 <- simTab[,"inf2"] # input function
  
  tr_waiting_t <- simTab[,"tr_waiting"] # input function
  cases_t <- simTab[,"cases"] # input function
  reports_t <- simTab[,"reports"] # input function
  
  waiting_local_t <- simTab[,"waiting_local"] # input function
  cases_local_t <- simTab[,"cases_local"] # input function
  reports_local_t <- simTab[,"reports_local"] # input function

  # scale transitions
  inf_rate <- (simzetaA/theta[["pop_travel"]])*dt
  inc_rate <- theta[["incubation"]]*2*dt
  rec_rate <- theta[["recover"]]*2*dt
  rep_rate_local <- theta[["report_local"]]*dt
  rep_rate <- theta[["report"]]*dt
  travel_frac <- travelF

  for(ii in seq((t_start+dt),t_end,dt) ){

    # transitions
    S_to_E1 <- susceptible_t*(infectious_t1+infectious_t2)*inf_rate # stochastic transmission

    # Delay until symptoms
    E1_to_E2 <- exposed_t1*inc_rate # as two compartments
    E2_to_I1 <- exposed_t2*inc_rate
    
    E1_to_E2_tr <- tr_exposed_t1*inc_rate # as two compartments
    E2_to_I1_tr <- tr_exposed_t2*inc_rate

    # Delay until recovery
    I1_to_I2 <- infectious_t1*rec_rate
    I2_to_R <- infectious_t2*rec_rate
    
    # Delay until reported
    W_to_Rep <- tr_waiting_t*rep_rate
    W_to_Rep_local <- waiting_local_t*rep_rate_local

    # Process model for SEIR
    susceptible_t <- susceptible_t - S_to_E1
    exposed_t1 <- exposed_t1 + S_to_E1*(1-travel_frac) - E1_to_E2
    exposed_t2 <- exposed_t2 + E1_to_E2 - E2_to_I1
    tr_exposed_t1 <- tr_exposed_t1 + S_to_E1*travel_frac - E1_to_E2_tr
    tr_exposed_t2 <- tr_exposed_t2 + E1_to_E2_tr - E2_to_I1_tr

    infectious_t1 <- infectious_t1 + E2_to_I1 - I1_to_I2
    infectious_t2 <- infectious_t2 + I1_to_I2 - I2_to_R

    # Case tracking
    waiting_local_t <- waiting_local_t + E2_to_I1 - W_to_Rep_local
    cases_local_t <- cases_local_t + E2_to_I1
    reports_local_t <- reports_local_t + W_to_Rep_local

    tr_waiting_t <- tr_waiting_t + E2_to_I1_tr - W_to_Rep
    cases_t <- cases_t + E2_to_I1_tr
    reports_t <- reports_t + W_to_Rep
  }

  simTab[,"sus"] <- susceptible_t # output
  simTab[,"exp1"] <- exposed_t1 # output
  simTab[,"exp2"] <- exposed_t2 # output
  simTab[,"tr_exp1"] <- tr_exposed_t1 # output
  simTab[,"tr_exp2"] <- tr_exposed_t2 # output
  simTab[,"inf1"] <- infectious_t1 # output
  simTab[,"inf2"] <- infectious_t2 # output
  
  simTab[,"tr_waiting"] <- tr_waiting_t # output
  simTab[,"cases"] <- cases_t # output
  simTab[,"waiting_local"] <- waiting_local_t # output
  simTab[,"reports"] <- reports_t # output
  simTab[,"cases_local"] <- cases_local_t # output
  simTab[,"reports_local"] <- reports_local_t # output

  simTab

}

# SMC function --------------------------------------------

smc_model <- function(theta,nn,dt=1){

  # nn = 100;   dt <- 0.25

  # Assumptions - using daily growth rate
  ttotal <- t_period
  t_length <- ttotal

  storeL <- array(0,dim=c(nn,t_length, length(theta_initNames)),dimnames = list(NULL,NULL,theta_initNames))

  # Add initial condition
  #storeL[,1,"exp1"] <- theta[["init_cases"]]
  #storeL[,1,"exp2"] <- theta[["init_cases"]]
  storeL[,1,"inf1"] <- theta[["init_cases"]]/2
  storeL[,1,"inf2"] <- theta[["init_cases"]]/2
  storeL[,1,"sus"] <- theta[["pop_travel"]] - theta[["init_cases"]]

  #simzeta <- matrix(rlnorm(nn*t_length, mean = -theta[["betavol"]]^2/2, sd = theta[["betavol"]]),ncol=ttotal)
  simzeta <- matrix(rnorm(nn*t_length, mean = 0, sd = theta[["betavol"]]),nrow=ttotal)
  simzeta[1,] <- exp(simzeta[1,])*theta[["beta"]] # define IC

  # Latent variables
  S_traj = matrix(NA,ncol=1,nrow=ttotal)
  C_local_traj = matrix(NA,ncol=1,nrow=ttotal)
  Rep_local_traj = matrix(NA,ncol=1,nrow=ttotal)
  C_traj = matrix(NA,ncol=1,nrow=ttotal)
  Rep_traj = matrix(NA,ncol=1,nrow=ttotal)
  I_traj = matrix(NA,ncol=1,nrow=ttotal)
  beta_traj = matrix(NA,ncol=1,nrow=ttotal);
  w <- matrix(NA,nrow=nn,ncol=ttotal); w[,1] <- 1  # weights
  W <- matrix(NA,nrow=nn,ncol=ttotal)
  A <- matrix(NA,nrow=nn,ncol=ttotal) # particle parent matrix
  l_sample <- rep(NA,ttotal)
  lik_values <- rep(NA,ttotal)

  # Iterate through steps

  for(tt in 2:ttotal){
    
    # DEBUG  tt=2

    # Add random walk on transmission ?
    simzeta[tt,] <- simzeta[tt-1,]*exp(simzeta[tt,])

    # travel restrictions in place?
    if(tt<wuhan_travel_time){travelF <- theta[["travel_frac"]]}else{travelF <- 0}

    # run process model
    storeL[,tt,] <- process_model(tt-1,tt,dt,theta,storeL[,tt-1,],simzeta[tt,],travelF)

    # calculate weights
    w[,tt] <- AssignWeights(data_list,storeL,nn,theta,tt)

    #c_local_val,c_val,rep_val,local_case_data_tt,case_data_tt,rep_data_tt,nn){

    # normalise particle weights
    sum_weights <- sum(w[1:nn,tt])
    W[1:nn,tt] <- w[1:nn,tt]/sum_weights

    # resample particles by sampling parent particles according to weights:
    A[, tt] <- sample(1:nn,prob = W[1:nn,tt],replace = T)

    # DEPRECATED
    # for (j in 1:nn){
    #   locs <- pickA[cumsum_W >= rand_vals[j]]; A[j, tt] <- locs[1]
    # }

    # Resample particles for corresponding variables
    storeL[,tt,] <- storeL[ A[, tt] ,tt,]
    simzeta[tt,] <- simzeta[tt, A[, tt]] #- needed for random walk on beta


  } # END PARTICLE LOOP
  

  # Estimate likelihood:
  for(tt in 1:ttotal){
    lik_values[tt] = log(sum(w[1:nn,tt])) # log-likelihoods
  }

  likelihood0 = -ttotal*log(nn)+ sum(lik_values) # log-likelihoods

  # Sample latent variables:
  locs <-  sample(1:nn,prob = W[1:nn,tt],replace = T)
  l_sample[ttotal] <- locs[1]
  C_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"cases"]
  C_local_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"cases_local"]
  Rep_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"reports"]
  Rep_local_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"reports_local"]
  S_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"sus"]
  I_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"inf1"]+storeL[l_sample[ttotal],ttotal,"inf2"]
  beta_traj[ttotal,] <- simzeta[ttotal,l_sample[ttotal]]

  for(ii in seq(ttotal,2,-1)){
    l_sample[ii-1] <- A[l_sample[ii],ii] # have updated indexing
    C_local_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"cases_local"]
    Rep_local_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"reports_local"]
    C_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"cases"]
    Rep_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"reports"]
    S_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"sus"]
    I_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"inf1"]+ storeL[l_sample[ii-1],ii-1,"inf2"]
    beta_traj[ii-1,] <- simzeta[ii-1,l_sample[ii-1]]
  }
  
  # DEBUG  plot(Rep_traj[,1]-C_traj[,1])


  return(list(S_trace=S_traj,C_local_trace=C_local_traj,Rep_local_trace=Rep_local_traj,C_trace=C_traj,Rep_trace=Rep_traj,I_trace=I_traj,beta_trace=beta_traj,lik=likelihood0 ))


}


# Likelihood calc for SMC --------------------------------------------

AssignWeights <- function(data_list,storeL,nn,theta,tt){

  # c_local_val=case_localDiff; c_val=caseDiff; rep_val=repDiff;

  # Gather data
  local_case_data_tt <- data_list$local_case_data_onset[tt]
  local_case_data_conf_tt <- data_list$local_case_data_conf[tt]
  case_data_tt <- data_list$int_case_onset[tt]
  rep_data_tt <- data_list$int_case_conf[tt,]
  
  # Scale for reporting lag
  case_data_tt_scale <- 1#data_list$int_case_onset_scale[tt] # deprecated
  local_case_data_tt_scale <- 1#data_list$local_case_data_onset_scale[tt] # deprecated

  # Gather variables
  case_localDiff <- storeL[,tt,"cases_local"] - storeL[,tt-1,"cases_local"]
  rep_local <- storeL[,tt,"reports_local"]
  caseDiff <- storeL[,tt,"cases"] - storeL[,tt-1,"cases"]
  repDiff <- storeL[,tt,"reports"] - storeL[,tt-1,"reports"]
  c_local_val <- pmax(0,case_localDiff)
  c_val <- pmax(0,caseDiff)
  rep_val <- pmax(0,repDiff) # NOTE CHECK FOR POSITIVITY
  r_local_rep_cum <- rep_local

  # Local confirmed cases (by onset)

  if(!is.na(local_case_data_tt)){
    expected_val <- c_local_val*theta[["onset_prop"]]*theta[["local_rep_prop"]]*local_case_data_tt_scale # scale by reporting proportion and known onsets
    loglikSum_local_onset <- dpois(local_case_data_tt,lambda = expected_val,log=T)
  }else{
    loglikSum_local_onset <- 0
  }

  # Local confirmed cases (by confirmation) -- HOLD OUT FOR NOW AS LIKELIHOOD LOW

  # if(!is.na(local_case_data_conf_tt)){
  #   expected_val <- r_local_rep_cum*theta[["local_rep_prop"]] # scale by reporting proportion and known onsets
  #   loglikSum_local_conf <- dpois(local_case_data_conf_tt,lambda = r_local_rep_cum,log=T)
  # }else{
  #   loglikSum_local_conf <- 0
  # }

  # - - -
  # International confirmed cases (by country)

  # Do location by location
  x_scaled <- rep_val
  x_expected <- sapply(x_scaled,function(x){x*as.numeric(top_risk$risk)}); #x_expected <- t(x_expected) # expected exported cases in each location

  # # Note here rows are particles, cols are locations.
  x_lam <- x_expected; dim(x_lam) <- NULL # flatten data on expectation
  y_lam <- rep(rep_data_tt,nn); #dim(y_lam) <- NULL

  # Calculate likelihood
  loglik <- dpois(y_lam,lambda=x_lam,log=T)
  loglikSum_int_conf <- rowSums(matrix(loglik,nrow=nn,byrow=T))


  # - - -
  # International onsets (total)
  if(!is.na(case_data_tt)){
    x_scaled <- c_val
    x_expected <- sapply(x_scaled,function(x){x*as.numeric(top_risk$risk)}) %>% t() # expected exported cases in each location
    x_lam <- rowSums(x_expected)
    y_lam <- case_data_tt

    # Calculate likelihood
    loglik <- dpois(y_lam,lambda=x_lam,log=T)
    loglikSum_inf_onset <- rowSums(matrix(loglik,nrow=nn,byrow=T))
  }else{
    loglikSum_inf_onset <- 0
  }

  # - - -
  # Tally up likelihoods
  loglikSum <- loglikSum_local_onset   + loglikSum_inf_onset #+ loglikSum_int_conf #+ loglikSum_local_conf
  exp(loglikSum) # convert to normal probability

}


# Likelihood calc for SMC with single international timeseries (DEPRECATED) --------------------------------------------

AssignWeights_single_int <- function(x_val,case_data_tt,nn){

  # case_data_tt <- case_time[ttotal]; x_val <- caseDiff

  # Do for single international timeseries
  x_scaled <- x_val #* passengers_daily / wuhan_area

  x_expected <- sapply(x_scaled,function(x){x*as.numeric(top_risk$risk)}) %>% t() # expected exported cases in each location

  x_lam <- rowSums(x_expected)
  y_lam <- case_data_tt #sum(case_data)

  # Calculation likelihood
  loglik <- dpois(y_lam,lambda=x_lam,log=T)

  loglikSum <- rowSums(matrix(loglik,nrow=nn,byrow=T))

  exp(loglikSum) # convert to normal probability

}


# Simple simulation function calc (DEPRECATED) --------------------------------------------

simple_sim <- function(){


  # Assumptions - using daily growth rate
  ttotal <- t_period
  nn <- 100
  dt <- 1
  t_length <- ttotal/dt

  storeL <- array(0,dim=c(nn,t_length, length(theta_initNames)),dimnames = list(NULL,NULL,theta_initNames))

  # Initial condition
  storeL[,1,"inf"] <- theta[["init_cases"]]
  simzeta <- matrix(rlnorm(nn*t_length, mean = -theta[["betavol"]]^2/2, sd = theta[["betavol"]]),ncol=ttotal)

  for(ii in 2:ttotal){
    storeL[,ii,] <- process_model(ii-1,ii,dt,theta,storeL[,ii-1,],simzeta[,ii-1])
  }

  # Calculate likelihood
  log_lik <- apply(storeL[,,"cases"],1,function(x){likelihood_calc(x,case_data_matrix)})

  #Calculate relative probability of curves
  relative_prob <- exp(log_lik)/max(exp(log_lik))

  par(mfrow=c(3,1),mar=c(2,3,1,1),mgp=c(2,0.7,0))


  # Plot outputs

  plot(date_range,storeL[1,,1],col="white",ylim=c(0,1e5),xlab="",ylab="cases in China")
  for(ii in 1:nn){
    lines(date_range,storeL[ii,,"inf"],col=rgb(0,0,1,relative_prob[ii]))
  }

  # Plot international cases
  plot(date_range,case_time,pch=19,ylab="international cases")

  # Plot daily growth rate
  plot(date_range,simzeta[1,],col="white",ylim=c(0,3),xlab="",ylab="daily growth rate")
  for(ii in 1:nn){
    lines(date_range,simzeta[ii,],col=rgb(0,0,1,relative_prob[ii]))
  }

  dev.copy(png,paste("plots/case_inference.png",sep=""),units="cm",width=10,height=15,res=150)
  dev.off()


}

MLE_check <- function(p_name = "local_rep_prop", theta_tab,nn=1e3){

  # theta_tab <- seq(0.001,0.01,0.001)
  store_lik <- NULL

  for(ii in 1:length(theta_tab)){

    theta[[p_name]] <- theta_tab[ii]

    # Run SMC and output likelihooda
    output_smc <- smc_model(theta,
                            nn=1e3 # number of particles
    )
    store_lik <- rbind(store_lik,c(theta_tab[ii],output_smc$lik))

  }
  
  colnames(store_lik) <- c("param","lik")
  store_lik <- as_tibble(store_lik)

}

numerical_solver <- function(r0, k){

  fun <- function (s) (1 + (r0/k)*(1 - s))^(-k) - s
  solutions <- rootSolve::multiroot(fun, c(0, 1))$root

  realistic_sol <- min(solutions)
  return(realistic_sol)

}

