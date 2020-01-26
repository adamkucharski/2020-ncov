# Process model for simulation --------------------------------------------

process_model <- function(t_start,t_end,dt,theta,simTab,simzetaA){
  
  # simTab <- storeL[,tt-1,]; t_start = 1; t_end = 2; dt = 1; simzetaA <- simzeta[,1]
  
  exposed_t1 <- simTab[,"exp1"] # input function
  exposed_t2 <- simTab[,"exp2"] # input function
  infectious_t1 <- simTab[,"inf1"] # input function
  infectious_t2 <- simTab[,"inf2"] # input function
  cases_t <- simTab[,"cases"] # input function
  reports_t <- simTab[,"reports"] # input function
  
  for(ii in seq((t_start+dt),t_end,dt) ){
    
    # transitions
    new_E1 <- (infectious_t1+infectious_t2)*theta[["beta"]]*simzetaA # stochastic transmission
    E1_to_E2 <- exposed_t1*theta[["incubation"]]*2 # as two compartments
    E2_to_I1 <- exposed_t2*theta[["incubation"]]*2
    I1_to_I2 <- theta[["recover"]]*infectious_t1*2
    I2_to_R <- theta[["recover"]]*infectious_t2*2
    I_to_C <- E2_to_I1
    C_to_Rep <- theta[["report"]]*cases_t
    
    # process model
    exposed_t1 <- exposed_t1 + new_E1 - E1_to_E2
    exposed_t2 <- exposed_t2 + E1_to_E2 - E2_to_I1
    infectious_t1 <- infectious_t1 + E2_to_I1 - I1_to_I2
    infectious_t2 <- infectious_t2 + I1_to_I2 - I2_to_R
    cases_t <- cases_t + I_to_C
    reports_t <- C_to_Rep
  }
  
  simTab[,"exp1"] <- exposed_t1 # output
  simTab[,"exp2"] <- exposed_t2 # output
  simTab[,"inf1"] <- infectious_t1 # output
  simTab[,"inf2"] <- infectious_t2 # output
  simTab[,"cases"] <- cases_t # output
  simTab[,"reports"] <- reports_t # output
  
  simTab
  
}

# SMC function --------------------------------------------

smc_model <- function(theta,nn){
  
  # nn = 100
  
  # Assumptions - using daily growth rate
  ttotal <- t_period
  dt <- 1
  t_length <- ttotal/dt
  
  storeL <- array(0,dim=c(nn,t_length, length(theta_initNames)),dimnames = list(NULL,NULL,theta_initNames))
  
  # Add initial condition
  storeL[,1,"inf1"] <- theta[["init_cases"]]
  storeL[,1,"inf2"] <- theta[["init_cases"]]
  #simzeta <- matrix(rlnorm(nn*t_length, mean = -theta[["betavol"]]^2/2, sd = theta[["betavol"]]),ncol=ttotal)
  simzeta <- matrix(rnorm(nn*t_length, mean = 0, sd = theta[["betavol"]]),nrow=ttotal)
  simzeta[1,] <- theta[["beta"]]*exp(simzeta[1,]) # ensure positive
  
  # PARTICLE FILTER GOES HERE
  # Latent variables
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
    
    # Add random walk on transmission ?
    simzeta[tt,] <- simzeta[tt-1,]*exp(simzeta[tt,])
    
    # run process model
    storeL[,tt,] <- process_model(tt-1,tt,dt,theta,storeL[,tt-1,],simzeta[tt,])
    
    # calculate weights
    caseDiff <- storeL[,tt,"reports"] - storeL[,tt-1,"reports"]
    w[,tt] <- AssignWeights(caseDiff,case_time[tt],nn) 
    
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
  Rep_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"reports"]
  I_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"inf1"]+storeL[l_sample[ttotal],ttotal,"inf2"]
  beta_traj[ttotal,] <- simzeta[ttotal,l_sample[ttotal]]
  
  for(ii in seq(ttotal,2,-1)){
    l_sample[ii-1] <- A[l_sample[ii],ii] # have updated indexing
    Rep_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"reports"]
    I_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"inf1"]+ storeL[l_sample[ii-1],ii-1,"inf2"]
    beta_traj[ii-1,] <- simzeta[ii-1,l_sample[ii-1]]
  }
 

  return(list(Rep_trace=Rep_traj,I_trace=I_traj,beta_trace=beta_traj,lik=likelihood0 ))
  
  
}


# Likelihood calc for SMC --------------------------------------------

AssignWeights <- function(x_val,case_data_tt,nn){
  
  # case_data_tt <- case_time[ttotal]; x_val <- caseDiff
  
  # x_val <- Rep_traj - c(0,head(Rep_traj,-1))
  
  # Do location by location
  # x_scaled <- x_val * passengers_daily / wuhan_area
  # x_expected <- sapply(x_scaled,function(x){x*as.numeric(top_risk$risk)}) %>% t() # expected exported cases in each location
  # # Note here rows are particles, cols are locations.
  # x_lam <- x_expected; dim(x_lam) <- NULL # flatten data on expectation
  # y_lam <- rep(case_data,nn); #dim(y_lam) <- NULL
  
  # Do for single international timeseries - probably more robust now
  x_scaled <- x_val * passengers_daily / wuhan_area
  
  x_expected <- sapply(x_scaled,function(x){x*as.numeric(top_risk$risk)}) %>% t() # expected exported cases in each location
  
  x_lam <- rowSums(x_expected)
  y_lam <- case_data_tt #sum(case_data)
  
  # Calculation likelihood
  loglik <- dpois(y_lam,lambda=x_lam,log=T)
  
  loglikSum <- rowSums(matrix(loglik,nrow=nn,byrow=T)) 
  
  exp(loglikSum) # convert to normal probability
  
}

# Likelihood calc (DEPRECATED) --------------------------------------------

likelihood_calc <- function(x_val,case_data_matrix){
  
  # Scale for air travel
  x_scaled <- x_val * passengers_daily / wuhan_area
  
  x_expected <- sapply(x_scaled,function(x){x*as.numeric(top_risk$risk)}) %>% t() # expected exported cases in each location
  
  #x_expected <- x_scaled[case_data$time]*case_data$export_probability # extract x values for observed data
  
  x_lam <- x_expected; dim(x_lam) <- NULL
  y_lam <- case_data_matrix; dim(y_lam) <- NULL
  
  # Define likelihood
  loglik <- dpois(y_lam,lambda=x_lam,log=T)
  
  sum(loglik)
  
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





