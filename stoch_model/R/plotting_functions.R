# Plot outputs from SMC --------------------------------------------

plot_outputs <- function(rep_plot,nn,cut_off){
  
  # rep_plot <- 10; nn <- 100; cut_off=5
  
  # Get median and 95%
  S_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  I_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  C_local_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  C_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  Rep_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  R0_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  
  for(kk in 1:rep_plot){
    output_smc <- smc_model(theta,nn)
    I_plot[,kk] <- output_smc$I_trace
    S_plot[,kk] <- output_smc$S_trace
    case_local_pos <- theta[["local_rep_prop"]]*(output_smc$C_local_trace - c(0,head(output_smc$C_local_trace,-1)))
    C_local_plot[,kk] <- case_local_pos #rpois(length(case_local_pos),lambda=case_local_pos)
    C_plot[,kk] <- output_smc$C_trace - c(0,head(output_smc$C_trace,-1))
    Rep_plot[,kk] <- output_smc$Rep_trace - c(0,head(output_smc$Rep_trace,-1)) # case difference
    R0_plot[,kk] <- (S_plot[,kk]/theta[["pop_travel"]])*output_smc$beta_trace/(theta[["recover"]])
  }
  
  # Calculate quantiles
  S_quantile <- apply(S_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))}) # thousands
  Inf_quantile <- apply(I_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})/1e3 # thousands

  # Local cases
  Case_local_quantile <- apply(C_local_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))}) 
  # International onset
  Case_quantile <- apply(C_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})
  # International confirmed
  Rep_quantile <- apply(Rep_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))}) 

  R0_quantile <- apply(R0_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})
  
  # Remove final few points (as estimation less reliable)
  S_quantileA <- S_quantile[,1:(ncol(R0_quantile)-cut_off)]/theta[["pop_travel"]]
  Inf_quantileA <- Inf_quantile[,1:(ncol(R0_quantile)-cut_off)]
  R0_quantileA <- R0_quantile[,1:(ncol(R0_quantile)-cut_off)]
  date_rangeA <- date_range[1:(length(date_range)-cut_off)]
  
  # - - - - - - - 
  # Calculate daily incidence
  #Case_diff_quantile <- Case_quantile[,1:ncol(Case_quantile)] - cbind(c(0,0,0),Case_quantile[,1:(ncol(Case_quantile)-1)])
  
  par(mfcol=c(3,2),mar=c(2,3,1,1),mgp=c(2,0.7,0))
  
  # Plot outputs
  a_col <- 0.4 # alpha
  xMin1 <- min(as.Date("2019-12-01"))
  xMin <- min(as.Date("2020-01-01"))
  xMax <- max(date_range)
  yMax <- max(Inf_quantileA[4,])
  
  # Plot outputs
  # Plot estimated local infections
  plot(date_rangeA,Inf_quantileA[1,],col="white",ylim=c(0,yMax),xlim=c(xMin1,xMax),xlab="",ylab="prevalence in Wuhan (thousands)")
  polygon(c(date_rangeA,rev(date_rangeA)),c(Inf_quantileA[2,],rev(Inf_quantileA[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_rangeA,rev(date_rangeA)),c(Inf_quantileA[1,],rev(Inf_quantileA[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_rangeA,Inf_quantileA[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")
  title(LETTERS[1],adj=0); letR = 2
  
  # Plot reproduction number
  plot(date_rangeA,R0_quantileA[1,],col="white",ylim=c(0,10),xlim=c(xMin1,xMax),xlab="",ylab="reproduction number")
  
  polygon(c(date_rangeA,rev(date_rangeA)),c(R0_quantileA[2,],rev(R0_quantileA[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_rangeA,rev(date_rangeA)),c(R0_quantileA[1,],rev(R0_quantileA[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_rangeA,R0_quantileA[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  lines(date_rangeA,1+0*R0_quantileA[3,],lty=2)
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,10),col="red")
  title(LETTERS[letR],adj=0); letR = letR + 1
  
  # Plot susceptibles
  # plot(date_rangeA,S_quantileA[1,],col="white",ylim=c(0,1.05),xlim=c(xMin,xMax),xlab="",ylab="propn susceptible")
  # polygon(c(date_rangeA,rev(date_rangeA)),c(S_quantileA[2,],rev(S_quantileA[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  # polygon(c(date_rangeA,rev(date_rangeA)),c(S_quantileA[1,],rev(S_quantileA[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  # lines(date_rangeA,S_quantileA[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  # 
  # lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")
  # title(LETTERS[2],adj=0)
  
  # Plot local cases onsets
  plot(date_range,Case_local_quantile[1,],col="white",ylim=c(0,30),xlim=c(xMin1,xMax),xlab="",ylab="local onsets")
  polygon(c(date_range,rev(date_range)),c(Case_local_quantile[2,],rev(Case_local_quantile[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(Case_local_quantile[1,],rev(Case_local_quantile[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,Case_local_quantile[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")

  points(case_data_wuhan$date,case_data_wuhan$number,pch=2)
  points(case_data_china$date,case_data_china$number)
  
  
  title(LETTERS[letR],adj=0); letR = letR +1
  box(lty = 1, col = 'black',lwd=2.5)
  
  # Plot international cases onsets
  plot(date_range,case_time,pch=19,ylim=c(0,8),xlim=c(xMin1,xMax),ylab="international onsets",col="white")
  
  polygon(c(date_range,rev(date_range)),c(fit_int_cases(Case_quantile[2,]),rev(fit_int_cases(Case_quantile[4,]))),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(fit_int_cases(Case_quantile[1,]),rev(fit_int_cases(Case_quantile[5,]))),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,fit_int_cases(Case_quantile[3,]),type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  points(date_range,case_data_onset_time)
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")
  title(LETTERS[letR],adj=0); letR = letR +1
  box(lty = 1, col = 'black',lwd=2.5)

  # Plot international cases confirmed
  plot(date_range,case_time,pch=19,ylim=c(0,8),xlim=c(xMin1,xMax),ylab="international confirmations",col="white")
  
  polygon(c(date_range,rev(date_range)),c(fit_int_cases(Rep_quantile[2,]),rev(fit_int_cases(Rep_quantile[4,]))),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(fit_int_cases(Rep_quantile[1,]),rev(fit_int_cases(Rep_quantile[5,]))),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,fit_int_cases(Rep_quantile[3,]),type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  points(date_range,case_time)
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")
  title(LETTERS[letR],adj=0); letR = letR +1

  dev.copy(png,paste("plots/case_inference.png",sep=""),units="cm",width=20,height=15,res=150)
  dev.off()
  
  
}

# Plot outputs from SMC --------------------------------------------

fit_int_cases <- function(x_val){
  
  x_scaled <- x_val #* passengers_daily / wuhan_area
  
  x_expected <- sapply(x_scaled,function(x){x*as.numeric(top_risk$risk)}) %>% t() # expected exported cases in each location
  x_lambda <- rowSums(x_expected)
  #rpois(length(x_lambda),lambda=x_lambda)
  x_lambda
  
}
