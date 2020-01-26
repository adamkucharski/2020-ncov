# Plot outputs from SMC --------------------------------------------

plot_outputs <- function(rep_plot,nn,cut_off){
  
  # rep_plot <- 10; nn <- 100; cut_off=5
  
  # Get median and 95%

  I_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  Rep_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  R0_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  
  for(kk in 1:rep_plot){
    output_smc <- smc_model(theta,nn)
    I_plot[,kk] <- output_smc$I_trace
    Rep_plot[,kk] <- output_smc$Rep_trace - c(0,head(output_smc$Rep_trace,-1)) # case difference
    R0_plot[,kk] <- output_smc$beta_trace/(theta[["recover"]]+theta[["incubation"]])
  }
  
  # Calculate quantiles
  Inf_quantile <- apply(I_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})/1e3 # thousands
  Rep_quantile <- apply(Rep_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))}) # thousands
  Rep_quantile_plot <- apply(Rep_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})/1e3 # thousands
  R0_quantile <- apply(R0_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})
  
  # Remove final few points (as estimation less reliable)
  R0_quantileA <- R0_quantile[,1:(ncol(R0_quantile)-cut_off)]
  date_rangeA <- date_range[1:(length(date_range)-cut_off)]
  
  # Calculate daily incidence
  #Case_diff_quantile <- Case_quantile[,1:ncol(Case_quantile)] - cbind(c(0,0,0),Case_quantile[,1:(ncol(Case_quantile)-1)])
  
  par(mfrow=c(3,1),mar=c(2,3,1,1),mgp=c(2,0.7,0))
  
  # Plot outputs
  a_col <- 0.4 # alpha
  xMin <- min(date_range)
  xMax <- max(date_range)

  # Plot estimated infections
  plot(date_range,Inf_quantile[1,],col="white",ylim=c(0,1e3),xlim=c(xMin,xMax),xlab="",ylab="prevalence in Wuhan (thousands)")
  polygon(c(date_range,rev(date_range)),c(Inf_quantile[2,],rev(Inf_quantile[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(Inf_quantile[1,],rev(Inf_quantile[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,Inf_quantile[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  
  
  # Plot international cases
  plot(date_range,case_time,pch=19,xlim=c(xMin,xMax),ylab="international cases",col="white")
  
  polygon(c(date_range,rev(date_range)),c(fit_int_cases(Rep_quantile[2,]),rev(fit_int_cases(Rep_quantile[4,]))),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(fit_int_cases(Rep_quantile[1,]),rev(fit_int_cases(Rep_quantile[5,]))),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,fit_int_cases(Rep_quantile[3,]),type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")

  points(date_range,case_time)
  
  # Plot reproduction number
  plot(date_rangeA,R0_quantileA[1,],col="white",ylim=c(0,5),xlim=c(xMin,xMax),xlab="",ylab="reproduction number")
  
  polygon(c(date_rangeA,rev(date_rangeA)),c(R0_quantileA[2,],rev(R0_quantileA[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_rangeA,rev(date_rangeA)),c(R0_quantileA[1,],rev(R0_quantileA[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_rangeA,R0_quantileA[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  lines(date_rangeA,1+0*R0_quantileA[3,],lty=2)
  
  dev.copy(png,paste("plots/case_inference.png",sep=""),units="cm",width=10,height=15,res=150)
  dev.off()
  
  
}

# Plot outputs from SMC --------------------------------------------

fit_int_cases <- function(x_val){
  
  x_scaled <- x_val * passengers_daily / wuhan_area
  
  x_expected <- sapply(x_scaled,function(x){x*as.numeric(top_risk$risk)}) %>% t() # expected exported cases in each location
  
  rowSums(x_expected)
  
}
