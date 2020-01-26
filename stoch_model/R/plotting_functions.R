# Plot outputs from SMC --------------------------------------------

plot_outputs <- function(rep_plot,nn){
  
  # Get median and 95%
  
  rep_plot <- 100 # how many SMC samples
  
  I_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  C_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  R0_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  
  for(kk in 1:rep_plot){
    output_smc <- smc_model(theta,nn)
    I_plot[,kk] <- output_smc$I_trace
    C_plot[,kk] <- output_smc$C_trace# - c(0,head(output_smc$C_trace,-1)) # case difference
    R0_plot[,kk] <- output_smc$beta_trace*theta[["beta"]]/(theta[["recover"]]+theta[["incubation"]])
  }
  
  # Calculate quantiles
  Inf_quantile <- apply(I_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})
  Case_quantile <- apply(C_plot,1,function(x){quantile(x,c(0.025,0.5,0.975))})
  R0_quantile <- apply(R0_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})
  
  par(mfrow=c(3,1),mar=c(2,3,1,1),mgp=c(2,0.7,0))
  
  # Plot outputs
  a_col <- 0.4 # alpha
  
  # Plot estimated csaes
  plot(date_range,Inf_quantile[1,],col="white",ylim=c(0,1e6),xlab="",ylab="prevalence in Wuhan")
  polygon(c(date_range,rev(date_range)),c(Inf_quantile[2,],rev(Inf_quantile[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(Inf_quantile[1,],rev(Inf_quantile[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,Inf_quantile[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  
  # Plot international cases
  plot(date_range,case_time,pch=19,ylab="international cases")
  
  # Plot daily growth rate
  plot(date_range,R0_quantile[1,],col="white",ylim=c(0,4),xlab="",ylab="reproduction number")
  
  polygon(c(date_range,rev(date_range)),c(R0_quantile[2,],rev(R0_quantile[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(R0_quantile[1,],rev(R0_quantile[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,R0_quantile[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  lines(date_range,1+0*R0_quantile[3,],lty=2)
  
  dev.copy(png,paste("plots/case_inference.png",sep=""),units="cm",width=10,height=15,res=150)
  dev.off()
  
  
}
