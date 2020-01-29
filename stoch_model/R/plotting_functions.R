# Helper functions --------------------------------------------

c.text<-function(x,sigF=3){
  bp1=signif(c(median(x),quantile(x,0.025),quantile(x,0.975)),sigF)
  paste(bp1[1]," (",bp1[2],"-",bp1[3],")",sep="")
}

c.nume<-function(x){
  bp1=c(median(x),quantile(x,0.025),quantile(x,0.975))
  as.numeric(bp1)
}

# Run SMC to get bootstrap estimates --------------------------------------------

run_fits <- function(rep_plot,nn,cut_off,dt,filename="1"){
  
  # rep_plot <- 10; nn <- 100; cut_off=0
  
  S_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  I_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  C_local_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  Rep_local_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  C_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  Rep_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  R0_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  
  for(kk in 1:rep_plot){
    output_smc <- smc_model(theta,nn,dt)
    I_plot[,kk] <- output_smc$I_trace
    S_plot[,kk] <- output_smc$S_trace
    case_local_pos <- theta[["local_rep_prop"]]*(output_smc$C_local_trace - c(0,head(output_smc$C_local_trace,-1)))
    C_local_plot[,kk] <- rpois(length(case_local_pos),lambda=case_local_pos)
    rep_local_pos <- theta[["local_rep_prop"]]*(output_smc$Rep_local_trace - c(0,head(output_smc$Rep_local_trace,-1)))
    Rep_local_plot[,kk] <- rpois(length(rep_local_pos),lambda=rep_local_pos)
    
    C_plot[,kk] <- fit_int_cases(output_smc$C_trace - c(0,head(output_smc$C_trace,-1)))
    Rep_plot[,kk] <- fit_int_cases(output_smc$Rep_trace - c(0,head(output_smc$Rep_trace,-1))) # case difference
    R0_plot[,kk] <- output_smc$beta_trace/(theta[["recover"]])
  }
  
  save(
    S_plot,
    I_plot,
    C_local_plot,
    Rep_local_plot,
    C_plot,
    Rep_plot,
    R0_plot,
  file=paste0("outputs/bootstrap_fit_",filename,".RData")) 
  
}

# Plot outputs from SMC --------------------------------------------

  
plot_outputs <- function(filename="1"){
  
  cut_off <- 0
  
  load(paste0("outputs/bootstrap_fit_",filename,".RData"))
  
  # Calculate quantiles
  S_quantile <- apply(S_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))}) # thousands
  Inf_quantile <- apply(I_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})/1e3 # thousands

  # Local cases
  Case_local_quantile <- apply(C_local_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))}) 
  Case_local_quantile_onset <- theta[["onset_prop"]]*Case_local_quantile
  
  Rep_local_quantile <- apply(Rep_local_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))}) 
  
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
  
  par(mfrow=c(3,2),mar=c(2,3,1,1),mgp=c(2,0.7,0))
  
  # Plot outputs
  a_col <- 0.4 # alpha
  xMin1 <- min(as.Date("2019-12-01"))
  xMin <- min(as.Date("2020-01-01"))
  xMax <- max(date_range)
  yMax <- max(Inf_quantileA[4,])
  
  # Plot outputs

  # Plot local cases onsets
  plot(date_range,Case_local_quantile_onset[1,],col="white",ylim=c(0,30),xlim=c(xMin1,xMax),xlab="",ylab="new onsets in Wuhan")
  polygon(c(date_range,rev(date_range)),c(Case_local_quantile_onset[2,],rev(Case_local_quantile_onset[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(Case_local_quantile_onset[1,],rev(Case_local_quantile_onset[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,Case_local_quantile_onset[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")
  text(labels="travel restrictions",x=wuhan_travel_restrictions-0.5,y=0.9*30,adj=1,col="red")
  
  text(labels="model estimate",x=xMin1,y=0.9*30,adj=0,col="blue")
  text(labels="fitted data",x=xMin1,y=0.8*30,adj=0,col="black")
  
  points(case_data_wuhan$date,case_data_wuhan$number,pch=17)
  points(case_data_china$date,case_data_china$number,pch=19)
  
  title(LETTERS[1],adj=0); letR = 2
  
  # - - -
  # Plot international cases onsets
  plot(date_range,case_time,pch=19,ylim=c(0,8),xlim=c(xMin1,xMax),ylab="new international onsets",col="white")
  
  polygon(c(date_range,rev(date_range)),c(Case_quantile[2,],rev(Case_quantile[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(Case_quantile[1,],rev(Case_quantile[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,Case_quantile[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  points(date_range,case_data_onset_time,pch=19)

  text(labels="model estimate",x=xMin1,y=0.9*8,adj=0,col="blue")
  text(labels="fitted data",x=xMin1,y=0.8*8,adj=0,col="black")
  text(labels="travel restrictions",x=wuhan_travel_restrictions-0.5,y=0.8*9,adj=1,col="red")
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")
  title(LETTERS[letR],adj=0); letR = letR +1
  
  # - - -
  # Plot local cases confirmed
  plot(date_range,case_time,pch=19,ylim=c(0,1000),xlim=c(xMin1,xMax),ylab="total cases confirmed in Wuhan",col="white")
  
  polygon(c(date_range,rev(date_range)),c(cumsum(Rep_local_quantile[2,]),rev(cumsum(Rep_local_quantile[4,]))),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(cumsum(Rep_local_quantile[1,]),rev(cumsum(Rep_local_quantile[5,]))),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,cumsum(Rep_local_quantile[3,]),type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  points(date_range,case_data_wuhan_conf_time)
  
  text(labels="model estimate",x=xMin1,y=0.9*1e3,adj=0,col="blue")
  text(labels="non-fitted data (used for validation)",x=xMin1,y=0.8*1e3,adj=0,col="black")
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")
  title(LETTERS[letR],adj=0); letR = letR +1
  
  # Plot international cases confirmed
  plot(date_range,case_time,pch=19,ylim=c(0,8),xlim=c(xMin1,xMax),ylab="new international cases confirmed",col="white")
  
  polygon(c(date_range,rev(date_range)),c(Rep_quantile[2,],rev(Rep_quantile[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(Rep_quantile[1,],rev(Rep_quantile[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,Rep_quantile[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  points(date_range,case_time,pch=1)
  
  text(labels="model estimate",x=xMin1,y=0.9*8,adj=0,col="blue")
  text(labels="non-fitted data (used for validation)",x=xMin1,y=0.8*8,adj=0,col="black")
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")
  title(LETTERS[letR],adj=0); letR = letR +1
  
  # Plot estimated local infections
  plot(date_rangeA,Inf_quantileA[1,],col="white",ylim=c(0,1.2*yMax),xlim=c(xMin1,xMax),xlab="",ylab="prevalence in Wuhan (thousands)")
  polygon(c(date_rangeA,rev(date_rangeA)),c(Inf_quantileA[2,],rev(Inf_quantileA[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_rangeA,rev(date_rangeA)),c(Inf_quantileA[1,],rev(Inf_quantileA[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_rangeA,Inf_quantileA[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")
  text(labels="model estimate",x=xMin1,y=1.1*yMax,adj=0,col="blue")
  
  title(LETTERS[letR],adj=0); letR = letR + 1
  
  # Plot reproduction number
  plot(date_rangeA,R0_quantileA[1,],col="white",ylim=c(0,10),xlim=c(xMin1,xMax),xlab="",ylab=expression(paste(R[0])))
  
  polygon(c(date_rangeA,rev(date_rangeA)),c(R0_quantileA[2,],rev(R0_quantileA[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_rangeA,rev(date_rangeA)),c(R0_quantileA[1,],rev(R0_quantileA[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_rangeA,R0_quantileA[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  lines(date_rangeA,1+0*R0_quantileA[3,],lty=2)
  
  text(labels="model estimate",x=xMin1,y=9,adj=0,col="blue")
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,10),col="red")
  title(LETTERS[letR],adj=0); letR = letR + 1
  
  # output figure
  dev.copy(png,paste("plots/cases_inference.png",sep=""),units="cm",width=20,height=15,res=150)
  #dev.copy(pdf,paste("plots/cases_inference.pdf",sep=""),width=8,height=8)
  dev.off()
  
  
}

# Plot dispersion --------------------------------------------

plot_dispersion <- function(filename="1"){
  
  # filename="1"
  
  load(paste0("outputs/bootstrap_fit_",filename,".RData"))
  
  # Extract R credible interval
  period_interest <- as.Date(c("2020-01-10","2020-01-14"))
  
  par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(2,0.7,0))
  
  for(ii in 1:2){
    
    index_pick <- match(period_interest,date_range)
    R0_all <- R0_plot[index_pick[ii],]
    dim(R0_all) <- NULL # collapse data
    
    R0_CrI <- quantile(R0_all,c(0.025,0.25,0.5,0.75,0.975))
    
    MERS_k <- 0.26
    SARS_k <- 0.16
    k_seq <- seq(0.01,0.5,0.01)
  
    # Outbreak calcs
    R0_med <- 1-sapply(k_seq,function(x){numerical_solver(R0_CrI[3],x)})
    R0_CrI_1 <- 1-sapply(k_seq,function(x){numerical_solver(R0_CrI[1],x)})
    R0_CrI_2 <- 1-sapply(k_seq,function(x){numerical_solver(R0_CrI[5],x)})
    R0_CrI_1_50 <- 1-sapply(k_seq,function(x){numerical_solver(R0_CrI[2],x)})
    R0_CrI_2_50 <- 1-sapply(k_seq,function(x){numerical_solver(R0_CrI[4],x)})
    
    # Plot results
    plot(k_seq,R0_med,type="l",ylim=c(0,1),xlab=c("extent of homogeneity in transmission (k)"),ylab="probability of large outbreak",col="white",xaxs="i",yaxs="i")
    polygon(c(k_seq,rev(k_seq)),c(R0_CrI_1,rev(R0_CrI_2)),lty=0,col=rgb(0,0.3,1,0.35))
    polygon(c(k_seq,rev(k_seq)),c(R0_CrI_1_50,rev(R0_CrI_2_50)),lty=0,col=rgb(0,0.3,1,0.35))
    lines(k_seq,R0_med,col="blue")
    
    lines(c(MERS_k,MERS_k),c(-1,10),lty=2); text(labels="MERS-CoV",x=1.02*MERS_k,y=0.7,adj=0,col="black")
    lines(c(SARS_k,SARS_k),c(-1,10),lty=2); text(labels="SARS",x=1.02*SARS_k,y=0.6,adj=0,col="black")
    title(LETTERS[ii],adj=0)

  }
  
  
  dev.copy(png,paste("plots/calc_1.png",sep=""),units="cm",width=20,height=10,res=150)
  dev.off()
  
  
}

# Output R0 estimates over time --------------------------------------------

r0_value_output <- function(filename="1"){
  
  # filename="1"
  
  load(paste0("outputs/bootstrap_fit_",filename,".RData"))

  
  period_interest <- as.Date("2020-01-01","")
  
  match(period_interest,date_range)
  
  R0_plot #TO DO
  

  
}

# Exporter case scaling function --------------------------------------------

fit_int_cases <- function(x_val){
  
  x_scaled <- x_val #* passengers_daily / wuhan_area
  
  x_expected <- sapply(x_scaled,function(x){x*as.numeric(top_risk$risk)}) %>% t() # expected exported cases in each location
  x_lambda <- rowSums(x_expected)
  #x_lambda
  rpois(length(x_lambda),lambda=x_lambda)
  


}
