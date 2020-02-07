# Helper functions --------------------------------------------

c.text<-function(x,sigF=3){
  bp1=signif(c(median(x),quantile(x,0.025),quantile(x,0.975)),sigF)
  paste(bp1[1]," (",bp1[2],"-",bp1[3],")",sep="")
}

c.nume<-function(x){
  bp1=c(median(x),quantile(x,0.025),quantile(x,0.975))
  as.numeric(bp1)
}

bin_conf <- function(x,n){
  htest <- binom.test(x,n)
  h_out <- c(x/n, htest$conf.int[1], htest$conf.int[2])
  h_out
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
  
  # filename="1"
  
  cut_off <- 0 #end_date - as.Date("2020-01-23")
  
  load(paste0("outputs/bootstrap_fit_",filename,".RData"))
  
  # Remove NA fits
  S_plot = S_plot[,!is.na(S_plot[t_period,])]
  I_plot = I_plot[,!is.na(I_plot[t_period,])]
  C_local_plot = C_local_plot[,!is.na(C_local_plot[t_period,])]
  Rep_local_plot = Rep_local_plot[,!is.na(Rep_local_plot[t_period,])]
  C_plot = C_plot[,!is.na(C_plot[t_period,])]
  Rep_plot = Rep_plot[,!is.na(Rep_plot[t_period,])]
  R0_plot = R0_plot[,!is.na(R0_plot[t_period,])]
  
  # Calculate quantiles
  S_quantile <- apply(S_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})/theta[["pop_travel"]] # proportion
  Inf_quantile <- apply(I_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})/theta[["pop_travel"]] # proportion

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
  S_quantileA <- S_quantile[,1:(ncol(R0_quantile)-cut_off)]
  Inf_quantileA <- Inf_quantile[,1:(ncol(R0_quantile)-cut_off)]
  R0_quantileA <- R0_quantile[,1:(ncol(R0_quantile)-cut_off)]
  date_rangeA <- date_range[1:(length(date_range)-cut_off)]
  Case_local_quantile_onsetA <- Case_local_quantile_onset[,1:(length(date_range)-cut_off)]
  
  # - - - - - - - 
  # Calculate daily incidence
  #Case_diff_quantile <- Case_quantile[,1:ncol(Case_quantile)] - cbind(c(0,0,0),Case_quantile[,1:(ncol(Case_quantile)-1)])
  
  par(mfrow=c(4,2),mar=c(2,3,1,1),mgp=c(2,0.7,0))
  
  # Plot outputs
  a_col <- 0.4 # alpha
  xMin1 <- min(as.Date("2019-12-01"))
  xMin <- min(as.Date("2020-01-01"))
  xMax <- end_date #max(date_range)
  
  # Plot local case onsets
  ym1 <- 40
  plot(date_rangeA,Case_local_quantile_onsetA[1,],col="white",ylim=c(0,ym1),xlim=c(xMin1,xMax),xlab="",ylab="new onsets in Wuhan")
  polygon(c(date_rangeA,rev(date_rangeA)),c(Case_local_quantile_onsetA[2,],rev(Case_local_quantile_onsetA[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_rangeA,rev(date_rangeA)),c(Case_local_quantile_onsetA[1,],rev(Case_local_quantile_onsetA[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_rangeA,Case_local_quantile_onsetA[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")
  text(labels="travel restrictions",x=wuhan_travel_restrictions+0.5,y=0.9*ym1,adj=0,col="red")
  
  text(labels="model estimate",x=xMin1,y=0.9*ym1,adj=0,col="blue")
  text(labels="fitted data",x=xMin1,y=0.8*ym1,adj=0,col="black")
  
  points(case_data_wuhan$date,case_data_wuhan$number,pch=17)
  points(case_data_china$date,case_data_china$number,pch=18)
  #points(date_range,case_data_wuhan_2_time,pch=19)
  
  title(LETTERS[1],adj=0); letR = 2

  
  # - - -
  # Plot international cases onsets
  ym1 <- 10
  plot(date_range,case_time,pch=19,ylim=c(0,ym1),xlim=c(xMin1,xMax),ylab="new international onsets",col="white")
  
  polygon(c(date_range,rev(date_range)),c(Case_quantile[2,],rev(Case_quantile[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(Case_quantile[1,],rev(Case_quantile[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,Case_quantile[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  points(date_range,case_data_onset_time,pch=19)

  text(labels="model estimate",x=xMin1,y=0.9*ym1,adj=0,col="blue")
  text(labels="fitted data",x=xMin1,y=0.8*ym1,adj=0,col="black")
  text(labels="travel restrictions",x=wuhan_travel_restrictions+0.5,y=0.8*ym1,adj=0,col="red")
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")
  title(LETTERS[letR],adj=0); letR = letR +1
  
  # - - -
  # Plot local cases confirmed
  # plot(date_range,case_time,pch=19,ylim=c(0,1000),xlim=c(xMin1,xMax),ylab="total cases confirmed in Wuhan",col="white")
  # 
  # polygon(c(date_range,rev(date_range)),c(cumsum(Rep_local_quantile[2,]),rev(cumsum(Rep_local_quantile[4,]))),lty=0,col=rgb(0,0.3,1,0.35))
  # polygon(c(date_range,rev(date_range)),c(cumsum(Rep_local_quantile[1,]),rev(cumsum(Rep_local_quantile[5,]))),lty=0,col=rgb(0,0.3,1,0.2))
  # lines(date_range,cumsum(Rep_local_quantile[3,]),type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  # points(date_range,case_data_wuhan_conf_time)
  # 
  # text(labels="model estimate",x=xMin1,y=0.9*1e3,adj=0,col="blue")
  # text(labels="non-fitted data (used for validation)",x=xMin1,y=0.8*1e3,adj=0,col="black")
  # 
  # lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")
  
  # Plot international cases confirmed
  ym1 <- 12
  plot(date_range,case_time,pch=19,ylim=c(0,ym1),xlim=c(xMin1,xMax),ylab="new international exports confirmed",col="white")
  
  polygon(c(date_range,rev(date_range)),c(Rep_quantile[2,],rev(Rep_quantile[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(Rep_quantile[1,],rev(Rep_quantile[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,Rep_quantile[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  points(date_range,case_time,pch=1)
  
  text(labels="model estimate",x=xMin1,y=0.95*ym1,adj=0,col="blue")
  text(labels="non-fitted data (used for validation)",x=xMin1,y=0.85*ym1,adj=0,col="black")
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")
  title(LETTERS[letR],adj=0); letR = letR +1
  
  # Plot international confirmations
  
  plot_international(Rep_plot) #- turn off?
  #plot.new()
  title(LETTERS[letR],adj=0); letR = letR +1
  
  
  # Plot estimated local infections
  yMax <- 0.15
  
  plot(date_rangeA,Inf_quantileA[1,],col="white",ylim=c(0,1.2*yMax),xlim=c(xMin1,xMax),xlab="",ylab="propn prevalence in Wuhan")
  polygon(c(date_rangeA,rev(date_rangeA)),c(Inf_quantileA[2,],rev(Inf_quantileA[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_rangeA,rev(date_rangeA)),c(Inf_quantileA[1,],rev(Inf_quantileA[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_rangeA,Inf_quantileA[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  
  # Japan flight 1
  CI_flight_japan_1 <- bin_conf(prop_flight_1_japan[1],prop_flight_1_japan[2])
  points(date_flights_out_1_japan,CI_flight_japan_1[1],pch=19)
  lines(c(date_flights_out_1_japan,date_flights_out_1_japan),c(CI_flight_japan_1[2],CI_flight_japan_1[3]))
  
  # Japan flight 2
  CI_flight_japan_2 <- bin_conf(prop_flight_2_japan[1],prop_flight_2_japan[2])
  points(date_flights_out_2_japan,CI_flight_japan_2[1],pch=19)
  lines(c(date_flights_out_2_japan,date_flights_out_2_japan),c(CI_flight_japan_2[2],CI_flight_japan_2[3]))
  
  # Japan flight 3
  CI_flight_japan_3 <- bin_conf(prop_flight_3_japan[1],prop_flight_3_japan[2])
  points(date_flights_out_3_japan,CI_flight_japan_3[1],pch=19)
  lines(c(date_flights_out_3_japan,date_flights_out_3_japan),c(CI_flight_japan_3[2],CI_flight_japan_3[3]))
  
  # Germany flight
  CI_flight <- bin_conf(prop_flight_2_germany[1],prop_flight_2_germany[2])
  points(date_flights_out_2_germany,CI_flight[1],pch=19)
  lines(c(date_flights_out_2_germany,date_flights_out_2_germany),c(CI_flight[2],CI_flight[3]))
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")
  text(labels="model estimate",x=xMin1,y=1.1*yMax,adj=0,col="blue")
  text(labels="fitted data from evacuation flights",x=xMin1,y=1*yMax,adj=0,col="black")
  
  title(LETTERS[letR],adj=0); letR = letR + 1
  
  # Plot susceptibles
# 
#   plot(date_rangeA,S_quantileA[1,],col="white",ylim=c(0,1),xlim=c(xMin1,xMax),xlab="",ylab="propn prevalence in Wuhan")
#   polygon(c(date_rangeA,rev(date_rangeA)),c(S_quantileA[2,],rev(S_quantileA[4,])),lty=0,col=rgb(0,0.3,1,0.35))
#   polygon(c(date_rangeA,rev(date_rangeA)),c(S_quantileA[1,],rev(S_quantileA[5,])),lty=0,col=rgb(0,0.3,1,0.2))
#   lines(date_rangeA,S_quantileA[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
#   
#   text(labels="model estimate",x=min(date_rangeA),y=9,adj=0,col="blue")
#   
#   lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,10),col="red")
#   title(LETTERS[letR],adj=0); letR = letR + 1
  
  # Plot case total predictions
  xMin1 <- min(as.Date("2020-01-23"))
  ym1 <- 5000
  plot(date_range,Rep_local_quantile[1,],col="white",ylim=c(0,ym1),xlim=c(xMin1,xMax),xlab="",ylab="new confirmed cases in Wuhan")
  polygon(c(date_range,rev(date_range)),c(Rep_local_quantile[2,],rev(Rep_local_quantile[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(Rep_local_quantile[1,],rev(Rep_local_quantile[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,Rep_local_quantile[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,1e6),col="red")
  text(labels="travel restrictions",x=wuhan_travel_restrictions+0.5,y=0.9*ym1,adj=0,col="red")
  
  points(cases_Wuhan$date,cases_Wuhan$new_case,pch=19)
  
  text(labels="model estimate",x=xMin1+1,y=0.65*ym1,adj=0,col="blue")
  text(labels="fitted data",x=xMin1+1,y=0.55*ym1,adj=0,col="black")
  
  title(LETTERS[letR],adj=0); letR = letR + 1
  
  # Plot reproduction number
  
  date_rangeB <- date_rangeA[date_rangeA>as.Date("2019-12-15")]
  R0_quantileB <- R0_quantileA[,date_rangeA>as.Date("2019-12-15")]
  xMax1 <- as.Date("2020-02-01") #xMax #xMax #
  
  plot(date_rangeB,R0_quantileB[1,],col="white",ylim=c(0,10),xlim=c(min(date_rangeB),xMax1),xlab="",ylab=expression(paste(R[t])))
  
  polygon(c(date_rangeB,rev(date_rangeB)),c(R0_quantileB[2,],rev(R0_quantileB[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_rangeB,rev(date_rangeB)),c(R0_quantileB[1,],rev(R0_quantileB[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_rangeB,R0_quantileB[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  lines(date_rangeB,1+0*R0_quantileB[3,],lty=2)
  
  text(labels="model estimate",x=min(date_rangeB),y=9,adj=0,col="blue")
  
  lines(c(wuhan_travel_restrictions,wuhan_travel_restrictions),c(0,10),col="red")
  title(LETTERS[letR],adj=0); letR = letR + 1
    
  # output figure
  dev.copy(png,paste("plots/cases_inference.png",sep=""),units="cm",width=20,height=20,res=150)
  #dev.copy(pdf,paste("plots/cases_inference.pdf",sep=""),width=8,height=8)
  dev.off()
  
  
}

# Plot international cumulative --------------------------------------------

plot_international <- function(Rep_plot){
  
  # filename="1"
  
  #load(paste0("outputs/bootstrap_fit_",filename,".RData"))
  
  # Get confirmed cases
  Rep_plot_total <- colSums(Rep_plot)
  #rep_CI <- apply(Rep_plot_total,1,function(x){quantile(x,c(0.025,0.5,0.975))})
  rep_CI <- quantile(Rep_plot_total,c(0.025,0.5,0.975))
  
  x_scaled <- rep_CI[2]
  x_expected <- sapply(x_scaled,function(x){x*as.numeric(top_risk$risk)}); #x_expected <- t(x_expected) # expected exported cases in each location
  
  # Tally country cases
  
  unique_country <- unique(case_data$country)
  
  case_export_vector <- rep(0,nrow(top_risk))
  for(ii in 1:nrow(top_risk)){
    
    match_ID <- match(top_risk[ii,]$label,unique_country)
    
    if(!is.na(match_ID)){
      total_exports <- case_data %>% filter(country==unique_country[match_ID]) %>% select(number) %>% sum()
    }else{
      total_exports <- 0
    }
    
    case_export_vector[ii] <- total_exports
  }

  
  # Plot expected vs observed
  #par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(2,0.7,0))
  par(mar=c(3,3,1,1))
  
  ymax1 <- 11
  plot(x_expected,case_export_vector,col="white",ylim=c(0,ymax1),xlim=c(0,ymax1),xlab="expected international exports from Wuhan",ylab="confirmed international exports")
  lines(c(-10,20),c(-10,20),lty=2)
  points(x_expected,case_export_vector,pch=1,col="blue",cex=0.7) #-- REMOVED TEMPORARILY
  #text(labels="Australia",x=1.2,y=6,col="blue",adj=0)
  #text(labels="USA",x=1.6,y=5,col="blue",adj=0)
  
  #text(labels="France",x=0.4,y=3,col="blue",adj=0)
  #text(labels="Thailand",x=3.9,y=3.7,col="blue",adj=0)
  
  par(mar=c(2,3,1,1))
  
  # dev.copy(png,paste("plots/export_plot.png",sep=""),units="cm",width=15,height=10,res=150)
  # dev.off()
  # 
}

# Plot dispersion --------------------------------------------

plot_dispersion <- function(filename="1"){
  
  # filename="1"
  
  load(paste0("outputs/bootstrap_fit_",filename,".RData"))
  
  # Extract R credible interval
  period_interest <- as.Date(c("2020-01-01","2020-01-23"))
  
  par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(2,0.7,0))
  
    index_pick <- match(period_interest,date_range)
    R0_all <- R0_plot[index_pick[1]:index_pick[2],]
    #dim(R0_all) <- NULL # collapse data
    
    med_R0 <- apply(R0_all,2,function(x){quantile(x,0.5)})
    
    R0_CrI <- quantile(med_R0,c(0.025,0.25,0.5,0.75,0.975))
    
    MERS_k <- 0.26
    SARS_k <- 0.16
    k_seq <- seq(0.01,0.5,0.01)
  
    # Outbreak calcs
    R0_med <- 1-sapply(k_seq,function(x){numerical_solver(R0_CrI[3],x)})
    R0_CrI_1 <- 1-sapply(k_seq,function(x){numerical_solver(R0_CrI[1],x)})
    R0_CrI_2 <- 1-sapply(k_seq,function(x){numerical_solver(R0_CrI[5],x)})
    R0_CrI_1_50 <- 1-sapply(k_seq,function(x){numerical_solver(R0_CrI[2],x)})
    R0_CrI_2_50 <- 1-sapply(k_seq,function(x){numerical_solver(R0_CrI[4],x)})
    
    # Estimate range
    c(R0_med[k_seq==SARS_k],R0_med[k_seq==MERS_k])
    
    # Plot results
    plot(k_seq,R0_med,type="l",ylim=c(0,1),xlab=c("extent of homogeneity in transmission"),ylab="probability of large outbreak with single case",col="white",xaxs="i",yaxs="i")
    polygon(c(k_seq,rev(k_seq)),c(R0_CrI_1,rev(R0_CrI_2)),lty=0,col=rgb(0,0.3,1,0.35))
    polygon(c(k_seq,rev(k_seq)),c(R0_CrI_1_50,rev(R0_CrI_2_50)),lty=0,col=rgb(0,0.3,1,0.35))
    lines(k_seq,R0_med,col="blue")
    
    lines(c(MERS_k,MERS_k),c(-1,10),lty=2); text(labels="MERS-CoV",x=1.02*MERS_k,y=0.7,adj=0,col="black")
    lines(c(SARS_k,SARS_k),c(-1,10),lty=2); text(labels="SARS",x=1.02*SARS_k,y=0.6,adj=0,col="black")
    title(LETTERS[1],adj=0)
    
    n_seq <- seq(0,10,1)
    ext_m <- 1-(1-R0_med[k_seq==SARS_k])^n_seq
    ext_1 <- 1-(1-R0_CrI_1[k_seq==SARS_k])^n_seq
    ext_2 <- 1-(1-R0_CrI_2[k_seq==SARS_k])^n_seq
    ext_11 <- 1-(1-R0_CrI_1_50[k_seq==SARS_k])^n_seq
    ext_22 <- 1-(1-R0_CrI_2_50[k_seq==SARS_k])^n_seq
    
    plot(n_seq,ext_m,type="l",ylim=c(0,1),xlim=c(-0.5,10.5),xlab=c("number of introductions"),ylab="probability of large outbreak",col="white",xaxs="i",yaxs="i")
    points(n_seq,ext_m,col="blue",pch=19)
    for(ii in 1:length(n_seq)){
      lines(c(n_seq[ii],n_seq[ii]),c(ext_1[ii],ext_2[ii]),col="blue")
    }
    
    
    title(LETTERS[2],adj=0)


  
  
  dev.copy(png,paste("plots/calc_1.png",sep=""),units="cm",width=20,height=10,res=150)
  dev.off()
  
  
}

# Output R0 estimates over time --------------------------------------------

r0_value_output <- function(filename="1"){
  
  # filename="1"
  
  load(paste0("outputs/bootstrap_fit_",filename,".RData"))
  
  # Extract R0 estimates

  period_interest <- as.Date("2020-01-15")
  
  R0_plot = R0_plot[,!is.na(R0_plot[t_period,])]

  med_R0 <- apply(R0_plot,1,function(x){quantile(x,c(0.5))})
  c.text(med_R0[match(period_interest,date_range)])
  
  #file_out <- as_tibble(cbind(date_range,c.text(t(med_R0))))
  
  # Median R0 range
  period_interest <- as.Date(c("2019-12-15","2020-01-23"))
  index_pick <- match(period_interest,date_range)
  R0_all <- R0_plot[index_pick[1]:index_pick[2],]
  
  med_R0 <- apply(R0_all,1,function(x){quantile(x,c(0.5))})
  
  c.text(med_R0)
  
  # - - - 
  # symptomatic cases in Wuhan
  
  period_interest <- as.Date("2020-01-23")
  
  c.text(I_plot[date_range==period_interest,])

  
}

# Exporter case scaling function --------------------------------------------

fit_int_cases <- function(x_val){
  
  x_scaled <- x_val #* passengers_daily / wuhan_area
  
  x_expected <- sapply(x_scaled,function(x){x*as.numeric(top_risk$risk)}) %>% t() # expected exported cases in each location
  x_lambda <- rowSums(x_expected)
  #x_lambda
  rpois(length(x_lambda),lambda=x_lambda)
  
}

# Calculate new oubtreak probability --------------------------------------

numerical_solver <- function(r0, k){
  
  fun <- function (s) {(1 + (r0/k)*(1 - s))^(-k) - s}
  solutions <- rootSolve::multiroot(fun, c(0, 1))$root
  
  realistic_sol <- min(solutions)
  return(realistic_sol)
  
}
