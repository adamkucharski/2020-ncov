
# Supplementary plots --------------------------------------------------------------

# - - - - - - - - - - 
# Parameter grid search for relative reporting outside Wuhan

par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(2,0.7,0))

s_out <- MLE_check(p_name = "local_rep_prop", theta_tab=seq(0.001,0.01,0.001),nn=1e3)

plot(s_out$param,s_out$lik,ylim=c(max(s_out$lik)-10,max(s_out$lik)),xlab="value",ylab="log likelihood",pch=19)
lines(c(0,1e2),c(1,1)*(max(s_out$lik)-1.92),lty=2)

# - - - - - - - - - - 
# Parameter grid search for beta volatility

par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(2,0.7,0))

s_out <- MLE_check(p_name = "betavol", theta_tab=seq(0.1,1,0.1),nn=1e3)

plot(s_out$param,s_out$lik,ylim=c(max(s_out$lik)-10,max(s_out$lik)),xlab="value",ylab="log likelihood",pch=19)
lines(c(0,1e2),c(1,1)*(max(s_out$lik)-1.92),lty=2)

# - - - - - - - - - - 
# Parameter grid search for effective population size

par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(2,0.7,0))

s_out <- MLE_check(p_name = "pop_travel", theta_tab=seq(800e5,1400e5,100e5),nn=1e3)

plot(s_out$param,s_out$lik,ylim=c(max(s_out$lik)-10,max(s_out$lik)),xlab="value",ylab="log likelihood",pch=19)
lines(c(0,1e2),c(1,1)*(max(s_out$lik)-1.92),lty=2)



# Sensitivity analysis --------------------------------------------------------------




#  Distributions of incubation and infectious periods and reporting.

#  Sensitivity to initial conditions and data used

#  Model fits for individual countries

# Model fits with fixed transmission rate

#  Model fits for individual countries
