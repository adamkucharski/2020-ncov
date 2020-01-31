
# Supplementary plots --------------------------------------------------------------

# - - - - - - - - - - 
# Parameter grid search for relative reporting outside Wuhan

par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(2,0.7,0))

s_out <- MLE_check(p_name = "local_rep_prop", theta_tab=seq(0.001,0.05,0.001),nn=2e3)

s_out[max(s_out$lik)==s_out$lik,] # maximum likelihood

plot(s_out$param,s_out$lik,ylim=c(max(s_out$lik)-5,max(s_out$lik)),xlab="value",ylab="log likelihood",pch=19)
lines(c(0,1e2),c(1,1)*(max(s_out$lik)-1.92),lty=2)

dev.copy(png,paste("plots/param_rel_1.png",sep=""),units="cm",width=10,height=10,res=150)
dev.off()

# - - - - - - - - - - 
# Parameter grid search for beta volatility

par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(2,0.7,0))

s_out <- MLE_check(p_name = "betavol", theta_tab=seq(0.01,1,0.01),nn=3e3)

plot(s_out$param,s_out$lik,ylim=c(max(s_out$lik)-5,max(s_out$lik)),xlab="value",ylab="log likelihood",pch=19)
lines(c(0,1e2),c(1,1)*(max(s_out$lik)-1.92),lty=2)

dev.copy(png,paste("plots/param_vol_1.png",sep=""),units="cm",width=10,height=10,res=150)
dev.off()



