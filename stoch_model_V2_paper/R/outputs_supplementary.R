
# Supplementary plots --------------------------------------------------------------

# - - - - - - - - - - 
# 1D parameter grid search for relative reporting outside Wuhan
# 
# par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(2,0.7,0))
# 
# s_out <- MLE_check(p_name = "local_rep_prop", theta_tab=seq(0.001,0.05,0.001),nn=2e3)
# 
# s_out[max(s_out$lik)==s_out$lik,] # maximum likelihood
# 
# plot(s_out$param,s_out$lik,ylim=c(max(s_out$lik)-5,max(s_out$lik)),xlab="value",ylab="log likelihood",pch=19)
# lines(c(0,1e2),c(1,1)*(max(s_out$lik)-1.92),lty=2)
# 
# dev.copy(png,paste("plots/param_rel_1.png",sep=""),units="cm",width=10,height=10,res=150)
# dev.off()
# 
# # - - - - - - - - - - 
# # 1D parameter grid search for beta volatility (DEPRECATED)
# 
# par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(2,0.7,0))
# 
# s_out <- MLE_check(p_name = "betavol", theta_tab=seq(0.01,1,0.01),nn=3e3)
# 
# plot(s_out$param,s_out$lik,ylim=c(max(s_out$lik)-5,max(s_out$lik)),xlab="value",ylab="log likelihood",pch=19)
# lines(c(0,1e2),c(1,1)*(max(s_out$lik)-1.92),lty=2)
# 
# dev.copy(png,paste("plots/param_vol_1.png",sep=""),units="cm",width=10,height=10,res=150)
# dev.off()
# 
# 
# # - - - - - - - - - - 
# # 2D parameter grid search for proportion confirmed and relative reporting outside Wuhan
# 
# # par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(2,0.7,0))
# 
# MLE_check_2D(p1_name = "local_rep_prop", p2_name = "confirmed_prop",
#                theta_tab1 = seq(0.002,0.03,0.002), theta_tab2 = seq(0.5,1,0.025),nn=1e3)
# 
# profile_plot(p1_name = "local_rep_prop", p2_name = "confirmed_prop", filename=1)
# 
# 

# # 2D parameter grid search for proportion confirmed, relative reporting outside Wuhan, and transmission vo

MLE_check_3D(p1_name = "local_rep_prop", 
             p2_name = "confirmed_prop",
             p3_name = "betavol",
             theta_tab1 = seq(0.001,0.015,0.001), 
             theta_tab2 = seq(0.5,1,0.05),
             theta_tab3 = seq(0.1,0.6,0.05),
             nn=2e3,
             filename=1)


profile_plot(p1_name = "local_rep_prop", 
             p2_name = "confirmed_prop", 
             p3_name = "betavol",
             filename=1)








