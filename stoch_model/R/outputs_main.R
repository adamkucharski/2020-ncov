# Run analysis--------------------------------------------------------------



# - - -
# Run bootstrap SMC 
aa <- Sys.time()

run_fits(rep_plot=100, # number of repeats
         nn=2e3, #number of particles
         dt=t_step,
         filename="1"
)

aa2 <- Sys.time() - aa

# Output plots
plot_outputs(filename="1")

#plot_dispersion(filename="1")


# Run models --------------------------------------------------------------



