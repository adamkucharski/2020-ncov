# Run analysis--------------------------------------------------------------



# - - -
# Run bootstrap SMC 
run_fits(rep_plot=100, # number of repeats
         nn=2e3, #number of particles
         dt=0.25,
         filename="1"
)

# Output plots
plot_outputs(filename="1") # Figure 2

# plot_dispersion(filename="1") # Figure 3


# Run models --------------------------------------------------------------



