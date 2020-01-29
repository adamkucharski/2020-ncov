# Run analysis--------------------------------------------------------------



# - - -
# Run bootstrap SMC 
run_fits(rep_plot=200, # number of repeats
         nn=2e3,#number of particles
         dt=0.25,
         filename="1"
)

# Output plots
plot_outputs(filename="1")


# Run models --------------------------------------------------------------



