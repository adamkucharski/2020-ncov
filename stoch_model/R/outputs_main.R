# Run analysis--------------------------------------------------------------



# - - -
# Run bootstrap SMC 
run_fits(rep_plot=1000, # number of repeats
         nn=1e3,#number of particles
         dt=0.25,
         filename="1"
)

# Output plots
plot_outputs(filename="1")


# Run models --------------------------------------------------------------



