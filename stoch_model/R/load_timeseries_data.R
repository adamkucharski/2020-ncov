# Timeseries data

# define values

start_date <- as.Date("2019-12-10")
end_date <- as.Date("2020-01-25")
date_range <- seq(start_date,end_date,1)

# load data
case_data$export_probability <- as.numeric(travel_data[match(case_data$location,travel_data$label),]$risk) # Add risk
case_data <- case_data[!is.na(case_data$export_probability),] # Only use available data

# tally cases
case_time <- rep(0,length(date_range))

for(ii in 1:length(date_range)){
  case_time[ii] = sum(case_data[case_data$date==date_range[ii],]$cases)
}

# shift data into weeks
t_period <- as.numeric(end_date-start_date)+1
case_data <- case_data %>% mutate(time = as.numeric(date - start_date + 1))

# compile matrix of cases in top 30 risk locations
n_risk <- 30
top_risk <- travel_data[1:n_risk,]
case_data_matrix <- matrix(0,nrow=t_period,ncol=n_risk)
match_list_cases <- match(case_data$location,top_risk$label)
case_data_matrix[case_data$time,match_list_cases] <- case_data$cases # add detected cases




