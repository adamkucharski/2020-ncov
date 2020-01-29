# Timeseries data

# Define values
omit_recent <- 5
omit_conf <- 0

start_date <- as.Date("2019-11-15")
end_date <- max(case_data_in$date) # omit recent day?
date_range <- seq(start_date,end_date,1)

# When restrictions started
wuhan_travel_restrictions <- as.Date("2020-01-23")
wuhan_travel_time <- as.numeric(wuhan_travel_restrictions - start_date + 1)


# Load international confirmation data --------------------------------------------
case_data <- case_data_in
cutoff_case_int <- max(case_data$date) - omit_conf # omit final days of time points

case_data$export_probability <- as.numeric(travel_data[match(case_data$country,travel_data$label),]$risk) # Add risk
case_data <- case_data[!is.na(case_data$export_probability),] # Only use available data

# tally cases
case_time <- rep(0,length(date_range))

for(ii in 1:length(date_range)){
  case_time[ii] = sum(case_data[case_data$date==date_range[ii],]$number)
}

case_time[date_range>cutoff_case_int] <- NA # NOTE not used

# shift data into weeks
t_period <- as.numeric(end_date-start_date)+1
case_data <- case_data %>% mutate(time = as.numeric(date - start_date + 1))

# compile matrix of cases in top 30 risk locations
n_risk <- 20
top_risk <- travel_data[1:n_risk,]

# Calculate exports by country
case_data_matrix <- matrix(0,nrow=t_period,ncol=n_risk)
match_list_cases <- match(case_data$country,top_risk$label)
for(ii in 1:nrow(case_data)){
  case_data_matrix[case_data[ii,]$time,match_list_cases[ii]] <- case_data[ii,]$number # add detected cases
}

# Load international onset data --------------------------------------------

case_data_onset <- international_onset_data_in
cutoff_time_int_onsets <- max(case_data_onset$date) - omit_recent # omit final days of time points

case_data_onset[case_data_onset$date>cutoff_time_int_onsets,"number"] <- NA
case_data_onset_time <- rep(0,length(date_range))

for(ii in 1:length(date_range)){
    case_data_onset_time[ii] <- sum(case_data_onset[case_data_onset$date==date_range[ii],]$number)
}

case_data_onset_time[date_range>cutoff_time_int_onsets] <- NA # omit final points


# Load China onset data --------------------------------------------

case_data_china <- china_onset_data_in
cutoff_time_china <- max(case_data_china$date) - omit_recent # omit final days of time points
case_data_china[case_data_china$date>cutoff_time_china,"number"] <- NA
case_data_china_time <- rep(0,length(date_range))

# ensure final points are omitted
for(ii in 1:length(date_range)){
  case_data_china_time[ii] <- sum(case_data_china[case_data_china$date==date_range[ii],]$number)
}
case_data_china_time[date_range>cutoff_time_china] <- NA # omit final points

# Load Wuhan early data --------------------------------------------

case_data_wuhan <- wuhan_onset_data_in
case_data_wuhan$number <- case_data_wuhan$number - case_data_wuhan$number_market # calculate non-market exposures

case_data_wuhan_time <- rep(0,length(date_range))

for(ii in 1:length(date_range)){
  case_data_wuhan_time[ii] = sum(case_data_wuhan[case_data_wuhan$date==date_range[ii],]$number)
}

# Quick plot
#plot(case_data_wuhan$date,case_data_wuhan$number,xlim=as.Date(c("2019-12-01","2020-01-25")),ylim=c(0,30)); points(case_data_china$date,case_data_china$number,col="blue")

# INITIAL APPROXIMATION -  ADD TOGETHER CHINA TIMESERIES

case_data_china_time <- case_data_china_time + case_data_wuhan_time

# Load Wuhan confirmed data --------------------------------------------

case_data_wuhan_conf <- wuhan_conf_data_in
case_data_wuhan_conf_time <- rep(0,length(date_range))
cutoff_time_wuhan <- max(case_data_wuhan_conf$date)

for(ii in 1:length(date_range)){
  case_data_wuhan_conf_time[ii] = sum(case_data_wuhan_conf[case_data_wuhan_conf$date==date_range[ii],]$number)
}

case_data_wuhan_conf_time[date_range>cutoff_time_wuhan] <- NA # omit all but single point

# Compile list of data to use:


data_list = list(local_case_data_onset = case_data_china_time, 
                 local_case_data_conf = case_data_wuhan_conf_time,
                 int_case_onset = case_data_onset_time,
                 int_case_conf = case_data_matrix)

#data_list = list(local_case_data_tt=case_data_china_time[tt],case_data_tt=case_data_onset_time[tt],rep_data_tt=case_data_matrix[tt,])




