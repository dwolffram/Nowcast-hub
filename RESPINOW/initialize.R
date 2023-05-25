#
# Initialization
#

# Show progress
message("Initialization")

# Load packages
library(tidyverse)
library(lubridate)
library(colorspace)
library(mgcv)

# Source functions
# walk(
#   .x = list.files(path = "Functions", full.names = TRUE),
#   .f = source)

source("RESPINOW/make_nowcast.R")
source("RESPINOW/repair_quantile_names.R")
source("RESPINOW/load_data.R")
source("RESPINOW/compute_nowcast.R")


# Set time locale to English
Sys.setlocale("LC_ALL", "C")

# Weeks start on Monday
options(lubridate.week.start = 1)

# Set maximum delay
max_delay <- 4L

# Number of Monte Carlo simulations for nowcast
n_sim <- 1000

# Set forecast_date
# forecast_date <- today()
# forecast_date <- floor_date(forecast_date, unit = "weeks", week_start = 7)
