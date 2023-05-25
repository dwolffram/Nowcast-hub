source(file = "RESPINOW/initialize.R")
forecast_date <- as.Date("2023-03-12")

# Show progress
message("Import raw data")

# Download raw reporting data from nowcast hub
reporting_data_raw <- read_csv(
  file = "https://raw.githubusercontent.com/KITmetricslab/RESPINOW-Data/main/data/Survstat/seasonal_influenza_reporting_triangle_survstat_preprocessed.csv",
  show_col_types = FALSE)

reporting_data_raw <- load_survstat("seasonal_influenza", forecast_date)

message("Tidy data")

reporting_data <- reporting_data_raw |>
  mutate(
    # Convert characters location and age_group to factors
    # This to ensure that the order is the same as in the raw data
    location = location |> fct_inorder(),
    age_group = age_group |> fct_inorder()) |>
  # Rename column value_>80d to value_81d
  rename(
    "value_11w" = "value_>10w") |>
  # Convert to tidy format
  # The delay comes from the number in the value_{**}d columns
  # The true number of hospitalisations by date and delay is called n_true
  pivot_longer(
    cols = starts_with("value"),
    names_to = "delay",
    names_pattern = "value_(\\d{1,2})w",
    names_transform = list(delay = as.integer),
    values_to = "n_true",
    values_transform = list(n_true = as.integer),
    values_drop_na = TRUE) |>
  # Relocate date before delay
  relocate(
    date, .before = delay) |>
  # Arrange records
  arrange(
    location, age_group, date, delay)


message("Set maximum delay")

# In the raw data the maximum delay is 81 days (see script tidy_data.R)
# where delay = 81 contains all delays >80 days
# Here we set a maximum delay to a given value: max_delay (see script initialize.R)
# and aggregate all delays >= max_delay into delay = max_delay
reporting_data <- reporting_data |>
  mutate(
    delay = if_else(
      condition = delay >= max_delay,
      true = max_delay,
      false = delay)) |>
  group_by(
    location, age_group, date, delay) |>
  summarise(
    n_true = sum(n_true),
    .groups = "drop")

# Complete reporting matrix by date and delay
# This sets n_true to NA outside the reporting triangle, as in the raw data
# We need this explicit NA's in function make_nowcast()
reporting_data <- reporting_data |>
  complete(
    location, age_group, date, delay)


message("Truncate reporting triangle")

# The reporting triangle is truncated up to forecast_date
# In a real-time setting this happens naturally
# n_rep is set to NA if date + delay > forecast_date
# If not, n_rep is equal to n_true
# This forced trunction is useful for evaluating the nowcast retrospectively
reporting_data <- reporting_data |>
  mutate(
    n_rep = if_else(
      condition = date + 7*delay > forecast_date,
      true = NA_integer_,
      false = n_true))

message("Make data for fitting")

reporting_fit_data <- reporting_data |>
  filter(
    # Only use dates back to 2*max(delay) from forecast_date
    date > forecast_date - 2*7*max_delay & date <= forecast_date) |>
  mutate(
    # Get the weekday of date and weekday of date_report
    date_report = date + 7*delay,
    # These variables are used in the model instead of date and delay:
    # date_trans is the number days since min(date)
    # delay_trans is the sqrt of delay, because this
    # almost results in a log-linear delay effect
    date_trans = as.numeric(date - min(date))/7,
    delay_trans = delay |> sqrt(),
    # I_delay_max is an indicator variable for the maximum delay category
    I_delay_max = (delay == max(delay)) |>
      as.integer() |>
      factor())

message("Split fit data")

# The dataset for fitting is split into three parts
# Each split gets its own model and nowcast
# Between the parentheses is the number of groups
# Excluded split: location not DE & age_group not 00+ (16 x 6 = 96)

# 1. location DE & agegroup 00+ (1 x 1 = 1)
reporting_fit_data_DE_00 <- reporting_fit_data |>
  filter(
    location == "DE" & age_group == "00+") |>
  droplevels()

# 2. location DE & age_group not 00+ (1 x 6 = 6)
reporting_fit_data_DE_not00 <- reporting_fit_data |>
  filter(
    location == "DE" & age_group != "00+") |>
  droplevels()

# 3. location not DE & age_group 00+ (16 x 1 = 16)
reporting_fit_data_notDE_00 <- reporting_fit_data |>
  filter(
    location != "DE" & age_group == "00+") |>
  droplevels()


# Show progress
message("Fit gam model DE 00")

reporting_gam_DE_00 <- bam(
  formula = n_rep ~
    s(date_trans, bs = "ps", k = 5) +
    s(delay_trans, bs = "ps", k = 5) +
    # s(day, bs = "re") +
    # s(day_report, bs = "re") +
    s(I_delay_max, bs = "re"),
  family = nb,
  data = reporting_fit_data_DE_00,
  discrete = TRUE,
  select = TRUE)

message("Make nowcast DE 00")

nowcast_data_DE_00 <- make_nowcast(
  fitted_gam = reporting_gam_DE_00,
  fit_data = reporting_fit_data_DE_00,
  forecast_date = forecast_date)


message("Fit gam model DE not 00")

reporting_gam_DE_not00 <- bam(
  formula = n_rep ~
    s(date_trans, bs = "ps", k = 5) +
    s(delay_trans, bs = "ps", k = 5) +
    # s(day, bs = "re") +
    # s(day_report, bs = "re") +
    s(I_delay_max, bs = "re") +
    s(age_group, bs = "re") +
    ti(date_trans, age_group, bs = c("ps", "re"), k = c(5, 6)) +
    ti(delay_trans, age_group, bs = c("ps", "re"), k = c(5, 6)) +
    # ti(day, age_group, bs = "re") +
    # ti(day_report, age_group, bs = "re") +
    ti(I_delay_max, age_group, bs = "re"),
  family = nb,
  data = reporting_fit_data_DE_not00,
  discrete = TRUE,
  select = TRUE)

message("Make nowcast DE not 00")

nowcast_data_DE_not00 <- make_nowcast(
  fitted_gam = reporting_gam_DE_not00,
  fit_data = reporting_fit_data_DE_not00,
  forecast_date = forecast_date)



message("Fit gam model not DE 00")

reporting_gam_notDE_00 <- bam(
  formula = n_rep ~
    s(date_trans, bs = "ps", k = 5) +
    s(delay_trans, bs = "ps", k = 5) +
    s(I_delay_max, bs = "re") +
    s(location, bs = "re") +
    ti(date_trans, location, bs = c("ps", "re"), k = c(5, 16)) +
    ti(delay_trans, location, bs = c("ps", "re"), k = c(5, 16)) +
    ti(I_delay_max, location, bs = "re"),
  family = nb,
  data = reporting_fit_data_notDE_00,
  discrete = TRUE,
  select = TRUE)

message("Make nowcast not DE 00")

nowcast_data_notDE_00 <- make_nowcast(
  fitted_gam = reporting_gam_notDE_00,
  fit_data = reporting_fit_data_notDE_00,
  forecast_date = forecast_date)


message("Combine nowcasts")

# We assume nowcast_data_DE_00 is always there
nowcast_data <- nowcast_data_DE_00

# The other two are optional
if (exists("nowcast_data_DE_not00")) {
  nowcast_data <- bind_rows(
    nowcast_data,
    nowcast_data_DE_not00)
}
if (exists("nowcast_data_notDE_00")) {
  nowcast_data <- bind_rows(
    nowcast_data,
    nowcast_data_notDE_00)
}

message("Export csv file")

# Specifications:
# https://github.com/KITmetricslab/hospitalization-nowcast-hub/wiki/Data-format

# Filename for given forecast_date
filename <- str_glue("RESPINOW/nowcasts_retrospective/{forecast_date}-spline_smoothing.csv")

# If the csv for the given forecast_date does not exists, create one

nowcast_data |>
  rename(
    target_end_date = date) |>
  filter(
    target_end_date >= forecast_date - 7*max_delay) |>
  mutate(
    forecast_date = forecast_date,
    target = str_glue("{(target_end_date - forecast_date)/7} week ahead inc hosp"),
    pathogen = "seasonal_influenza") |>
  pivot_longer(
    cols = starts_with("N_"),
    names_sep = "_",
    names_to = c(NA, "type", "quantile"),
    values_to = "value") |>
  select(
    location, age_group, forecast_date, target_end_date, target, type, quantile, value, pathogen) |>
  write_csv(
    file = filename,
    quote = "none")
