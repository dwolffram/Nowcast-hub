source("RESPINOW/initialize.R")

forecast_dates <- seq(from = as.Date("2022-11-06"),
                      to = as.Date("2023-05-07"),
                      by = 7)

for (d in as.list(forecast_dates)) {
  print(d)
  tryCatch(compute_nowcast("seasonal_influenza", d),
           error = function(e) {print(paste0("Error for date: ", d, "."))})
}

for (d in as.list(forecast_dates)) {
  print(d)
  tryCatch(compute_nowcast("rsv_infection", d),
           error = function(e) {print(paste0("Error for date: ", d, "."))})
}

compute_nowcast("seasonal_influenza", "2023-04-02")
