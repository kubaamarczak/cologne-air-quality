## ##-##-##-##-##-##-##-##-##-##-##-##-##-##-##-##-##-##-##-##-##-##-##-##-##-##
## statistical analysis of air quality in cologne with traffic and background
## stations.

## 4 stations (background: Chorweiler, Rodenkirchen;
##             traffic: Clevischer Ring, Turiner Straße)

## analysis based on https://luftdaten.umweltbundesamt.de
## -##-##-##-##-##-##-##-##-##-##-##-##-##-##-##-##-##-##-


## -----------------------------------------------------------------------------
## libraries
## ---------

library("arrow")
library("here")
library("dplyr")
library("lubridate")
library("mgcv")
library("gratia")
library("ggplot2")
library("patchwork")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## read and manipulate parquet file from fetch
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## read
df <- arrow::read_parquet(
  here::here(
    "data",
    "processed",
    "cologne_measures_2016_2024.parquet"
  )
)

## set timezone
df$date_start <- lubridate::with_tz(df$date_start, "UTC")

## station codes
stations <- c(1108, 1114, 1254, 1255)

## add necessary columns for analysis
df <- df |>
  dplyr::mutate(
    station_type = factor(station,
      levels = stations,
      labels = c(
        "background", "background",
        "traffic", "traffic"
      )
    ),
    hour = hour(date_start),
    weekend = wday(date_start, week_start = 1) >= 6,
    month = month(date_start),
    year = year(date_start),
    time_index = as.numeric(difftime(date_start,
      min(date_start),
      units = "days"
    ))
  ) |>
  dplyr::filter(!is.na(value))

## peek in the data frame
df

## component 1: PM₁₀, component 5: NO₂; scope 2: 1SMW

## filter for NO₂
df_no2 <- df |>
  dplyr::filter(component == 5)

## data frame to work with
df_no2

## limiting on NO₂ for this analysis

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## model building
## ~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## tweedie distribution
## ~~~~~~~~~~~~~~~~~~~~

## -----------------------------------------------------------
## function for diagnostics
## ------------------------

diagnostics <- function(model) {
  ## summary
  print(summary(model))

  ## diagnostic plots
  plot(model, pages = 1, residuals = TRUE)
  plot(model, pages = 1, seWithMean = TRUE)

  ## GAM check
  gam.check(model)

  ## acf (corrected/normal residuals depending on ar model or not)
  r <- if (!is.null(model$std.rsd)) model$std.rsd else residuals(model)
  acf(r, main = "ACF of residuals")

  ## returning the model itself
  invisible(model)
}

## -----------------------------------------------------------
## simple model with hourly cycles
## -------------------------------

simple_model <- gam(value ~ s(hour, bs = "cc", k = 20),
  method = "REML",
  family = tw(link = "log"),
  knots = list(hour = c(0, 24)),
  data = df_no2
)

diagnostics(simple_model)

## -----------------------------------------------------------
## including specified station type (bg, tr)
## -----------------------------------------

station_type_model <- bam(
  value ~ station_type +
    s(hour, by = station_type, bs = "cc", k = 24),
  method = "fREML",
  family = tw(link = "log"),
  knots = list(hour = c(0, 24)),
  data = df_no2
)

diagnostics(station_type_model)

AIC(simple_model, station_type_model)

## -----------------------------------------------------------
## including monthly cycles
## ------------------------

monthly_model <- bam(
  value ~ station_type +
    s(hour, by = station_type, bs = "cc", k = 24) +
    s(month, bs = "cc", k = 12),
  method = "fREML",
  family = tw(link = "log"),
  knots = list(hour = c(0, 24), month = c(0.5, 12.5)),
  data = df_no2
)

diagnostics(monthly_model)

AIC(simple_model, station_type_model, monthly_model)

## -----------------------------------------------------------
## including interaction between station type and weekend or not
## -------------------------------------------------------------

## modeling interaction
df_no2$type_day <- interaction(df_no2$station_type, df_no2$weekend)

interaction_model <- bam(
  value ~ type_day + s(hour, by = type_day, bs = "cc", k = 24)
  + s(month, bs = "cc", k = 12),
  method = "fREML",
  family = tw(link = "log"),
  knots = list(hour = c(0, 24), month = c(0.5, 12.5)),
  data = df_no2
)

diagnostics(interaction_model)

AIC(simple_model, station_type_model, monthly_model, interaction_model)

## -----------------------------------------------------------
## including time since beginning of measurement
## ---------------------------------------------

time_model <- bam(
  value ~ type_day + s(hour, by = type_day, bs = "cc", k = 24) +
    s(month, bs = "cc", k = 12) + s(time_index, k = 20),
  method = "fREML",
  family = tw(link = "log"),
  knots = list(hour = c(0, 24), month = c(0.5, 12.5)),
  data = df_no2
)

diagnostics(time_model)

AIC(
  simple_model, station_type_model, monthly_model, interaction_model,
  time_model
)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## AR(1) model
## ~~~~~~~~~~~

## -----------------------------------------------------------
## residuals for AR(1) model fit
## -----------------------------

## residuals
r <- residuals(time_model)

## exact value of lag 1 (rho)
rho_e <- acf(r, plot = FALSE)$acf[2] # lag 1
rho_e

## -----------------------------------------------------------
## altering df for analysis with AR(1)
## -----------------------------------

df_no2 <- df_no2 |>
  dplyr::filter(!is.na(value)) |>
  dplyr::arrange(station, date_start) |>
  dplyr::group_by(station) |>
  dplyr::mutate(
    gap = as.numeric(difftime(date_start,
      dplyr::lag(date_start),
      units = "hours"
    )),
    ar_start = is.na(gap) | gap > 1
  ) |>
  dplyr::ungroup()

## -----------------------------------------------------------
## fitting first AR(1) model
## -------------------------

ar1_model <- bam(
  value ~ type_day + s(hour, by = type_day, bs = "cc", k = 24) +
    s(month, bs = "cc", k = 12) + s(time_index, k = 20),
  rho = rho_e,
  AR.start = ar_start,
  discrete = TRUE,
  method = "fREML",
  family = tw(link = "log"),
  knots = list(hour = c(0, 24), month = c(0.5, 12.5)),
  data = df_no2
)

diagnostics(ar1_model)

## -----------------------------------------------------------
## grid search for optimal rho
## ---------------------------

## rho values to search on
rho_grid <- seq(0.70, 0.95, by = 0.05)

## preparing result df
grid_search_tw <- data.frame(
  rho = rho_grid,
  lag1 = NA_real_,
  freml = NA_real_
)

## iterating over rho values and getting lag1 and freml values
for (i in seq_along(rho_grid)) {
  r <- rho_grid[i]

  mdl <- bam(
    value ~ type_day + s(hour, by = type_day, bs = "cc", k = 24) +
      s(month, bs = "cc", k = 12) + s(time_index, k = 20),
    rho = r,
    AR.start = ar_start,
    discrete = TRUE,
    method = "fREML",
    family = tw(link = "log"),
    knots = list(hour = c(0, 24), month = c(0.5, 12.5)),
    data = df_no2
  )

  grid_search_tw$lag1[i] <- acf(mdl$std.rsd, plot = FALSE)$acf[2]
  grid_search_tw$freml[i] <- mdl$gcv.ubre
  cat("rho =", r, "finished\n")
}

grid_search_tw

## -----------------------------------------------------------
## function for grid search visualization
## --------------------------------------

grid_search_vis <- function(data, x = rho, y, title = NULL, y_axis = NULL,
                            vline = NULL) {
  p <- ggplot(data, aes(x = {{ x }}, y = {{ y }})) +
    geom_line() +
    labs(x = expression(rho), y = y_axis, title = title) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(), axis.line = element_line("black"),
      axis.ticks = element_line("black")
    )

  if (!is.null(vline)) {
    p <- p + geom_vline(
      xintercept = vline,
      linetype = "dashed"
    )
  }
  p
}

grid_search_vis(data = grid_search_tw, y = freml, y_axis = "fREML") /
  grid_search_vis(data = grid_search_tw, y = lag1, y_axis = "Lag-1")

## no minimum in the freml graph
## => AR(1) is estimated for tweedie

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## gaussian distribution
## ~~~~~~~~~~~~~~~~~~~~~

## -----------------------------------------------------------------------------
## gaussian family model
## ---------------------

## AR(1) is exact for the gaussian model, approximated for tweedie
## => finding rho more precisely is the goal

## replacing ≤0 values with minimal values in df
min_pos <- min(df_no2$value[df_no2$value > 0])
df_no2$value_pos <- pmax(df_no2$value, min_pos)

## fitting new model w/o AR(1)
gaussian_model <- bam(
  log(value_pos) ~ type_day +
    s(hour, by = type_day, bs = "cc", k = 24) +
    s(month, bs = "cc", k = 12) + s(time_index, k = 20),
  method = "fREML",
  family = gaussian(),
  knots = list(hour = c(0, 24), month = c(0.5, 12.5)),
  data = df_no2
)

diagnostics(gaussian_model)

## -----------------------------------------------------------
## getting new rho for gaussian model
## ----------------------------------

r_gau <- residuals(gaussian_model)

rho_gau <- acf(r_gau, plot = FALSE)$acf[2] # lag 1
rho_gau

## -----------------------------------------------------------
## fit AR(1) model with gaussian family
## ------------------------------------

ar1_model_gaussian <- bam(
  log(value_pos) ~ type_day +
    s(hour, by = type_day, bs = "cc", k = 24) +
    s(month, bs = "cc", k = 12) + s(time_index, k = 20),
  rho = rho_gau,
  AR.start = ar_start,
  discrete = TRUE,
  method = "fREML",
  family = gaussian(),
  knots = list(hour = c(0, 24), month = c(0.5, 12.5)),
  data = df_no2
)

diagnostics(ar1_model_gaussian)

AIC(gaussian_model, ar1_model_gaussian)

## -----------------------------------------------------------
## grid search for optimal rho again
## ---------------------------------

## preparing result df
grid_search_ga <- data.frame(
  rho = rho_grid,
  lag1 = NA_real_,
  freml = NA_real_
)

## iterating over rho values and getting lag1 and freml values
for (i in seq_along(rho_grid)) {
  r <- rho_grid[i]

  mdl <- bam(
    log(value_pos) ~ type_day + s(hour, by = type_day, bs = "cc", k = 24) +
      s(month, bs = "cc", k = 12) + s(time_index, k = 20),
    rho = r,
    AR.start = ar_start,
    discrete = TRUE,
    method = "fREML",
    family = gaussian(),
    knots = list(hour = c(0, 24), month = c(0.5, 12.5)),
    data = df_no2
  )

  grid_search_ga$lag1[i] <- acf(mdl$std.rsd, plot = FALSE)$acf[2]
  grid_search_ga$freml[i] <- mdl$gcv.ubre
  cat("rho =", r, "fertig\n")
}

grid_search_ga

grid_search_vis(data = grid_search_ga, y = freml, y_axis = "fREML") /
  grid_search_vis(data = grid_search_ga, y = lag1, y_axis = "Lag-1")

## grid search results in a minimum at rho = 0.9 for freml
## => performing more granulated search around rho = 0.9

## -----------------------------------------------------------
## finer search in target area
## ---------------------------

## rho values to search on
rho_grid2 <- seq(0.86, 0.94, by = 0.01)

## preparing result df
grid_search_ga2 <- data.frame(
  rho = rho_grid2,
  lag1 = NA_real_,
  freml = NA_real_
)

## iterating over rho values and getting lag1 and freml values
for (i in seq_along(rho_grid2)) {
  r2 <- rho_grid2[i]

  mdl <- bam(
    log(value_pos) ~ type_day + s(hour, by = type_day, bs = "cc", k = 24) +
      s(month, bs = "cc", k = 12) + s(time_index, k = 20),
    rho = r2,
    AR.start = ar_start,
    discrete = TRUE,
    method = "fREML",
    family = gaussian(),
    knots = list(hour = c(0, 24), month = c(0.5, 12.5)),
    data = df_no2
  )

  grid_search_ga2$lag1[i] <- acf(mdl$std.rsd, plot = FALSE)$acf[2]
  grid_search_ga2$freml[i] <- mdl$gcv.ubre
  cat("rho =", r2, "fertig\n")
}

grid_search_ga2

grid_search_vis(
  data = grid_search_ga2, y = freml, y_axis = "fREML",
  vline = 0.9
) /
  grid_search_vis(
    data = grid_search_ga2, y = lag1, y_axis = "Lag-1",
    vline = 0.9
  )

grid_search_vis(
  data = rbind(grid_search_ga, grid_search_ga2),
  y = freml, y_axis = "fREML", vline = 0.9
) /
  grid_search_vis(
    data = rbind(grid_search_ga, grid_search_ga2),
    y = lag1, y_axis = "Lag-1", vline = 0.9
  )

## rho = 0.90 seems to be the best fit (min(fREML))
## => rho = 0.9 is chosen for both models (gaussian, tweedie) for comparison

## -----------------------------------------------------------------------------
## function for model fitting
## --------------------------

model_fit <- function(data, rho, family) {
  if (family == "gaussian") {
    mdl <- bam(
      log(value_pos) ~ type_day +
        s(hour, by = type_day, bs = "cc", k = 24) +
        s(month, bs = "cc", k = 12) + s(time_index, k = 20),
      rho = rho,
      AR.start = data$ar_start,
      discrete = TRUE,
      method = "fREML",
      family = gaussian(),
      knots = list(hour = c(0, 24), month = c(0.5, 12.5)),
      data = data
    )
  } else if (family == "tweedie") {
    mdl <- bam(
      value ~ type_day + s(hour, by = type_day, bs = "cc", k = 24) +
        s(month, bs = "cc", k = 12) + s(time_index, k = 20),
      rho = rho,
      AR.start = data$ar_start,
      discrete = TRUE,
      method = "fREML",
      family = tw(link = "log"),
      knots = list(hour = c(0, 24), month = c(0.5, 12.5)),
      data = data
    )
  } else {
    stop("unknown family")
  }

  diagnostics(mdl)
}

## -----------------------------------------------------------------------------
## fitting final models for comparison with rho = 0.9 from gaussian
## ----------------------------------------------------------------

model_fit(df_no2, rho = 0.9, family = "gaussian")
final_model <- model_fit(df_no2, rho = 0.9, family = "tweedie")
final_model

## important values
## reference value: background / weekday
exp(3.113434) # 22.49817
## CI
cat("[", exp(3.113434 - 1.96 * 0.005087), ", ", 
      exp(3.113434 + 1.96 * 0.005087), "]", sep = "")
## [22.27497, 22.72361]

## traffic / weekday
exp(0.584178) # 1.793516
## CI
cat("[", exp(0.584178 - 1.96 * 0.006847), ", ", 
      exp(0.584178 + 1.96 * 0.006847), "]", sep = "")
## [1.769608, 1.817748]

## background / weekend
exp(-0.234284) # 0.7911371
## CI
cat("[", exp(-0.234284 - 1.96 * 0.007539), ", ", 
      exp(-0.234284 + 1.96 * 0.007539), "]", sep = "")
## [0.7795329, 0.8029141]

## traffic / weekend
exp(0.365318) # 1.440972
## CI
cat("[", exp(0.365318 - 1.96 * 0.008322), ", ", 
      exp(0.365318 + 1.96 * 0.008322), "]", sep = "")
## [1.417659, 1.464669]

## -----------------------------------------------------------------------------
## readme plots
## ------------

## acf plots
png(here::here("figures", "acf_comparison.png"),
    width = 10, height = 4, units = "in", res = 150)

par(mfrow = c(1, 2))
acf(residuals(time_model), main = "w/o AR(1)")
acf(ar1_model$std.rsd,      main = "w/ AR(1), ρ = 0.90")

dev.off()

## tweedie effect plots (actual)
final_tw <- model_fit(df_no2, rho = 0.9, family = "tweedie")

p <- gratia::draw(final_tw, ncol = 3)
ggsave(here::here("figures", "partial_effects_tweedie_ar1.png"),
       p, width = 10, height = 6, dpi = 150)

## tweedie effect plots (multiplicative)
p <- gratia::draw(final_tw, fun = exp, ncol = 3) &
  geom_hline(yintercept = 1, linetype = "dashed")

ggsave(here::here("figures", "partial_effects_multiplicative.png"),
       p, width = 10, height = 6, dpi = 150)

## grid search
grid_all <- rbind(grid_search_ga, grid_search_ga2) |>
  dplyr::distinct(rho, .keep_all = TRUE) |>
  dplyr::arrange(rho)

p_rho <- (grid_search_vis(grid_all, y = freml, y_axis = "fREML", vline = 0.9) & 
            theme(axis.title.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.text.x = element_blank())) /
  grid_search_vis(grid_all, y = lag1,  y_axis = "Lag-1-ACF", vline = 0.9)

ggsave(here::here("figures", "rho_grid_search.png"),
       p_rho, width = 6, height = 6, dpi = 150)

## -----------------------------------------------------------------------------

