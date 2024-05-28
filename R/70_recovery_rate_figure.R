source("helper_functions.R")
cr_check_packages()

refit_models <- FALSE

### Production of Figure 7 changes in the rate of recovery through time

## Load data
{
  manual <- read.csv(
    file = "../data/processed/manual recovery indiv manta reefs_2022.csv",
    strip.white = TRUE
  )
}

## Data processing
{
  ## data wrangling
  ## calculate decline and recovery of hard coral cover to historical and prior benchmark
  ## calculate rate of hard coral recovery to historical and prior benchmark
  ## assign categorical variable of recovery (yes or no)
  manual <- manual |>
     mutate(
       abs.decline.bench = ifelse(Disturb.year <= Bench.year, NA,
         bench.cover - Disturb.cover),
       rel.decline.bench = ifelse(Disturb.year <= Bench.year, NA,
         abs.decline.bench / bench.cover * 100),
       abs.decline.prior = prior.cover - Disturb.cover,
       rel.decline.prior = abs.decline.prior / prior.cover * 100,
       rel.perc.benchmark.recovery = ifelse(recovery.year == Bench.year, NA,
         recovery.cover / bench.cover * 100),
       rel.perc.prior.recovery = ifelse(recovery.year == prior.year, NA,
         recovery.cover / prior.cover * 100),
       years.benchmark.recovery = (recovery.year - Bench.year),
       years.prior.recovery = recovery.year - prior.year,
       recovery.interval = recovery.year - Disturb.year,
       recovery.rate = (recovery.cover - Disturb.cover) /
     (recovery.year - Disturb.year),
     diff.2.baseline = recovery.cover - bench.cover
     ) |>
     mutate(
       Region = factor(Region, levels = c("Northern", "Central", "Southern")),
       Reef_name = factor(Reef_name),
       Disturbance = factor(Disturbance,
         levels = c("Bleaching", "COTS", "Cyclone", "Multiple", "Unknown")),
       Disturbance2 = factor(Disturbance2,
         levels = c("Bleaching", "COTS", "Cumulative", "Cyclone", "Unknown")),
       Dist1 = factor(Dist1,
         levels = c("Bleaching", "COTS", "Cyclone", "Unknown", "Multiple")
       )
     )

  manual <- manual |>
     dplyr::mutate(
       recovered.prior = ifelse(rel.perc.prior.recovery >= 100, "yes", "no"),
       recovered.bench = ifelse(rel.perc.benchmark.recovery >= 100, "yes", "no")
     )
}


## recovery vs year of recovery
{
  recov.rate.summary <-
    manual |>
    group_by(recovery.year) |>
    summarise(mean = mean(recovery.rate))
  
}
## Preliminary Plot
{
  recovery.rate.plot <-
    manual |>
    ggplot(aes(x = recovery.year, y = recovery.rate)) +
    geom_point(stat = "summary", alpha = 0.6,
      fun = "mean", size = 2.5) +
    geom_errorbar(stat = "summary", fun.data = mean_se,
      width = 0, alpha = 0.6) +
    theme_classic()

  recovery.rate.plot
}


## Fit model
{
  if (refit_models) {
    manual |>
      filter(recovery.year != "NA") |>
      droplevels() |> 
      summarise(
        Mean = mean(log(recovery.rate)),
        Median = median(log(recovery.rate)),
        SD = sd(log(recovery.rate)),
        SD_slope = sd(log(recovery.rate)) / sd(recovery.year)
      )
    priors <- prior("normal(-1, 1)", class = "Intercept") +
      prior("normal(0,5)", class = "b") +
      prior(gamma(0.01, 0.01), class = "shape") +
      prior(student_t(3, 0, 1), class = "sd")

    form <- bf(recovery.rate ~ poly(recovery.year - 1980, 3) * prior.cover + (1 | Disturbance) + (1 | Reef_name),
      family = Gamma(link = "log")
    )

    recovery.brm<-brm(form, 
      data=manual %>% filter(recovery.year!='NA'),
      prior=priors,
      iter = 2000, warmup = 500, thin = 5,
      chains = 3, cores = 3,
      control = list(adapt_delta = 0.95, max_treedepth = 20),
      backend = "cmdstanr",
      save_pars=save_pars(all = TRUE))

    saveRDS(recovery.brm, file = "../data/modelled/recovery.brm.RData")
  }
  recovery.brm <- readRDS(file = "../data/modelled/recovery.brm.RData")
}
## MCMC diagnostics
{
  stan_trace(recovery.brm$fit)
  stan_ac(recovery.brm$fit)
  stan_ess(recovery.brm$fit)
  stan_rhat(recovery.brm$fit)
}
## Model validation
{
  recovery.brm |>
    pp_check(type = "dens_overlay", ndraws = 100) +
    scale_x_log10() +
    recovery.brm |>
    pp_check(type = "loo_pit_overlay", ndraws = 100)

  resids <- make_brms_dharma_res(
    recovery.brm,
    integerResponse = FALSE
  )
  wrap_elements(~ testUniformity(resids)) +
    wrap_elements(~ plotResiduals(resids,
      form = factor(rep(1, nrow(recovery.brm$data))))) +
    wrap_elements(~ plotResiduals(resids)) +
    wrap_elements(~ testDispersion(resids))
}
## Model Summary
{
  summary(recovery.brm)
  ## conditional_effects(recovery.brm)
}
## Posteriors
{
  predgrid <-
    with(
      manual |>
        filter(recovery.year != "NA"),
      seq(min(recovery.year), max(recovery.year), length = 1000)
    )

  recovery.fitted.trend <-
    recovery.brm |>
    emmeans(~recovery.year,
      at = list(recovery.year = predgrid),
      type = "response"
    ) |>
    as.data.frame()
}
## Plot
{
  recovery.rate.thru.time.plot <-
    recovery.rate.plot +
    geom_ribbon(data = recovery.fitted.trend,
      aes(y = response, x = recovery.year,
        ymin = lower.HPD, ymax = upper.HPD),
      alpha = 0.2) +
    geom_line(data = recovery.fitted.trend,
      aes(y = response, x = recovery.year), colour = "blue") +
    scale_x_continuous("Year of recovery") +
    scale_y_continuous("Recovery rate (% per year)") +
    theme_classic(base_size = 10) +
    theme(
      axis.title = element_text(size = rel(1.2), face = "bold"),
      axis.text = element_text(size = rel(1.2))
    )

  recovery.rate.thru.time.plot

  ggsave(
    filename = "../outputs/figures/figure_7.pdf",
    recovery.rate.thru.time.plot,
    height = 20/1.6, width = 20, units = "cm"
  )

  ggsave(
    filename = "../outputs/figures/figure_7.png",
    recovery.rate.thru.time.plot,
    height = 20/1.6, width = 20, units = "cm", dpi = 600
  )
}
