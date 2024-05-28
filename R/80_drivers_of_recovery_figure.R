source("helper_functions.R")
cr_check_packages()

refit_models <- FALSE


###### GAM model - recovery to historical benchmark - autocorrelative term included

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

## More processing
{
  dat1 <- manual |>
    filter(!is.na(rel.perc.benchmark.recovery) & years.benchmark.recovery > 0) |>
    mutate(
      dat1Reef_name = factor(Reef_name),
      recovery.year = as.numeric(recovery.year)
    )
}

## Fit gam model
{
  model1.1 <- gam(
    rel.perc.benchmark.recovery ~ s(log(Disturb.cover)) +
      s(Reef_name, bs = "re") +
      s(bench.cover) +
      s(years.benchmark.recovery),
    cor = corAR1(form = ~ recovery.year | Reef_name),
    # s(recovery.year),
    data = dat1,
    REML = T,
    family = "scat"
  )
}

## model diagnostics
{
  resids <- simulateResiduals(model1.1)
  res <- recalculateResiduals(resids, group = dat1$recovery.year)
  testTemporalAutocorrelation(res, time = unique(dat1$recovery.year)) #### no evidence of AC


  plot(model1.1, pages = 1, scale = 0, seWithMean = T, shift = coef(model1.1)[1])
  gam.check(model1.1)
  mgcv::concurvity(model1.1)
}
## Summarise model
{
  summary(model1.1)
}

## pretty plots
##dist.cover
{
  gam.pred.distcover.bench <- emmeans(model1.1,
    specs = "log(Disturb.cover)",
    type = "response",
    rg.limit = 13000,
    at = list(`log(Disturb.cover)` = seq(min(log(manual$Disturb.cover)),
      max(log(manual$Disturb.cover)),
      length = 100
    ))
  ) |>
    as.data.frame() |>
    mutate(Disturb.cover = exp(`log(Disturb.cover)`))

  ## plot

  dist.cover.bench.gamplot <- ggplot(gam.pred.distcover.bench, aes(x = Disturb.cover, y = emmean)) +
    geom_line() +
    geom_line(aes(y = lower.CL), linetype = "dashed") +
    geom_line(aes(y = upper.CL), linetype = "dashed") +
    geom_rug(data = dat1, aes(x = Disturb.cover, y = rel.perc.benchmark.recovery), sides = "b", length = unit(1, "mm")) +
    scale_x_log10("Post-disturbance cover (%)", labels = function(x) x * 100) +
    scale_y_continuous("Recovery to historical benchmark (%)", limits = c(60, 120)) +
    ggtitle("b)") +
    theme_classic()

  dist.cover.bench.gamplot
}

#bench.cover
{
  gam.pred.benchcover.bench <- emmeans(model1.1, ~bench.cover,
    type = "response",
    rg.limit = 13000,
    at = list(bench.cover = seq(min(manual$bench.cover),
      max(manual$bench.cover),
      length = 100
    ))
  ) %>%
    as.data.frame()

  bench.cover.bench.gamplot <- ggplot(gam.pred.benchcover.bench, aes(x = bench.cover, y = emmean)) +
    geom_line() +
    geom_line(aes(y = lower.CL), linetype = "dashed") +
    geom_line(aes(y = upper.CL), linetype = "dashed") +
    geom_rug(data = dat1, aes(x = bench.cover, y = rel.perc.benchmark.recovery), sides = "b", length = unit(1, "mm")) +
    scale_x_continuous("Historic benchmark cover (%)", labels = function(x) x * 100) +
    scale_y_continuous("Recovery to historical benchmark (%)", limits = c(60, 120)) +
    ggtitle("d)") +
    theme_classic()

  bench.cover.bench.gamplot
}

## historical recovery
{
  newdata <- manual |> # filter(!is.na(rel.perc.benchmark.recovery)) |>
    filter(!is.na(rel.perc.benchmark.recovery) & years.benchmark.recovery > 0) |>
    mutate(years.benchmark.recovery = as.numeric(years.benchmark.recovery)) |>
    tidyr::expand(years.benchmark.recovery)

  gam.pred.recoveryyear.bench <- emmeans(model1.1, ~years.benchmark.recovery, type = "response", at = newdata) |>
    as.data.frame()

  recovery.year.bench.gamplot <- ggplot(gam.pred.recoveryyear.bench, aes(x = years.benchmark.recovery, y = emmean)) +
    geom_line() +
    geom_line(aes(y = lower.CL), linetype = "dashed") +
    geom_line(aes(y = upper.CL), linetype = "dashed") +
    geom_rug(data = dat1, aes(x = years.benchmark.recovery, y = rel.perc.benchmark.recovery), sides = "b", length = unit(1, "mm")) +
    scale_x_continuous("Years for recovery", limits = c(0, 40)) + # ,labels=function(x)x*1)+
    scale_y_continuous("Recovery to historical benchmark (%)", limits = c(60, 100)) +
    ggtitle("f)") +
    theme_classic()

  recovery.year.bench.gamplot
}

## GAM model - recovery to prior benchmark - autocorrelative term included
{
  dat2 <- manual |>
    filter(!is.na(rel.perc.prior.recovery) & years.prior.recovery > 0) |>
    mutate(
      dat1Reef_name = factor(Reef_name),
      recovery.year = as.numeric(recovery.year)
    )

  model2.1 <- gam(
    rel.perc.prior.recovery ~ s(log(Disturb.cover)) +
      s(Reef_name, bs = "re") +
      s(prior.cover) +
      s(years.prior.recovery),
    cor = corAR1(form = ~ recovery.year | Reef_name),
    data = dat1, REML = T, family = "scat"
  )
}
## model diagnostics
{
  plot(model2.1, pages = 1, scale = 0, seWithMean = T, shift = coef(model2.1)[1])
  gam.check(model2.1)
  mgcv::concurvity(model2.1)

  resids <- simulateResiduals(model2.1)
  res <- recalculateResiduals(resids, group = dat2$recovery.year)
  testTemporalAutocorrelation(res, time = unique(dat2$recovery.year))
}
## Summarise model
{
  summary(model2.1)
}
## Plot
{
  gam.pred.distcover.prior <- emmeans(model2.1,
    specs = "log(Disturb.cover)",
    type = "response",
    rg.limit = 13000,
    at = list(`log(Disturb.cover)` = seq(min(log(manual$Disturb.cover)),
      max(log(manual$Disturb.cover)),
      length = 100
    ))
  ) %>%
    as.data.frame() %>%
    mutate(Disturb.cover = exp(`log(Disturb.cover)`))

  dist.cover.prior.gamplot <- ggplot(gam.pred.distcover.prior, aes(x = Disturb.cover, y = emmean)) +
    geom_line() +
    geom_line(aes(y = lower.CL), linetype = "dashed") +
    geom_line(aes(y = upper.CL), linetype = "dashed") +
    geom_rug(data = dat1, aes(x = Disturb.cover, y = rel.perc.benchmark.recovery), sides = "b", length = unit(1, "mm")) +
    scale_x_log10("Post-disturbance cover (%)", labels = function(x) x * 100) +
    scale_y_continuous("Recovery to prior benchmark (%)", limits = c(50, 140)) +
    ggtitle("a)") +
    theme_classic()

  dist.cover.prior.gamplot
}
{
  
  gam.pred.benchcover.prior <- emmeans(model2.1, ~prior.cover,
    type = "response",
    rg.limit = 13000,
    at = list(prior.cover = seq(min(manual$prior.cover), max(manual$prior.cover),
      length = 100
    ))
  ) |>
    as.data.frame()

  prior.cover.bench.gamplot <- ggplot(gam.pred.benchcover.prior, aes(x = prior.cover, y = emmean)) +
    geom_line() +
    geom_line(aes(y = lower.CL), linetype = "dashed") +
    geom_line(aes(y = upper.CL), linetype = "dashed") +
    geom_rug(data = dat1, aes(x = prior.cover, y = rel.perc.benchmark.recovery), sides = "b", length = unit(1, "mm")) +
    scale_x_continuous("Prior benchmark cover (%)", labels = function(x) x * 100) +
    scale_y_continuous("Recovery to prior benchmark (%)", limits = c(50, 140)) +
    ggtitle("c)") +
    theme_classic()

  prior.cover.bench.gamplot

}
{
  newdata1 <- manual |>
    filter(!is.na(rel.perc.prior.recovery)) |>
    mutate(years.prior.recovery = as.numeric(years.prior.recovery)) |>
    tidyr::expand(years.prior.recovery)

  gam.pred.recoveryyear.prior <- emmeans(model2.1, ~years.prior.recovery, type = "response", at = newdata1) |>
    as.data.frame()

  recovery.year.prior.gamplot <- ggplot(gam.pred.recoveryyear.prior, aes(x = years.prior.recovery, y = emmean)) +
    geom_line() +
    geom_line(aes(y = lower.CL), linetype = "dashed") +
    geom_line(aes(y = upper.CL), linetype = "dashed") +
    geom_rug(data = dat1, aes(x = years.prior.recovery, y = rel.perc.benchmark.recovery), sides = "b", length = unit(1, "mm")) +
    scale_x_continuous("Years for recovery", limits = c(0, 40)) + # ,labels=function(x)x*100)+
    scale_y_continuous("Recovery to prior benchmark (%)", limits = c(60, 140)) +
    ggtitle("e)") +
    theme_classic()

  recovery.year.prior.gamplot
}
## Put them all together
{
  gam.plots <-
    (dist.cover.prior.gamplot / prior.cover.bench.gamplot / recovery.year.prior.gamplot) |
    (dist.cover.bench.gamplot / bench.cover.bench.gamplot / recovery.year.bench.gamplot)
  gam.plots <- gam.plots & theme_classic(base_size = 9)

   ggsave(
     filename = "../outputs/figures/figure_8.png",
     gam.plots,
     height = 20, width = 18, units = "cm", dpi = 600
   )
   ggsave(
     filename = "../outputs/figures/figure_8.pdf",
     gam.plots,
     height = 20, width = 18, units = "cm"
   )
}


