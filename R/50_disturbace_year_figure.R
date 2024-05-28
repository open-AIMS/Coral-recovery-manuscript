source("helper_functions.R")
cr_check_packages()

refit_models <- FALSE

## Year of disturbance component
{
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
  ## Fit model
  {
    if (refit_models) {
      manual |> summarise(
        Mean = mean(log(Disturb.cover)),
        Median = median(log(Disturb.cover)),
        SD = sd(log(Disturb.cover)),
        SD_slope = sd(log(Disturb.cover)) / sd(Disturb.year)
      )
      priors <- prior("normal(-2, 1)", class = "Intercept") +
        prior("normal(0,2)", class = "b") +
        prior(student_t(3, 0, 1), class = "sd") +
        prior(gamma(0.01, 0.01), class = "shape")

      form <- bf(Disturb.cover ~ scale(Disturb.year) + (1 | Disturbance) + (1 | Reef_name),
        family = Gamma(link = "log")
      )
      post.dist.cover.brm <- brm(
        form = form,
        prior = priors,
        data = manual,
        iter = 2000, warmup = 500, thin = 5,
        chains = 3, cores = 3,
        control = list(adapt_delta = 0.95, max_treedepth = 20),
        backend = "cmdstanr"
      )
      saveRDS(post.dist.cover.brm, file = "../data/modelled/post.dist.cover.brm.RData")
    }
    post.dist.cover.brm <- readRDS(file = "../data/modelled/post.dist.cover.brm.RData")
  }
  ## MCMC diagnostics
  {
    stan_trace(post.dist.cover.brm$fit)
    stan_ac(post.dist.cover.brm$fit)
    stan_ess(post.dist.cover.brm$fit)
    stan_rhat(post.dist.cover.brm$fit)
  }
  ## Model validation
  {
    post.dist.cover.brm |>
      pp_check(type = "dens_overlay", ndraws = 100) +
      scale_x_log10() +
      post.dist.cover.brm |>
      pp_check(type = "loo_pit_overlay", ndraws = 100)

    resids <- make_brms_dharma_res(
      post.dist.cover.brm,
      integerResponse = FALSE
    )
    wrap_elements(~ testUniformity(resids)) +
      wrap_elements(~ plotResiduals(resids,
        form = factor(rep(1, nrow(post.dist.cover.brm$data))))) +
      wrap_elements(~ plotResiduals(resids)) +
      wrap_elements(~ testDispersion(resids))
  }
  ## Model Summary
  {
    summary(post.dist.cover.brm)
    conditional_effects(post.dist.cover.brm)
  }
  ## Alternative summary
  {
    post.dist.cover.brm |>
      as_draws_df() |>
      mutate(across(starts_with("b_"), exp)) |>
      dplyr::select(-starts_with("r_"),
        -starts_with("lp"),
        -starts_with("Intercept")) |>
      summarise_draws(median, HDInterval::hdi )
  }
  ## Calculating percent change for one unit change in Disturbance year
  {
    post.dist.cover.brm |>
      emmeans(~Disturb.year, at = list(Disturb.year = c(2020, 2021)), type = "response") |>
      gather_emmeans_draws() |>
      group_by(.draw) |>
      summarise(.value = 100 * (exp(diff(.value)) - 1)) |>
      summarise_draws(median, HDInterval::hdi, P = ~ mean(.x < 0))
  }
  ## old way - note, this will express change per 1 unit of scaled
  ## disturbance year rather than 1 unit of Disturbance year
  {
    ###extracting for linear
    mcmc=post.dist.cover.brm %>% as.matrix()
    post.dist.cover.slope<-median_hdci(exp(mcmc[,2]))

    ### calculating 'slope' (percentage change +/- CIs)
    (1-post.dist.cover.slope[1])*100
    (1-post.dist.cover.slope[2])*100
    (1-post.dist.cover.slope[3])*100
    ### exceedence probabilty
    sum(mcmc[,2]<0)/length(mcmc[,2])
  }
  ## Plot trend
  {
    predgrid <- with(
      manual,
      seq(min(Disturb.year), max(Disturb.year), length = 1000)
    )

    post.cover.fitted.trend <-
      post.dist.cover.brm |>
      emmeans(~Disturb.year,
        at = list(Disturb.year = predgrid),
        type = "response"
      ) |>
      as.data.frame()

    plot3 <-
      manual |>
      ggplot(aes(x = Disturb.year, y = Disturb.cover)) +
      geom_point(shape = 16, alpha = 0.6, stat = "summary", fun.y = "mean") +
      geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.1) +
      scale_y_continuous("Post disturbance coral cover (%)",
        breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
        labels = c(0, 10, 20, 30, 40, 50)
      ) +
      scale_x_continuous("Year of disturbance") +
      scale_size_area(name = "Pre-disturbance\n cover") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 0))

    plot3 <-
      plot3 +
      geom_ribbon(data = post.cover.fitted.trend,
        aes(y = response, x = Disturb.year,
          ymin = lower.HPD, ymax = upper.HPD), alpha = 0.2) +
      geom_line(data = post.cover.fitted.trend,
        aes(y = response, x = Disturb.year), colour = "blue") +
      ggtitle("a)")

    plot3
  }
}

##### panel b) has the relative loss of coral increased through time
{
  ## Fit model
  {
    if (refit_models) {
      manual |> summarise(
        Mean = mean(log(rel.decline.prior)),
        Median = median(log(rel.decline.prior)),
        SD = sd(log(rel.decline.prior)),
        SD_slope = sd(log(rel.decline.prior)) / sd(Disturb.year)
      )
      priors <- prior("normal(4, 1)", class = "Intercept") +
        prior("normal(0,2)", class = "b") +
        prior(student_t(3, 0, 1), class = "sd") +
        prior(gamma(0.01, 0.01), class = "shape")

      form <- bf(rel.decline.prior ~ scale(Disturb.year) + (1 | Disturbance) + (1 | Reef_name),
        family = Gamma(link = "log")
      )

      rel.coral.loss.brm <- brm(
        form = form,
        prior = priors,
        data = manual,
        family = Gamma(link = "log"),
        chains = 3, cores = 3,
        iter = 2000, warmup = 500, thin = 5,
        control = list(adapt_delta = 0.95, max_treedepth = 20),
        backend = "cmdstanr",
        save_pars = save_pars(all = TRUE)
      )
      saveRDS(rel.coral.loss.brm, file = "../data/modelled/rel.coral.loss.brm.RData")
    }
    rel.coral.loss.brm <- readRDS(file = "../data/modelled/rel.coral.loss.brm.RData")
  }
  ## MCMC diagnostics
  {
    stan_trace(rel.coral.loss.brm$fit)
    stan_ac(rel.coral.loss.brm$fit)
    stan_ess(rel.coral.loss.brm$fit)
    stan_rhat(rel.coral.loss.brm$fit)
  }
  ## Model validation
  {
    rel.coral.loss.brm |>
      pp_check(type = "dens_overlay", ndraws = 100) +
      scale_x_log10() +
      rel.coral.loss.brm |>
      pp_check(type = "loo_pit_overlay", ndraws = 100)

    resids <- make_brms_dharma_res(
      rel.coral.loss.brm,
      integerResponse = FALSE
    )
    wrap_elements(~ testUniformity(resids)) +
      wrap_elements(~ plotResiduals(resids,
        form = factor(rep(1, nrow(rel.coral.loss.brm$data))))) +
      wrap_elements(~ plotResiduals(resids)) +
      wrap_elements(~ testDispersion(resids))
  }
  ## Model Summary
  {
    summary(rel.coral.loss.brm)
    conditional_effects(rel.coral.loss.brm)
  }
  ## Alternative summary
  {
    rel.coral.loss.brm |>
      as_draws_df() |>
      mutate(across(starts_with("b_"), exp)) |>
      dplyr::select(-starts_with("r_"),
        -starts_with("z_"),
        -starts_with("lp"),
        -starts_with("Intercept")) |>
      summarise_draws(median, HDInterval::hdi )
  }
  ## Calculating percent change for one unit change in Disturbance year
  {
    rel.coral.loss.brm |>
      emmeans(~Disturb.year, at = list(Disturb.year = c(2020, 2021)), type = "response") |>
      gather_emmeans_draws() |>
      group_by(.draw) |>
      summarise(.value = 100 * (exp(diff(.value)) - 1)) |>
      summarise_draws(median, HDInterval::hdi, P = ~ mean(.x > 0))
  }
  ## old way - note, this will express change per 1 unit of scaled
  ## disturbance year rather than 1 unit of Disturbance year
  {
    ###extracting for linear
    mcmc=rel.coral.loss.brm %>% as.matrix()
    rel.coral.loss.brm.slope<-median_hdci(exp(mcmc[,2]))

    ### calculating 'slope' (percentage change +/- CIs)
    (1-rel.coral.loss.brm.slope[1])*100
    (1-rel.coral.loss.brm.slope[2])*100
    (1-rel.coral.loss.brm.slope[3])*100
    ### exceedence probabilty
    sum(mcmc[,2]<0)/length(mcmc[,2])
  }
  ## Plot trend
  {
    predgrid <- with(
      manual,
      seq(min(Disturb.year), max(Disturb.year), length = 1000)
    )

    rel.coral.loss.fitted.trend <-
      rel.coral.loss.brm |>
      emmeans(~Disturb.year,
        at = list(Disturb.year = predgrid),
        type = "response"
      ) |>
      as.data.frame()

    rel.coral.loss.thru.time.plot <-
      manual |> ggplot(aes(x = Disturb.year, y = rel.decline.prior)) +
      geom_point(shape = 16, alpha = 0.6,
        stat = "summary", fun.y = "mean") +
      geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.1) +
      scale_x_continuous("Year of disturbance") +
      scale_y_continuous("Relative coral loss (%)") +
      ggtitle("b)") +
      theme_classic()

    rel.coral.loss.thru.time.plot <-
      rel.coral.loss.thru.time.plot +
      geom_ribbon(data = rel.coral.loss.fitted.trend,
        aes(y = response, x = Disturb.year,
          ymin = lower.HPD, ymax = upper.HPD), alpha = 0.2) +
      geom_line(data = rel.coral.loss.fitted.trend,
        aes(y = response, x = Disturb.year), colour = "blue") +
      ggtitle("b)")

    rel.coral.loss.thru.time.plot
  }
}

### panel c) disturbance number and interval per decade ------------------------------
{
  ## Preparations
  {
    # number of reefs with disturbances - redoing as average number of reefs per year as per reviewer 2
    manual$Dist.Decade <- factor(manual$Dist.Decade,
      levels = c("86 to 90", "91 to 00", "01 to 10", "11 to 20")
    )

    reefs.with.disturbance.dec <- # length(unique(manual$Reef_name))
      manual |>
      group_by(Dist.Decade, Disturb.year) |>
      dplyr::summarise(length(unique(Reef_name)))

    colnames(reefs.with.disturbance.dec)[3] <- c("Reefs.with.disturbance")
    colnames(reefs.with.disturbance.dec)[2] <- c("year")
  }
  ## Get and process data
  {
    ### number of reefs surveyed
    manta <- read.csv(
      file = "../data/processed/manta tow by reef 2021.csv",
      strip.white = T
    )
    head(manta)

    manta <- manta |>
      filter(REPORT_YEAR > 1985) |>
      filter(SECTOR != "TS") |>
      mutate(cREPORT_YEAR = factor(REPORT_YEAR)) |>
      mutate(SECTOR = factor(SECTOR)) |>
      mutate(Dist.Decade = factor(Dist.Decade,
        levels = c("86 to 90", "91 to 00", "01 to 10", "11 to 20"))) |>
      mutate(Region = factor(Region))
  }
  ## More preparations
  {
    reefs.surveyed.dec <-
      manta |>
      group_by(Dist.Decade, REPORT_YEAR) |>
      dplyr::summarise(length(unique(REEF_NAME))) |>
      as.data.frame() #|> rename(reefs.surveyed=length(unique(REEF_NAME)))
    colnames(reefs.surveyed.dec)[2] <- c("year")
    colnames(reefs.surveyed.dec)[3] <- c("Reefs.surveyed")

    reefs.with.disturbance.dec <-
      reefs.with.disturbance.dec |>
      left_join(reefs.surveyed.dec) |>
      as.data.frame()

    reefs.with.disturbance.dec <-
      reefs.with.disturbance.dec |>
      mutate(Dist.Decade = factor(Dist.Decade,
        levels = c("86 to 90", "91 to 00", "01 to 10", "11 to 20"))) |>
      mutate(perc.reefs.disturb = Reefs.with.disturbance / Reefs.surveyed * 100)

    reefs.with.disturbance.dec.ave <-
      reefs.with.disturbance.dec |>
      group_by(Dist.Decade) |>
      summarise(
        average = mean(perc.reefs.disturb),
        sd = sd(perc.reefs.disturb),
        se = sd / sqrt(length(perc.reefs.disturb))
      )
  }
  ## Figure
  {
    perc.reefs.disturb <-
      reefs.with.disturbance.dec.ave |>
      filter(Dist.Decade != "NA") |> 
      ggplot(aes(x = Dist.Decade, y = average)) +
        geom_bar(label = "n", stat = "identity", colour = "black") +
      geom_errorbar(stat = "identity",
        aes(ymin = average - se, ymax = average + se), width = 0) +
        geom_text(x = 1, y = -0.5, label = "(348)", size = 2) +
        geom_text(x = 2, y = -0.5, label = "(309)", size = 2) +
        geom_text(x = 3, y = -0.5, label = "(258)", size = 2) +
        geom_text(x = 4, y = -0.5, label = "(204)", size = 2) +
      scale_y_continuous("Percent survey reefs with disturbances",
        limits = c(0, 20)) +
        scale_x_discrete("Decade of survey") +
        ggtitle("c)") +
        theme_classic() +
        theme(
          legend.position = c(0.275, 0.95),
          legend.text = element_text(size = 7.5),
          legend.key.size = unit(0.5, "cm")
        )

    perc.reefs.disturb
  }
}

### panel d) interval between disturbances
{
  ## Bleaching
  {
    ## Get data
    {
      interval <- read.csv(
        file = "../data/processed/220509 bleaching intervals.csv",
        strip.white = TRUE
      )
      head(interval)
      interval <-
        interval |>
        mutate(
          Dist_decade = factor(Dist_decade,
            levels = c("91 to 00", "01 to 10", "11 to 20")),
          Reef = factor(Reef)
        )
    }
    ## Fit model
    {
      if (refit_models) {
        interval |> group_by(Dist_decade) |> summarise(
          Mean = mean(log(Bleach_int)),
          Median = median(log(Bleach_int)),
          SD = sd(log(Bleach_int))
        )
        priors <- prior("normal(2, 2)", class = "Intercept") +
          prior("normal(0,2)", class = "b") +
          prior(student_t(3, 0, 2), class = "sd") +
          prior(gamma(0.01, 0.01), class = "shape")

        form <- bf(Bleach_int~Dist_decade+(1|Reef),
          family=negbinomial(link = "log", link_shape = "log")
        )
        
        bleach.interval.brm <- brm(
          form = form,
          data = interval,
          chains = 3, cores = 3,
          iter = 2000, warmup = 500, thin = 5,
          control = list(adapt_delta = 0.95, max_treedepth = 20),
          backend = "cmdstanr",
          seed = 3
        )
        saveRDS(bleach.interval.brm, file = "../data/modelled/bleach.interval.brm.RData")
      }
      bleach.interval.brm <- readRDS(file = "../data/modelled/bleach.interval.brm.RData")
    }
    ## MCMC diagnostics
    {
      stan_trace(bleach.interval.brm$fit)
      stan_ac(bleach.interval.brm$fit)
      stan_ess(bleach.interval.brm$fit)
      stan_rhat(bleach.interval.brm$fit)
    }
    ## Model validation
    {
      bleach.interval.brm |>
        pp_check(type = "dens_overlay", ndraws = 100) +
        scale_x_log10() +
        bleach.interval.brm |>
        pp_check(type = "loo_pit_overlay", ndraws = 100)

      resids <- make_brms_dharma_res(
        bleach.interval.brm,
        integerResponse = FALSE
      )
      wrap_elements(~ testUniformity(resids)) +
        wrap_elements(~ plotResiduals(resids,
          form = factor(rep(1, nrow(bleach.interval.brm$data))))) +
        wrap_elements(~ plotResiduals(resids)) +
        wrap_elements(~ testDispersion(resids))
    }
    ## Model Summary
    {
      summary(bleach.interval.brm)
      conditional_effects(bleach.interval.brm)
    }
    ## Marginal means
    {
      bleach.interval.est <-
        bleach.interval.brm |> emmeans(~Dist_decade, type = "response") |>
        as.data.frame() |>
        rename(Decade = Dist_decade)
      ## 80s data is missing so add dummy variables
      Decade <- "85 to 90"
      prob <- 0
      lower.HPD <- 0
      upper.HPD <- 0

      eighties <- cbind(Decade, prob, lower.HPD, upper.HPD) |>
        as.data.frame() |>
        mutate(
          Decade = factor(Decade),
          prob = as.numeric(prob),
          lower.HPD = as.numeric(lower.HPD),
          upper.HPD = as.numeric(upper.HPD)
        )

      bleach.interval.est <- rbind(eighties, bleach.interval.est)

      ### exceedence probabilty
      mcmc <- bleach.interval.brm |> as.matrix()
      sum(mcmc[, 2] < 0) / length(mcmc[, 2])

      ## Alternative calculations for boxplots
      bleach.interval.box <-
        bleach.interval.brm |>
        emmeans(~Dist_decade, type = "response") |>
        gather_emmeans_draws() |> 
        mutate(.value = exp(.value)) |> 
        rename(Decade = Dist_decade) |>
        dplyr::select(-.draw, -.chain, -.iteration) |> 
        summarise_draws(median,
          q = ~ boxplot.stats(.x)$stats) 
    }
    ## Plot
    {
      modelled.bleach.interval.plot <-
        bleach.interval.est |> ggplot(aes(x = Decade, y = prob)) +
        geom_bar(stat = "identity", fill = "firebrick3") +
        geom_errorbar(
          stat = "identity",
          aes(ymin = lower.HPD, ymax = upper.HPD), width = 0
        ) +
        ggtitle("d)") +
        scale_y_continuous("Years between disturbance") +
        scale_x_discrete("Decade of first bleaching",
          limits = c("85 to 90", "91 to 00", "01 to 10", "11 to 20"),
          labels = c("85 to 90*", "91 to 00**", "01 to 10", "11 to 20")
        ) +
        theme_classic()
      modelled.bleach.interval.plot  
    }
  }
  ## Cyclones
  {
    ## Get data
    {
      load(file = "../data/processed/cyclone.interval.RData")
    }
    ## Fit model
    {
      if (refit_models) {
        cyclones.interval |> group_by(Decade) |> summarise(
          Mean = mean(log(Interval)),
          Median = median(log(Interval)),
          SD = sd(log(Interval))
        )
        priors <- prior("normal(2, 2)", class = "Intercept") +
          prior("normal(0,2)", class = "b") +
          prior(student_t(3, 0, 2), class = "sd") +
          prior(gamma(0.01, 0.01), class = "shape")

        form <- bf(Interval ~ 0 + Decade + (1 | REEF_NAME),
          family = negbinomial(link = "log", link_shape = "log")
        )
        
        cyclones.interval.brm <- brm(
          form = form,
          data = cyclones.interval,
          chains = 3, cores = 3,
          iter = 2000, warmup = 500, thin = 5,
          control = list(adapt_delta = 0.95, max_treedepth = 20),
          backend = "cmdstanr"
        )
        saveRDS(cyclones.interval.brm, file = "../data/modelled/cyclones.interval.brm.RData")
      }
      cyclones.interval.brm <- readRDS(file = "../data/modelled/cyclones.interval.brm.RData")
    }
    ## MCMC diagnostics
    {
      stan_trace(cyclones.interval.brm$fit)
      stan_ac(cyclones.interval.brm$fit)
      stan_ess(cyclones.interval.brm$fit)
      stan_rhat(cyclones.interval.brm$fit)
    }
    ## Model validation
    {
      cyclones.interval.brm |>
        pp_check(type = "dens_overlay", ndraws = 100) +
        scale_x_log10() +
        cyclones.interval.brm |>
        pp_check(type = "loo_pit_overlay", ndraws = 100)

      resids <- make_brms_dharma_res(
        cyclones.interval.brm,
        integerResponse = FALSE
      )
      wrap_elements(~ testUniformity(resids)) +
        wrap_elements(~ plotResiduals(resids,
          form = factor(rep(1, nrow(cyclones.interval.brm$data))))) +
        wrap_elements(~ plotResiduals(resids)) +
        wrap_elements(~ testDispersion(resids))
    }
    ## Model Summary
    {
      summary(cyclones.interval.brm)
      conditional_effects(cyclones.interval.brm)
    }
    ## Plot
    {
      cyclones.summary <-
        cyclones.interval.brm |>
        emmeans(~Decade, type = "response") |>
        as.data.frame()

      cyclones.summary <-
        cyclones.summary |>
        mutate(Decade = factor(Decade,
          levels = sort(as.numeric(as.character(unique(Decade)))),
          labels = c("85 to 90", "91 to 00", "01 to 10", "11 to 20")
        )) |>
        arrange(Decade)

      modelled.cyclone.interval.plot <-
        cyclones.summary |> ggplot(aes(x = Decade, y = prob)) +
        geom_bar(stat = "identity", fill = "blue3") +
        geom_errorbar(stat = "identity",
          aes(ymin = lower.HPD, ymax = upper.HPD), width = 0) +
        ggtitle("d)") +
        scale_y_continuous("Years between disturbance",
          breaks = c(0, 2, 4, 6, 8, 10)) +
        scale_x_discrete("Decade of first cyclone",
          breaks = c("85 to 90", "91 to 00", "01 to 10", "11 to 20"),
          labels = c("85 to 90", "91 to 00", "01 to 10", "11 to 20")
        ) +
        theme_classic()
      modelled.cyclone.interval.plot

      cyclones.interval.box <-
        cyclones.interval.brm |>
        emmeans(~Decade, type = "response") |>
        gather_emmeans_draws() |> 
        mutate(.value = exp(.value)) |> 
        dplyr::select(-.draw, -.chain, -.iteration) |> 
        summarise_draws(median,
          q = ~ boxplot.stats(.x)$stats) 
    }
  }
  ## COTS
  {
    ## Get data
    {
      cots.interval.outbreak <-
        read.csv(
          file = "../data/processed/cots.interval.outbreak.csv",
          strip.white = TRUE
        )
    }
    ## Fit model
    {
      if (refit_models) {
        cots.interval.outbreak <- cots.interval.outbreak %>%
          mutate(
            REEF_NAME = factor(REEF_NAME),
            Zone = factor(Zone),
            Decade = factor(Decade)
          )

        cots.interval.outbreak |>
          group_by(Decade) |>
          summarise(
            Mean = mean(log(Interval)),
            Median = median(log(Interval)),
            SD = sd(log(Interval))
          )
        priors <- prior("normal(3, 1)", class = "Intercept") +
          prior("normal(0,2)", class = "b") +
          prior(student_t(3, 0, 1), class = "sd") +
          prior(gamma(0.01, 0.01), class = "shape")

        form <- bf(Interval ~ 0 + Decade + (Decade | REEF_NAME),
          family = negbinomial(link = "log", link_shape = "log")
        )

        cots.interval.outbreak.brm <- brm(
          form = form,
          data = cots.interval.outbreak,
          chains = 3, cores = 3,
          iter = 2000, warmup = 500, thin = 5,
          control = list(adapt_delta = 0.95, max_treedepth = 20),
          backend = "cmdstanr"
        )
        saveRDS(cots.interval.outbreak.brm, file = "../data/modelled/cots.interval.outbreak.brm.RData")
      }
      cots.interval.outbreak.brm <- readRDS(file = "../data/modelled/cots.interval.outbreak.brm.RData")
    }
    ## MCMC diagnostics
    {
      stan_trace(cots.interval.outbreak.brm$fit)
      stan_ac(cots.interval.outbreak.brm$fit)
      stan_ess(cots.interval.outbreak.brm$fit)
      stan_rhat(cots.interval.outbreak.brm$fit)
    }
    ## Model validation
    {
      cots.interval.outbreak.brm |>
        pp_check(type = "dens_overlay", ndraws = 100) +
        scale_x_log10() +
        cots.interval.outbreak.brm |>
        pp_check(type = "loo_pit_overlay", ndraws = 100)

      resids <- make_brms_dharma_res(
        cots.interval.outbreak.brm,
        integerResponse = FALSE
      )
      wrap_elements(~ testUniformity(resids)) +
        wrap_elements(~ plotResiduals(resids,
          form = factor(rep(1, nrow(cots.interval.outbreak.brm$data))))) +
        wrap_elements(~ plotResiduals(resids)) +
        wrap_elements(~ testDispersion(resids))
    }
    ## Model Summary
    {
      summary(cots.interval.outbreak.brm)
      conditional_effects(cots.interval.outbreak.brm)
    }
    ## Plot
    {
      ave.cots.interval <-
        cots.interval.outbreak.brm |>
        emmeans(~Decade, type = "response") |>
        as.data.frame() |>
        mutate(Decade = factor(Decade,
          levels = c(1980, 1990, 2000, 2010)
        ))
      cots2010 <- c(2010, 0, 0, 0)
      ave.cots.interval <- rbind(ave.cots.interval, cots2010) |>
        mutate(Decade = factor(Decade,
          levels = c(1980, 1990, 2000, 2010),
          labels = c("85 to 90", "91 to 00", "01 to 10", "11 to 20")
        ))

      modelled.cots.interval.plot <-
        ave.cots.interval |> ggplot(aes(x = Decade, y = prob)) +
        geom_bar(stat = "identity", fill = "seagreen") +
        geom_errorbar(
          stat = "identity",
          aes(ymin = lower.HPD, ymax = upper.HPD), width = 0
        ) +
        ggtitle("d)") +
        scale_y_continuous("Disturbance interval (years)",
          breaks = c(0, 2, 4, 6, 8, 10)
        ) +
        scale_x_discrete("Decade of first cots outbreak") +
        theme_classic()
      modelled.cots.interval.plot

      
      cots.interval.outbreak.box <-
        cots.interval.outbreak.brm |>
        emmeans(~Decade, type = "link") |>
        gather_emmeans_draws() |>
        mutate(.value = exp(.value)) |> 
        dplyr::select(-.draw, -.chain, -.iteration) |> 
        summarise_draws(median,
          q = ~ boxplot.stats(.x)$stats) 
    }
  }
  ## Combined interval plot
  {
    bleach.interval.est <-
      bleach.interval.est |>
      mutate(Disturbance = c("Bleaching", "Bleaching", "Bleaching", "Bleaching"))
    ave.cots.interval <- ave.cots.interval |>
      mutate(Disturbance = c("COTS", "COTS", "COTS", "COTS")) |>
      mutate(Decade = factor(Decade))
    cyclones.summary <-
      cyclones.summary |>
      mutate(Disturbance = c("Cyclone", "Cyclone", "Cyclone", "Cyclone"))

    dist_int <- rbind(
      bleach.interval.est,
      ave.cots.interval,
      cyclones.summary) |>
      as.data.frame()

    modelled.combined.interval.plot <-
      dist_int |>
      ggplot(aes(x = Decade, y = prob, fill = Disturbance)) +
      geom_bar(stat = "identity", width = 0.8,
        position = position_dodge(width = 0.8),
        colour = "black") +
      geom_errorbar(stat = "identity",
        aes(ymin = lower.HPD, ymax = upper.HPD), width = 0,
        position = position_dodge(width = 0.8)) +
      ggtitle("d)") +
      scale_y_continuous("Disturbance interval (years)") + # ,breaks=c(0,2,4,6,8,10))+
      scale_x_discrete("Decade of first disturbance") + # ,labels=c("85 to 90","91 to 00","01 to 10","11 to 20"))+
      # scale_fill_manual(values=c('firebrick3','seagreen','blue3'))+
      ## scale_fill_brewer(palette = "YlGnBu", type = "qual", direction = 1) +
      scale_fill_manual(breaks = c("Bleaching", "COTS", "Cyclone"),
        values = c(bleaching_palette[5], cots_palette[2], cyclone_palette[3]))+
      geom_text(aes(x = 0.75, y = 0, label = "**"), color = "black", size = 5) +
      geom_text(aes(x = 4, y = 0, label = "**"), color = "black", size = 5) +
      theme_classic() +
      theme(
        legend.position = c(0.8, 0.9),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.width = unit(0.25, "cm"),
        legend.key.height = unit(0.25, "cm")
      )
    modelled.combined.interval.plot

    ## cvd_grid(modelled.combined.interval.plot)
  }
  ## Alternative figure
  {
    dist_int <- 
      bleach.interval.box |>
      mutate(Disturbance = "Bleaching") |>
      rbind(cots.interval.outbreak.box |>
              mutate(Disturbance = "COTS")) |>
      rbind(cyclones.interval.box |>
              mutate(Disturbance = "Cyclone")) |> 
      mutate(Decade = case_when(
        Decade == "1980" ~ "85 to 90",
        Decade == "1990" ~ "91 to 00",
        Decade == "2000" ~ "01 to 10",
        Decade == "2010" ~ "11 to 20",
        Decade == "2020" ~ "21 to 30",
        .default = Decade
      )) |>
      filter(Decade != "21 to 30") |>
        mutate(
          Decade = factor(Decade,
            levels = c("85 to 90", "91 to 00", "01 to 10", "11 to 20")
          ),
          Disturbance = factor(Disturbance)
        ) |>
      ungroup()
    
    dist_int <-
      dist_int |>
      full_join(dist_int |>
                  expand(Decade, Disturbance))
    dist_int |> mutate(across(everything(), ~ replace(is.na, 0)))
    
    replace(dist_int, is.na(dist_int), 0)
    
      na_replace(0)
      replace(is.na(.), 0)
      replace_na(list(everything() = 0))
    

    modelled.combined.interval.plot <-
      dist_int |>
      ggplot(aes(x = Decade, fill = Disturbance)) +

  geom_boxplot(aes(x = Decade, y = q.3,
                   lower = q.2, upper = q.4,
                   middle = q.3,
                   ymin = q.1, ymax = q.5),
               stat = "Identity") +
      ggtitle("d)") +
      scale_y_continuous("Disturbance interval (years)") + # ,breaks=c(0,2,4,6,8,10))+
      scale_x_discrete("Decade of first disturbance") + # ,labels=c("85 to 90","91 to 00","01 to 10","11 to 20"))+
      # scale_fill_manual(values=c('firebrick3','seagreen','blue3'))+
      ## scale_fill_brewer(palette = "YlGnBu", type = "qual", direction = 1) +
      scale_fill_manual(breaks = c("Bleaching", "COTS", "Cyclone"),
        values = c(bleaching_palette[5], cots_palette[2], cyclone_palette[3]))+
      geom_text(aes(x = 0.75, y = 0, label = "**"), color = "black", size = 5) +
      geom_text(aes(x = 4, y = 0, label = "**"), color = "black", size = 5) +
      theme_classic() +
      theme(
        legend.position = c(0.8, 0.9),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.width = unit(0.25, "cm"),
        legend.key.height = unit(0.25, "cm")
      )
    modelled.combined.interval.plot

    ## cvd_grid(modelled.combined.interval.plot)
  }
  
}

##  combine into multi-panel fig for publications
 {
   multipanel.fig <-
     (plot3 | rel.coral.loss.thru.time.plot) /
     (perc.reefs.disturb | modelled.combined.interval.plot)

   multipanel.fig

   ggsave(
     filename = "../outputs/figures/figure_5.png",
     multipanel.fig,
     height = 20, width = 18, units = "cm", dpi = 600
   )
   ggsave(
     filename = "../outputs/figures/figure_5.pdf",
     multipanel.fig,
     height = 20, width = 18, units = "cm"
   )
   
 }
