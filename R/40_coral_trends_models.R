source("helper_functions.R")
cr_check_packages()

## NOTE: this script is provided so as to illustrate the steps taken
## to model coral cover. Access to the input data is not provided
## within this repository. Access to the primary data can be requested
## by emailing the senior author directly.

## LoadData
{
  load("../data/processed/manta.tow.RData")
  load("../data/processed/manta.sum.RData")
}

## Fit models
{
  {
    ## generate a list of reefs we are using to help accumulate other
    ## sources of associated data
    all.reefs <- manta.sum |>
      dplyr::select(P_CODE.mod, REEF_NAME, REEF_ID, Latitude, Longitude) |>
      group_by(REEF_NAME, REEF_ID) |>
      summarize(across(c(Latitude, Longitude), mean)) |>
      as.data.frame()
    write.csv(all.reefs,
      file = "../data/all.reefs_3Zone.csv",
      quote = FALSE, row.names = FALSE
    )
    ## Genuine stan cannot handle proportional data for binomial families
    ## (particularly when weights are applied). A work-around is to
    ## multiple the proportion by the weights and convert this into an integer
    dat.all <- manta.sum |>
      mutate(Location = Region) |>
      dplyr:::select(Cover, REEF_NAME, Tows, P_CODE.mod, Location,
        REPORT_YEAR) |>
      mutate(Year = factor(REPORT_YEAR),
        N = length(unique(REEF_NAME))) |>
      ungroup() |>
      mutate(
        Cvr1 = as.integer(as.vector(Cover) * Tows),
        Cvr0 = Tows - Cvr1
      )
  }
  ## Northern
  {
    ## Northern data
    {
      dat.all.northern = dat.all %>%
        filter(Location=='Northern GBR') %>%
        droplevels %>%
        mutate(P_CODE.mod=factor(ifelse(is.na(P_CODE.mod),'Other',P_CODE.mod))) %>%
        group_by(REEF_NAME) %>%
        mutate(W=mean(Tows, na.rm=TRUE)) %>%
        ungroup %>%
        mutate(W1=W/sum(W)) %>%
        group_by(Year) %>%
        mutate(W2=Tows/sum(Tows)) %>%
        ungroup
      save(dat.all.northern, file='../data/modelled/dat.all.northern.RData')
      ## Tow level data
      manta.tow.northern = manta.tow %>%
        filter(Region=='Northern GBR') %>%
        droplevels %>%
        mutate(oLIVE_CORAL=factor(LIVE_CORAL,
          levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'),
          ordered=TRUE),
          nREEF_NAME=as.numeric(as.factor(REEF_NAME)),
          nLIVE_CORAL=as.numeric(oLIVE_CORAL),
          REEF_YEAR = interaction(REEF_NAME, Year)
        )
      save(manta.tow.northern, file='../data/modelled/manta.tow.northern.RData')
    }
    ## ---- Northern.BRMS.tow.beta ry disp **
    {
      if ("BRMS beta ry disp" %in% models & "northern" %in% zone) {
        cat("Fitting brms ry disp for Northern\n\n")
        mod.northern_brms.beta.ry.disp <- brm(bf(Cover ~ Year + (1 | REEF_NAME / REEF_YEAR), phi ~ 0 + Year),
          data = manta.tow.northern,
          family = Beta(link = "logit"),
          iter = 1e4,
          warmup = 5e3,
          thin = 5,
          chains = 4, cores = 4,
          prior = prior(normal(0, 3), class = "b") +
            prior(normal(0, 3), class = "Intercept") +
            prior(gamma(2, 1), class = "sd") #+
          ## prior(gamma(2, 1), class = "phi")
        )
        ## ---- Northern.BRMS.tow.beta disp diagnostics
        {
          ## sampling diagnostics
          pdf(file = "../output/figures/traceplots_northern_brms.beta.ry.disp.pdf")
          rstan::traceplot(mod.northern_brms.beta.ry.disp$fit)
          dev.off()

          ## density overlay
          pdf(file = "../output/figures/density_northern_brms.beta.ry.disp.pdf")
          mod.northern_brms.beta.ry.disp %>% bayesplot::pp_check(type = "dens_overlay", ndraws = 100)
          dev.off()

          ## DHARMa residuals
          preds <- mod.northern_brms.beta.ry.disp %>%
            posterior_predict(nsamples = 250, summary = FALSE)
          mod.resids <- createDHARMa(
            simulatedResponse = t(preds),
            observedResponse = manta.tow.northern$Cover,
            fittedPredictedResponse = apply(preds, 2, median),
            integerResponse = FALSE
          )
          pdf(file = "../output/figures/DHARMa_northern_brms.beta.ry.disp.pdf")
          mod.resids %>% plot()
          dev.off()
          save(mod.resids, file = paste0("../data/modelled/resids.northern_brms.beta.ry.disp.RData"))
        }
        ## ----end
        dat.northern_brms.beta.ry.disp = emmeans(mod.northern_brms.beta.ry.disp,
          ~Year,
          type = "response"
        ) %>%
          as.data.frame()
        save(mod.northern_brms.beta.ry.disp, dat.northern_brms.beta.ry.disp, file = paste0("../data/modelled/mod.northern_brms.beta.ry.disp.RData"))
        rm(list = c("dat.northern_brms.beta.ry.disp", "mod.northern_brms.beta.ry.disp"))
        gc()
      }
    }
  }
  ## Central
  {
    # Central Data
    {
      dat.all.central = dat.all %>%
        filter(Location == "Central GBR") %>%
        droplevels() %>%
        mutate(P_CODE.mod = factor(ifelse(is.na(P_CODE.mod), "Other", P_CODE.mod))) %>%
        group_by(REEF_NAME) %>%
        mutate(W = mean(Tows, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(W1 = W / sum(W)) %>%
        group_by(Year) %>%
        mutate(W2 = Tows / sum(Tows)) %>%
        ungroup()
      save(dat.all.central, file = "../data/modelled/dat.all.central.RData")

      ## Tow level data
      manta.tow.central = manta.tow %>%
        filter(Region == "Central GBR") %>%
        droplevels() %>%
        mutate(
          oLIVE_CORAL = factor(LIVE_CORAL,
            levels = c("0", "1L", "1", "1U", "2L", "2", "2U", "3L", "3", "3U", "4L", "4", "4U", "5L", "5", "5U"),
            ordered = TRUE
          ),
          nREEF_NAME = as.numeric(as.factor(REEF_NAME)),
          nLIVE_CORAL = as.numeric(oLIVE_CORAL),
          REEF_YEAR = interaction(REEF_NAME, Year)
        )
      save(manta.tow.central, file = "../data/modelled/manta.tow.central.RData")
    }
    ## Central.BRMS.tow.beta ry disp **
    {
      if ('BRMS beta ry disp' %in% models & 'central' %in% zone) {
        cat('Fitting brms ry disp for Central\n\n')
        mod.central_brms.beta.ry.disp <- brm(bf(Cover ~ Year + (1|REEF_NAME/REEF_YEAR), phi~0+Year),
          data=manta.tow.central,
          family=Beta(link='logit'),
          iter=1e4,
          warmup=5e3,
          thin=5,
          chains=4, cores=4,
          prior = prior(normal(0, 3), class = "b") +
            prior(normal(0, 3), class = "Intercept") +
            prior(gamma(2, 1), class = "sd") #+
          ## prior(gamma(2, 1), class = "phi")
        )
        dat.central_brms.beta.ry.disp = emmeans(mod.central_brms.beta.ry.disp, ~Year, type='response') %>%
          as.data.frame()
        ## ---- Central.BRMS.tow.beta ry disp diagnostics
        {
          ## sampling diagnostics
          pdf(file = '../output/figures/traceplots_central_brms.beta.ry.disp.pdf')
          rstan::traceplot(mod.central_brms.beta.ry.disp$fit)
          dev.off()
          
          ## density overlay
          pdf(file = '../output/figures/density_central_brms.beta.ry.disp.pdf')
          mod.central_brms.beta.ry.disp %>% bayesplot::pp_check(type = "dens_overlay", ndraws = 100)
          dev.off()
          
          ## DHARMa residuals
          preds <- mod.central_brms.beta.ry.disp %>%
            posterior_predict(ndraws = 250, summary = FALSE)
          mod.resids <- createDHARMa(
            simulatedResponse = t(preds),
            observedResponse = manta.tow.central$Cover,
            fittedPredictedResponse = apply(preds, 2, median),
            integerResponse = FALSE
          )
          pdf(file = '../output/figures/DHARMa_central_brms.beta.ry.disp.pdf')
          mod.resids %>% plot()
          dev.off()
          save(mod.resids, file=paste0('../data/modelled/resids.central_brms.beta.ry.disp.RData'))
        }
        ## ----end
        save(mod.central_brms.beta.ry.disp, dat.central_brms.beta.ry.disp, file=paste0('../data/modelled/mod.central_brms.beta.ry.disp.RData'))
        rm(list=c('dat.central_brms.beta.ry.disp', 'mod.central_brms.beta.ry.disp'))
        gc()
      }
    }
  }
  ## Southern
  {
    ## Southern Data
    {
      dat.all.southern = dat.all %>%
        filter(Location=='Southern GBR') %>%
        droplevels %>%
        mutate(P_CODE.mod=factor(ifelse(is.na(P_CODE.mod),'Other',P_CODE.mod))) %>%
        group_by(REEF_NAME) %>%
        mutate(W=mean(Tows, na.rm=TRUE)) %>%
        ungroup %>%
        mutate(W1=W/sum(W)) %>%
        group_by(Year) %>%
        mutate(W2=Tows/sum(Tows)) %>%
        ungroup
      save(dat.all.southern, file='../data/modelled/dat.all.southern.RData')
      
      ## Tow level data
      manta.tow.southern = manta.tow %>%
        filter(Region=='Southern GBR') %>%
        droplevels %>%
        mutate(oLIVE_CORAL=factor(LIVE_CORAL,
          levels=c('0','1L','1','1U','2L','2','2U','3L','3','3U','4L','4','4U','5L','5','5U'),
          ordered=TRUE),
          nREEF_NAME=as.numeric(as.factor(REEF_NAME)),
          nLIVE_CORAL=as.numeric(oLIVE_CORAL),
          REEF_YEAR = interaction(REEF_NAME, Year)
        )
      save(manta.tow.southern, file='../data/modelled/manta.tow.southern.RData')
    }
    ## ---- Southern.BRMS.tow.beta disp ry **
    {
      if ('BRMS beta ry disp' %in% models & 'southern' %in% zone) {
        priors <- prior(normal(0, 3), class = "b") +
          prior(normal(0, 3), class = "Intercept") +
          prior(gamma(2, 1), class = "sd")
        ## The above priors where 0,1  0,1  2,1
        ## might like to try 0,2 0,1.5 2,1
        mod.southern_brms.beta.ry.disp <- brm(bf(Cover ~ Year + (1|REEF_NAME/REEF_YEAR), phi~0+Year),
          data=manta.tow.southern,
          family=Beta(link='logit'),
          iter=1e4,
          warmup=5e3,
          thin=5,
          chains=4, cores=4,
          prior = priors
          ## prior = prior(normal(0, 3), class = "b") +
          ##     prior(normal(0, 3), class = "Intercept") +
          ##     prior(gamma(2, 1), class = "sd") #+
          ## ## prior(gamma(2, 1), class = "phi")
        )
        dat.southern_brms.beta.ry.disp = emmeans(mod.southern_brms.beta.ry.disp, ~Year, type='response') %>%
          as.data.frame()
        ## ---- Southern.BRMS.tow.beta ry disp diagnostics
        {
          ## sampling diagnostics
          pdf(file = '../output/figures/traceplots_southern_brms.beta.ry.disp.pdf')
          rstan::traceplot(mod.southern_brms.beta.ry.disp$fit)
          dev.off()
          
          ## density overlay
          pdf(file = '../output/figures/density_southern_brms.beta.ry.disp.pdf')
          mod.southern_brms.beta.ry.disp %>% bayesplot::pp_check(type = "dens_overlay", ndraws = 100)
          dev.off()
          
          ## DHARMa residuals
          preds <- mod.southern_brms.beta.ry.disp %>%
            posterior_predict(ndraws = 250, summary = FALSE)
          mod.resids <- createDHARMa(
            simulatedResponse = t(preds),
            observedResponse = manta.tow.southern$Cover,
            fittedPredictedResponse = apply(preds, 2, median),
            integerResponse = FALSE
          )
          pdf(file = '../output/figures/DHARMa_southern_brms.beta.ry.disp.pdf')
          mod.resids %>% plot()
          dev.off()
          save(mod.resids, file=paste0('../data/modelled/resids.southern_brms.beta.ry.disp.RData'))
        }
        ## ----end
        save(mod.southern_brms.beta.ry.disp, dat.southern_brms.beta.ry.disp, file=paste0('../data/modelled/mod.southern_brms.beta.ry.disp.RData'))
        rm(list=c('dat.southern_brms.beta.ry.disp', 'mod.southern_brms.beta.ry.disp'))
        gc()
      }
    }
  }
}


