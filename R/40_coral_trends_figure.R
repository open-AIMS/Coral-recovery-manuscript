source("helper_functions.R")
cr_check_packages()

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


## Plots
{
  ## BRMS tow, beta ry disp
  load(file='../data/modelled/mod.northern_brms.beta.ry.disp.RData')
  load(file='../data/modelled/mod.central_brms.beta.ry.disp.RData')
  load(file='../data/modelled/mod.southern_brms.beta.ry.disp.RData')
  



  mceiling <- function(x,base){
    base*ceiling(x/base)
  }
  final_year_seq <- mceiling(final_year,5)

  ## Number of reefs
  nd <- manta.tow %>% group_by(Region) %>%
    summarise(Year=mean(range(as.numeric(as.character(Year)))),
      N=paste0('N=',length(unique(REEF_NAME)))) %>%
    bind_rows(manta.tow %>% summarise(Year=mean(range(as.numeric(as.character(Year)))),
      N=paste0('N=', length(unique(REEF_NAME)))) %>%
        mutate(Region='GBR'))


  model_source = 'brms.beta.ry.disp'
  include_n=FALSE
  include_gbr <- FALSE


  ## ---- defineData
  {
    
    ##dat.gbr <- sym(paste0('dat.gbr_',model_source))
    dat.northern <- sym(paste0('dat.northern_',model_source))
    dat.central <- sym(paste0('dat.central_',model_source))
    dat.southern <- sym(paste0('dat.southern_',model_source))
    
    newdata =
      #dat.gbr %>% eval %>% mutate(Region='GBR') %>%
      #rbind(dat.northern %>% eval %>% mutate(Region='Northern GBR')) %>%
      dat.northern %>% eval %>% mutate(Region='Northern GBR') %>%
      rbind(dat.central %>% eval %>% mutate(Region='Central GBR')) %>%
      rbind(dat.southern %>% eval %>% mutate(Region='Southern GBR')) %>%
      mutate(Region=factor(Region, levels=unique(Region))) %>%
      rename_with(recode, lower.HPD = 'lower', upper.HPD='upper',
        lower.CL = 'lower', upper.CL = 'upper',
        conf.low = 'lower', conf.high = 'upper',
        mean='response', estimate='response')
    if (!include_gbr) newdata <- newdata %>% filter(Region!='GBR') %>% droplevels()
    write_csv(newdata, file=paste0('../data/modelled/modelled_',model_source,'.csv'))
    ## ----end
  }


  {
    ## ---- initalPlot
    g1<-ggplot(newdata, aes(y = response, x = as.numeric(as.character(Year))))+
      geom_blank(aes(y=0.10,x=1995))+geom_blank(aes(y=0.35,x=1995))+
      facet_wrap(~Region, nrow = 3, scales='free_x',
        labeller = labeller(Region = setNames(paste0("\n", levels(newdata$Region),"\n"), levels(newdata$Region))))+
      geom_blank()+
      geom_pointrange(aes(ymin=lower, ymax=upper))+
      geom_line(aes(x = as.numeric(as.character(Year))), color='blue') +
      scale_y_continuous(expression(Coral~cover~('%')),labels=function(x) x*100, expand=c(0,0),limits=c(0,0.50)) +
      scale_x_continuous('',breaks=seq(1985,final_year_seq,by=5), limits=c(1985,final_year))+
      theme_classic()+
      theme(strip.background=element_rect(fill=hues[2], color='black', size=0.5),
        panel.background=element_rect(color='black'),
        axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
        axis.text.x=element_text(size=rel(1.5)), axis.title.x=element_blank(),
        axis.text.y=element_text(size=rel(1.5)),
        panel.grid.minor=element_line(size=0.1,color=NA),
        panel.grid.major=element_line(size=0.1,color='gray70'),
        panel.grid.minor.x=element_line(size=0,color='white',linetype=NULL),
        panel.grid.major.x=element_line(size=0,color='white',linetype=NULL),
        strip.text=element_text(margin=margin(t=0.1, b=0.1,unit='lines'),size=20,lineheight=0.5, face='bold',hjust=0.6,vjust=0.5),
        plot.margin=unit(c(0,0,2,0),'pt'),
        panel.spacing.x=unit(10,'pt'))
    if (include_n)
      g1 <- g1 + geom_text(data=nd, aes(y=Inf,x=Year, label=N), vjust=1.2)
    g1
    ggsave(file=paste0('../output/figures/threePanels.Bars_Stacked_',model_source,'_',ifelse(include_n,'with_n',''),'.pdf'), g1, width=5, height=9, units='in',dpi=300)
    ## ----end
    ## ---- addBanner
    if(!include_gbr) gt1=gt2; gt2=gt3; gt3=gt4;
    gT <- ggplot_gtable(ggplot_build(g1))
    facets <- grep("strip-t-1-1", gT$layout$name)
    gg <- with(gT$layout[facets,],
      gtable_add_grob(gT, ggplotGrob(gt1),t=t, l=5, b=b, r=6, name="pic_predator"))
    facets <- grep("strip-t-1-2", gT$layout$name)
    gg <- with(gg$layout[facets,],
      gtable_add_grob(gg, ggplotGrob(gt2),t=t, l=5, b=b, r=6, name="pic_predator"))
    facets <- grep("strip-t-1-3", gT$layout$name)
    gg <- with(gg$layout[facets,],
      gtable_add_grob(gg, ggplotGrob(gt3),t=t, l=5, b=b, r=6, name="pic_predator"))
    grid.draw(gg)
    ## ----end
    ## ---- savePlot
    save(gg, file=paste0('../data/spatial/threePanels.Bars_Stacked_',model_source,'.RData'))
    ggsave(file=paste0('../output/figures/threePanels.Bars_Stacked_',model_source,'_',ifelse(include_n,'with_n',''),'.pdf'), gg, width=5, height=9, units='in',dpi=300)
    ggsave(file=paste0('../output/figures/threePanels.Bars_Stacked_',model_source,'_',ifelse(include_n,'with_n',''),'.png'), gg, width=5, height=9, units='in',dpi=300)
    png(file=paste0('../output/figures/threePanels.Bars_Stacked_',model_source,'_',ifelse(include_n,'with_n',''),'.png'), width=5, height=9, units='in', res=300)
    grid.draw(gg)
    dev.off()
    ## ----end
    ## ----end





    if (!include_gbr) newdata <- newdata %>% filter(Region!='GBR') %>% droplevels()


    {
      g1 <-
        ggplot(newdata, aes(y = response, x = as.numeric(as.character(Year)))) +
        facet_wrap(~Region, nrow = 3, scales="free_x",
          labeller = labeller(Region = setNames(paste0("\n",
            levels(newdata$Region),"\n"),
            levels(newdata$Region)))) +
        geom_blank() +
        geom_pointrange(aes(ymin = lower, ymax = upper)) +
        geom_line(aes(x = as.numeric(as.character(Year))), color = "blue") +
        scale_y_continuous(expression(Coral~cover~("%")),
          labels = function(x) x * 100,
          expand = c(0, 0),
          limits = c(0, 0.50)) +
        scale_x_continuous("",breaks = seq(1985, final_year_seq, by = 5),
          limits = c(1985, final_year)) +
        theme_classic()+
        theme(strip.background = element_rect(fill = hues[2],
          color = "black", size = 0.5),
          panel.background = element_rect(color = "black"),
          axis.title.y = element_text(size = rel(1.5),
            margin = margin(r = 1,unit = "lines")),
          axis.text.x = element_text(size = rel(1.5)),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = rel(1.5)),
          panel.grid.minor = element_line(size = 0.1, color = NA),
          panel.grid.major = element_line(size = 0.1, color = "gray70"),
          panel.grid.minor.x = element_line(size = 0, color = "white",
            linetype = NULL),
          panel.grid.major.x = element_line(size = 0,color = "white",
            linetype = NULL),
          strip.text = element_text(margin = margin(t = 0.1, b = 0.1, unit = "lines"),
            size = 20,
            lineheight = 0.5,
            face = "bold",
            hjust = 0.6,
            vjust = 0.5),
          plot.margin = unit(c(0, 0, 2, 0), "pt"),
          panel.spacing.x = unit(10, "pt"))
      if (include_n)
        g1 <- g1 + geom_text(data = nd, aes(y = Inf,x = Year, label = N), vjust = 1.2)
      g1
    }
    {
      
      if(!include_gbr) gt1=gt2; gt2=gt3; gt3=gt4;
      gT <- ggplot_gtable(ggplot_build(g1))
      facets <- grep("strip-t-1-1", gT$layout$name)
      gg <- with(gT$layout[facets,],
        gtable_add_grob(gT, ggplotGrob(gt1),t=t, l=5, b=b, r=6, name="pic_predator"))
      facets <- grep("strip-t-1-2", gT$layout$name)
      gg <- with(gg$layout[facets,],
        gtable_add_grob(gg, ggplotGrob(gt2),t=t, l=5, b=b, r=6, name="pic_predator"))
      facets <- grep("strip-t-1-3", gT$layout$name)
      gg <- with(gg$layout[facets,],
        gtable_add_grob(gg, ggplotGrob(gt3),t=t, l=5, b=b, r=6, name="pic_predator"))
      grid.draw(gg)
      ## ----end
      ## ---- savePlot
      save(gg, file=paste0("../data/spatial/threePanels.Bars_Stacked_",model_source,".RData"))

      ggsave(file = paste0("../output/figures/threePanels.Bars_Stacked_", model_source, "_", ifelse(include_n, "with_n", ""), ".pdf"),
        gg,
        width = 5, height = 9, units = "in", dpi =30o)
      ggsave(
        file = paste0("../output/figures/threePanels.Bars_Stacked_", model_source, "_", ifelse(include_n, "with_n", ""), ".png"),
        gg,
        width = 5, height = 9, units = "in", dpi = 300
      )

      png(
        file = paste0("../output/figures/threePanels.Bars_Stacked_", model_source, "_", ifelse(include_n, "with_n", ""), ".png"),
        width = 5, height = 9, units = "in", res = 300
      )
      grid.draw(gg)
      dev.off()
    }

  }
