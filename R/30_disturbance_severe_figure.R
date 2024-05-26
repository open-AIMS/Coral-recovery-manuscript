source("helper_functions.R")
cr_check_packages()

finalYear <- 2022

## Start with table of disturbances
## ---- LoadData
{
  load(file='../data/modelled/bleaching.full_3Zone.RData')
  load(file='../data/modelled/cots.full_3Zone.RData')
  load(file='../data/modelled/cyclones.full_3Zone.RData')
  load(file='../data/processed/all.reefs.cyclones.RData')
  load(file = "../data/modelled/cots.sum.all_3Zone.RData")
}
## ----end


labs <- cots.sum.all |>
  pull(Location) |>
  levels()
labs.shorter <- gsub(' GBR', '', labs)

make_all_banners()

## New analysis
## ---- helperFunctions
{
  ## ---- calcFreqs
  {
    calcFreqs <- function(mod) {
      l1 <- emmeans(mod,
        specs = "time", by = "Zone",
        type = "response", at = list(time = 0)
      )
      list(Intercepts = l1)
      ##ideally, it would be good to be able to nominate the family to use here
      ##since the backtransform for slopes would just be exp rather than the
      ##inverse of a logit.  Add 1 to the slopes and intervals....
      ##SLOPES ARE ON A log odd-ratio scale.  WE SHOULD exp them so that they represent
      ## the factor change per year (or if them multiply by 100, the percent change per year)
      #lt=lstrends(mod, specs="Zone", var="time")
      lt = emtrends(mod, specs = "Zone", var = "time")
      l2 = test(lt)
      list(
        Intercepts = as.data.frame(summary(l1)) |>
          full_join(test(l1)) |>
          mutate(p.value = round(p.value, 3)),
        slopes = as.data.frame(summary(lt)) |>
          full_join(l2) |>
          mutate(p.value = round(p.value, 3))
      )
    }
  }
  ## ----end
  ## ---- calcFreqs.matrix
  {
    calcFreqs.matrix <- function(mod) {
      dat <- recover.data.glmmTMB(mod)
      form <- formula(delete.response(terms(mod)))
      ## Slopes - actually rates (change in probability of being impacted per year)
      newdata <- data.frame(
        time = 1,
        Zone = factor(levels(dat$Zone), levels = levels(dat$Zone))
      )
      Xmat <- model.matrix(form, data = newdata)
      Xmat[, c(1, 3, 4)] = 0
      coefs <- fixef(mod)[[1]]
      (fit = as.vector(coefs %*% t(Xmat)))
      SE = sqrt(diag(Xmat %*% vcov(mod)[[1]] %*% t(Xmat)))
      q = qnorm(0.975) #asymptotic (z test)
      l1 = data.frame(
        fit = exp(fit), lower = exp(fit - q * SE),
        upper = exp(fit + q * SE)
      )
      ## Intercepts - probabilty of being impacted at time 0 (1985)
      newdata <- data.frame(
        time = 0,
        Zone = factor(levels(dat$Zone), levels = levels(dat$Zone))
      )
      Xmat <- model.matrix(form, data = newdata)
      coefs <- fixef(mod)[[1]]
      fit <- as.vector(coefs %*% t(Xmat))
      SE <- sqrt(diag(Xmat %*% vcov(mod)[[1]] %*% t(Xmat)))
      q <- qnorm(0.975) #asymptotic (z test)
      l2 <- data.frame(
        fit = binomial()$linkinv(fit),
        lower = binomial()$linkinv(fit - q * SE),
        upper = binomial()$linkinv(fit + q * SE)
      )
      list(Intercept = l2, Slope = l1)
    }
  }
  ## ----end
  ## ---- ACF.glmmTMB
  {
    ACF.glmmTMB <- function (object, maxLag,
                             resType = c("pearson", "response", 
                               "deviance","raw"),
                             re=names(object$modelInfo$reTrms$cond$flist[1]),
                             ...) 
    {
      resType <- match.arg(resType)
      res <- resid(object, type = resType)
      res <- split(res,object$modelInfo$reTrms$cond$flist[[re]])
      if (missing(maxLag)) {
        maxL <- 30
        maxLag <- min(c(
          maxL = max(lengths(res)) - 1,
          as.integer(10 * log10(maxL + 1))
        ))
      }
      val <- lapply(res, function(el, maxLag) {
        N <- maxLag + 1L
        tt <- double(N)
        nn <- integer(N)
        N <- min(c(N, n <- length(el)))
        nn[1:N] <- n + 1L - 1:N
        for (i in 1:N) {
          tt[i] <- sum(el[1:(n - i + 1)] * el[i:n])
        }
        array(c(tt, nn), c(length(tt), 2))
      }, maxLag = maxLag)
      val0 <- rowSums(sapply(val, function(x) x[, 2]))
      val1 <- rowSums(sapply(val, function(x) x[, 1]))/val0
      val2 <- val1/val1[1L]
      z <- data.frame(lag = 0:maxLag, ACF = val2)
      attr(z, "n.used") <- val0
      class(z) <- c("ACF", "data.frame")
      z
    }
  }
  ## ----end
  ## ---- recover.data.glmmTMB
  {
  recover.data.glmmTMB <- function(object, ...) {
    fcall <- getCall(object)
    recover_data(fcall, delete.response(terms(object)),
      attr(model.frame(object), "na.action"), ...)
  }
  ## ----end
    ## ---- lsm.basis.glmmTMB
    lsm.basis.glmmTMB <- function (object, trms, xlev, grid, vcov.,
                                   mode = "asymptotic", component="cond", ...) {
      if (mode != "asymptotic") stop("only asymptotic mode is available")
      if (component != "cond") stop("only tested for conditional component")
      if (missing(vcov.)) 
        V <- as.matrix(vcov(object)[[component]])
      else V <- as.matrix(.my.vcov(object, vcov.))
      dfargs = misc = list()
      if (!is.null(object$modelInfo$family)) {
        fam = object$modelInfo$family$family
        misc$tran = object$modelInfo$family$link
        misc$inv.lbl = "response"
        if (!is.na(pmatch(fam, "binomial"))) 
          misc$inv.lbl = "prob"
        else if (!is.na(pmatch(fam, "poisson"))) 
          misc$inv.lbl = "rate"
      }
      #misc = lsmeans:::.std.link.labels(object$modelInfo$family, misc)
      if (mode == "asymptotic") {
        dffun = function(k, dfargs) NA
      }
      ## use this? misc = .std.link.labels(family(object), misc)
      contrasts = attr(model.matrix(object), "contrasts")
      m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
      X = model.matrix(trms, m, contrasts.arg = contrasts)
      bhat = fixef(object)[[component]]
      if (length(bhat) < ncol(X)) {
        kept = match(names(bhat), dimnames(X)[[2]])
        bhat = NA * X[1, ]
        bhat[kept] = fixef(object)[[component]]
        modmat = model.matrix(trms, model.frame(object), contrasts.arg = contrasts)
        nbasis = estimability::nonest.basis(modmat)
      }
      else nbasis = estimability::all.estble
      list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
        dfargs = dfargs, misc = misc,...)
    }
  }
  ## ----end
  ## ---- plotEffects
  {
    plotEffects <- function(mod, firstYear,individualReefs=NULL,
                            ribbonFillColor='blue', points = TRUE,
                            type = 1, maxYear = (finalYear)) {
      dat=recover.data.glmmTMB(mod)#mod$frame
      tt <- terms(mod)
      Terms <- delete.response(tt)

      newdata = NULL
      for (z in unique(dat$Zone)) {
        dat1 = dat %>% filter(Zone==z)
        pts = with(dat1,
          expand.grid(Zone=z,
            time=seq(min(time), max(time), len=100)))
        newdata = rbind(newdata,pts)
      }
      newdata = newdata %>%
        mutate(Zone=factor(Zone,
          levels=c('Northern GBR','Central GBR','Southern GBR')))
      m <- model.frame(Terms, newdata)
      Xmat <- model.matrix(Terms, m)
      coefs = fixef(mod)[[1]]
      fit=as.vector(coefs %*% t(Xmat))
      SE=sqrt(diag(Xmat %*% vcov(mod)[[1]] %*% t(Xmat)))
      q=qnorm(0.975) #asymptotic (z test)
      newdata = cbind(newdata, data.frame(fit=binomial()$linkinv(fit),
        lower=binomial()$linkinv(fit-q*SE),
        upper=binomial()$linkinv(fit+q*SE))) %>%
        mutate(Date=firstYear + time)

      if (!is.null(individualReefs)) {
        individualReefs = individualReefs %>%
          mutate(Date=firstYear + time, Time=factor(time))
        if (type ==1) {  
          g1=ggplot(newdata, aes(y=fit, x=Date)) +
            geom_blank()
          
          if(points) g1 <- g1 + geom_point(data=individualReefs, aes(y=fit), color='grey')
          g1 <- g1 + 
            geom_boxplot(data=individualReefs, aes(y=fit, group=Time), outlier.shape=NA) +
            facet_grid(~Zone) +
            geom_ribbon(aes(ymin=lower, ymax=upper), fill=ribbonFillColor, alpha=0.5, color=NA) +
            geom_line(color=ribbonFillColor) +
            scale_x_continuous('', expand=c(0,0)) + #, limits=c(1984,maxYear))+
            coord_cartesian(xlim=c(1984, (finalYear+1))) +
            scale_y_continuous('Pr(impact)', limits=c(0,1.00))+
            theme_classic()+
            theme(strip.background = element_blank(), plot.margin=unit(c(0,2,0,1),'lines'),
              panel.spacing=unit(1,'lines'))
        } else {  ## if need to transpose
          individualReefs <- individualReefs %>% mutate(Col = 1) %>%
            bind_rows(individualReefs %>% mutate(Col = 2)) %>%
            bind_rows(individualReefs %>% mutate(Col = 3))
          
          g1=ggplot(newdata, aes(y=fit, x=Date)) +
            geom_blank()
          if(points) g1 <- g1 + geom_point(data=individualReefs, aes(y=fit), color='grey')
          g1 <- g1 + 
            geom_boxplot(data=individualReefs, aes(y=fit, group=Time), outlier.shape=NA) +
            facet_grid(Zone~Col) +
            geom_ribbon(aes(ymin=lower, ymax=upper), fill=ribbonFillColor, alpha=0.5, color=NA) +
            geom_line(color=ribbonFillColor) +
            scale_x_continuous('', expand=c(0,0)) + #, limits=c(1984,2020))+
            coord_cartesian(xlim=c(1984, (finalYear+1))) +
            scale_y_continuous('Pr(impact)', limits=c(0,1.00))+
            theme_classic()+
            theme(strip.background = element_blank(), plot.margin=unit(c(0,2,0,1),'lines'),
              panel.spacing=unit(1,'lines'))
        }
      } else {
        g1=ggplot(newdata, aes(y=fit, x=Date)) +
          geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, color=NA) +
          geom_line() +
          facet_grid(~Zone) +
          scale_y_continuous('Pr(impact)')+
          theme_classic()
      }
      g1
    }
  }
  ## ----end
  ## ---- modelEachReef
  {
    modelEachReef = function(mod,form=CYCLONEany~time, dat=mod$frame) {
      df = list()
      for (i in 1:length(unique(dat$REEF_NAME))) {
        d1=dat %>% filter(REEF_NAME==unique(dat$REEF_NAME)[i]) %>% droplevels
        if (nrow(d1)>3) {
          newmod = glm(form,data=d1, family=binomial)
          d1=data.frame(REEF_NAME=unique(d1$REEF_NAME),
            Zone=unique(d1$Zone),
            time=seq(min(d1$time),max(d1$time),by=1))
          df[[i]]=cbind(d1, fit=predict(newmod, newdata=d1, type='response'))
        }
      }
      df = do.call('rbind',df)
    }
  }
  ## ----end
  ## ---- modelEachreefar1
  {
    modelEachReefAR1 = function(mod,form=CYCLONEany~time) {
      dat=recover.data.glmmTMB(mod) %>% full_join(mod$frame)
      wch=grep('poly',colnames(dat))
      dat=dat[,-wch]
      df = list()
      for (i in 1:length(unique(dat$REEF_NAME))) {
        d1=dat %>% dplyr::filter(REEF_NAME==unique(dat$REEF_NAME)[i]) %>% droplevels
        newmod = glmmTMB(form,data=d1, family=binomial)
        d1 = data.frame(
          REEF_NAME = unique(d1$REEF_NAME),
          Zone = unique(d1$Zone), time = seq(min(d1$time), max(d1$time), by = 1)
        )
        df[[i]]=cbind(d1, fit=predict(newmod, newdata=d1, type='response'))
      }
      df = do.call('rbind',df)
    }
  }
  ## ----end
  ## ---- format_plots
  format_plots <- function(gplot, sublabel) {
    g <-
      gplot +
      geom_text(data = data.frame(x = 1985, y = 1,
        Zone = factor("Northern GBR",
          levels = c("Northern GBR", "Central GBR", "Southern GBR")), l = sublabel),
        aes(y = y, x = x, label = l),
        hjust = 0, vjust = 1) +
      scale_x_continuous("", breaks = seq(1985, 2020, by = 5),
        position = "bottom", limits = c(1984.5, (finalYear + 1))) +
      coord_cartesian(xlim = c(1985, finalYear)) +
      scale_y_continuous(expression(phantom("(")), lim = c(0, 1.00)) +
      facet_grid(Zone ~Col, scales="fixed",
        labeller = labeller(Zone = setNames(paste0("", labs.shorter, "\n"), labs))) +
      theme(
        panel.border = element_rect(fill = NA, color = "black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = rel(1.0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = rel(1.2)),
        panel.grid.minor = element_line(size = 0.1, color = NA),
        panel.grid.major = element_line(size = 0.1, color = "gray70"),
        panel.grid.minor.x = element_line(size = 0.1, color = NA,
          linetype = "dashed"),
        panel.grid.major.x = element_line(size = 0.1, color = "gray70",
          linetype = "dashed"),
        plot.margin = unit(c(0, 5, 5, 6), "pt"),
        strip.text.x = element_blank(),
        strip.text = element_text(margin = margin(l = 0.5, r = 0.5,
          unit = "lines"), size = 20,
          lineheight = 0.5, face = "bold",
          hjust = 0.9, vjust = 0),
        strip.background = element_rect(fill = hues[2],
          color = "black", size = 0.5),
        plot.title.position = "panel")

    return(ggplot_gtable(ggplot_build(g)))
  }
  ## ----end
  ## ---- combine_panels
  combine_panels <- function(g.bleaching, g.cots, g.cyclones) {
    gplot = g.bleaching
    panels_src <- grep("panel-1-1", g.cots$layout$name)
    panels_target <- grep("panel-1-2", gplot$layout$name)
    gplot$grobs[panels_target] <- g.cots$grobs[panels_src]
    panels_src <- grep("panel-1-1", g.cyclones$layout$name)
    panels_target <- grep("panel-1-3", gplot$layout$name)
    gplot$grobs[panels_target] <- g.cyclones$grobs[panels_src]
    gplot %>% grid.draw()

    panels_src <- grep("panel-2-1", g.cots$layout$name)
    panels_target <- grep("panel-2-2", gplot$layout$name)
    gplot$grobs[panels_target] <- g.cots$grobs[panels_src]
    panels_src <- grep("panel-2-1", g.cyclones$layout$name)
    panels_target <- grep("panel-2-3", gplot$layout$name)
    gplot$grobs[panels_target] <- g.cyclones$grobs[panels_src]

    panels_src <- grep("panel-3-1", g.cots$layout$name)
    panels_target <- grep("panel-3-2", gplot$layout$name)
    gplot$grobs[panels_target] <- g.cots$grobs[panels_src]
    panels_src <- grep("panel-3-1", g.cyclones$layout$name)
    panels_target <- grep("panel-3-3", gplot$layout$name)
    gplot$grobs[panels_target] <- g.cyclones$grobs[panels_src]

    return(gplot)
  }

  ## ----end

  ## ---- add_strip_thumbs
  add_strip_thumbs <- function(gplot) {
    facets <- grep("strip-r-1", gplot$layout$name)
    gplot <- with(
      gplot$layout[facets, ],
      gtable_add_grob(gplot,
        ggplotGrob(banner_thumb_northern_v),
        t = 6, l = 10, b = 8, r = 10, name = "pic_predator"
      )
    )

    facets <- grep("strip-r-2", gplot$layout$name)
    gplot <- with(
      gplot$layout[facets, ],
      gtable_add_grob(gplot,
        ggplotGrob(banner_thumb_central_v),
        t = 10, l = 10, b = 10, r = 10, name = "pic_predator"
      )
    )
    facets <- grep("strip-r-3", gplot$layout$name)
    gplot <- with(
      gplot$layout[facets, ],
      gtable_add_grob(gplot,
        ggplotGrob(banner_thumb_southern_v),
        t = 12, l = 10, b = 12, r = 10, name = "pic_predator"
      )
    )
    return(gplot)
  }

## ----end

}
## ----end

dist.table <- list()

## Fit models
{
  ## Cyclones
  {
    ## Cyclone data
    {
      cyclones.full <-
        all.reefs |>
        left_join(cyclones.full) |>
        filter(!is.na(CYCLONEcat)) |>
        filter(REPORT_YEAR < (finalYear + 1))
      ## spread the data so that for each reef/year there is a binary response for each category
      cyclones.binary <- cyclones.full |>
        mutate(CYCLONEcat = as.factor(CYCLONEcat)) |>
        bind_cols() %>% 
        data.frame(model.matrix(~ -1 + CYCLONEcat, data = .)) |>
        mutate(
          CYCLONEany = ifelse(CYCLONEcat0 == 1, 0, 1),
          CYCLONEsevere = ifelse((CYCLONEcat2 + CYCLONEcat3) < 1, 0, 1),
          time = REPORT_YEAR - min(REPORT_YEAR),
          Zone = factor(Zone, levels = c("Northern GBR", "Central GBR", "Southern GBR"))
        )
    }
    ## Severe Cyclone models
    {
      dist.table[['Cyclones']] <- list()

      cyclones.severe.glmmTMB <-
        glmmTMB(CYCLONEsevere ~ time * Zone + (1 | REEF_NAME),
          data = cyclones.binary,
          family = binomial()
        )

      plot(ACF(cyclones.severe.glmmTMB, resType = "pearson"), alpha = 0.05)
      dist.table[['Cyclones']][['Intercept']] <-
        calcFreqs(cyclones.severe.glmmTMB)[[1]] %>%
        dplyr::select(-time) %>%
        mutate(Disturbance = 'Cyclones', Stat = 'Intercept') %>%
        dplyr::select(Disturbance, Stat, everything())
      dist.table[['Cyclones']][['Slope']] <- calcFreqs(cyclones.severe.glmmTMB)[[2]] %>%
        mutate(Disturbance = 'Cyclones', Stat = 'Slope') %>%
        dplyr::select(Disturbance, Stat, everything())
      dist.table[['Cyclones']][['PercentChange']] <-
        ((calcFreqs.matrix(cyclones.severe.glmmTMB)[[2]]-1)*100) %>%
        mutate(Disturbance = 'Cyclones', Stat = 'PercentChange') %>%
        dplyr::select(Disturbance, Stat, everything())  
      dist.table[['Cyclones']][['InterceptProb']] <-
        calcFreqs.matrix(cyclones.severe.glmmTMB)[[1]] %>%
        mutate(Disturbance = 'Cyclones', Stat = 'InterceptProb') %>%
        dplyr::select(Disturbance, Stat, everything())  

      df <- modelEachReef(cyclones.severe.glmmTMB,
        form = CYCLONEsevere ~ time
      )
      d1_2021 <- plotEffects(cyclones.severe.glmmTMB,
        firstYear = min(cyclones.binary$REPORT_YEAR),
        individualReefs = df,
        ribbonFillColor = brewer_pal(palette = "Blues")(4)[4],
        points = FALSE, type = 2
      )
      g.severe.cyclones_2021 <- d1_2021
      save(g.severe.cyclones_2021, file = '../data/modelled/g.severe.cyclones_2021.RData')
    }
  }

  ## COTS
  {
    ## COTS data
    {
      cots.full <-
        all.reefs |>
        left_join(cots.full |>
                    distinct()) |>
        filter(!is.na(COTScat)) |>
        filter(REPORT_YEAR < (finalYear + 1))
      ## spread the data so that for each reef/year there is a binary response for each category
      ##For some reason, some reefs (e.g. '16017S') have two different COTScat in a particular year (1994)
      ## To correct for this, I will give them the max category.
      cots.binary <-
        cots.full |>
        mutate(COTScat = factor(COTScat, levels = c("Zero", "NO", "IO", "AO"))) |>
        group_by(REEF_NAME, REEF_ID, Zone, Latitude, Longitude, REPORT_YEAR) |>
        summarize(COTScat = levels(COTScat)[max(as.numeric(COTScat))]) |>
        ungroup() |>
        mutate(COTScat = factor(COTScat, levels = c("Zero", "NO", "IO", "AO"))) |>
        bind_cols() %>%
        data.frame(model.matrix(~-1 + COTScat, data = .)) |>
        mutate(COTSany = ifelse(COTScatZero == 1, 0, 1),
          COTSsevere = ifelse((COTScatIO + COTScatAO) < 1, 0, 1),
          time = REPORT_YEAR - min(REPORT_YEAR),
          Time = as.factor(time),
          Zone = factor(Zone, levels = c("Northern GBR", "Central GBR", "Southern GBR")))
      save(cots.binary, file = "../data/modelled/cots.binary.RData")
    }
    ## CotsAO models
    {
      dist.table[['COTS']] <- list()

      ## cots.catAO.glmmTMB <-
      ##   glmmTMB(COTScatAO ~ time * Zone + (1 | REEF_NAME), # + ar1(-1+time|REEF_NAME),
      ##     data = cots.binary,
      ##     family = binomial()
      ##   )
      ## plot(ACF(cots.catAO.glmmTMB, resType = "pearson"), alpha = 0.05)

      ## dist.table[["COTS"]][["Intercept"]] <-
      ##   calcFreqs(cots.catAO.glmmTMB)[[1]] |>
      ##   dplyr::select(-time) |>
      ##   mutate(Disturbance = "COTS", Stat = "Intercept") |>
      ##   dplyr::select(Disturbance, Stat, everything())
      ## dist.table[["COTS"]][["Slope"]] <-
      ##   calcFreqs(cots.catAO.glmmTMB)[[2]] |>
      ##   mutate(Disturbance = "COTS", Stat = "Slope") |>
      ##   dplyr::select(Disturbance, Stat, everything())
      ## dist.table[["COTS"]][["PercentChange"]] <-
      ##   ((calcFreqs.matrix(cots.catAO.glmmTMB)[[2]]-1)*100) |>
      ##   mutate(Disturbance = "COTS", Stat = "PercentChange") |>
      ##   dplyr::select(Disturbance, Stat, everything())  
      ## dist.table[["COTS"]][["InterceptProb"]] <-
      ##   calcFreqs.matrix(cots.catAO.glmmTMB)[[1]] |>
      ##   mutate(Disturbance = "COTS", Stat = "InterceptProb") |>
      ##   dplyr::select(Disturbance, Stat, everything())  

      ## df <- modelEachReef(cots.catAO.glmmTMB,
      ##   form = COTScatAO ~ time
      ## )
      ## d1_2021 <- plotEffects(cots.catAO.glmmTMB,
      ##   firstYear = min(cots.binary$REPORT_YEAR),
      ##   individualReefs = df,
      ##   brewer_pal(palette = "Greens")(3)[3],
      ##   points = FALSE, type = 2
      ## )
      ## g.AO.cots_2021 <- d1_2021
      ## save(g.AO.cots_2021, file='../data/modelled/g.AO.cots_2021.RData')

      
      cots.severe.glmmTMB <-
        glmmTMB(COTSsevere ~ time * Zone + (1 | REEF_NAME), # + ar1(-1+time|REEF_NAME),
          data = cots.binary,
          family = binomial()
        )
      plot(ACF(cots.severe.glmmTMB, resType = "pearson"), alpha = 0.05)

      dist.table[["COTS"]][["Intercept"]] <-
        calcFreqs(cots.severe.glmmTMB)[[1]] |>
        dplyr::select(-time) |>
        mutate(Disturbance = "COTS", Stat = "Intercept") |>
        dplyr::select(Disturbance, Stat, everything())
      dist.table[["COTS"]][["Slope"]] <-
        calcFreqs(cots.severe.glmmTMB)[[2]] |>
        mutate(Disturbance = "COTS", Stat = "Slope") |>
        dplyr::select(Disturbance, Stat, everything())
      dist.table[["COTS"]][["PercentChange"]] <-
        ((calcFreqs.matrix(cots.severe.glmmTMB)[[2]]-1)*100) |>
        mutate(Disturbance = "COTS", Stat = "PercentChange") |>
        dplyr::select(Disturbance, Stat, everything())  
      dist.table[["COTS"]][["InterceptProb"]] <-
        calcFreqs.matrix(cots.severe.glmmTMB)[[1]] |>
        mutate(Disturbance = "COTS", Stat = "InterceptProb") |>
        dplyr::select(Disturbance, Stat, everything())  

      df <- modelEachReef(cots.severe.glmmTMB,
        form = COTSsevere ~ time
      )
      d1_2021 <- plotEffects(cots.severe.glmmTMB,
        firstYear = min(cots.binary$REPORT_YEAR),
        individualReefs = df,
        brewer_pal(palette = "Greens")(3)[3],
        points = FALSE, type = 2
      )
      g.severe.cots_2021 <- d1_2021
      save(g.severe.cots_2021, file='../data/modelled/g.severe.cots_2021.RData')
    }
  }

  ## Bleaching
  {
    ## Bleaching data
    {
      ## mike supplied some bleaching data to fill in the gaps between aerial surveys
      load(file <- "../data/modelled/bleaching.full_3Zone.RData")
      bleaching.full <-
        bleaching.full_3Zone |>
        filter(!is.na(Zone)) |>
        droplevels() |>
        filter(REPORT_YEAR < (finalYear + 1)) |>
        group_by(REEF_ID) |>
        arrange(REPORT_YEAR) |>
        ungroup() |>
        mutate(Zone = factor(Zone,
          levels = c("Northern","Central","Southern"),
          labels = c("Northern GBR","Central GBR","Southern GBR")))

      bleaching.binary <- bleaching.full |>
        filter(!is.na(BLEACHINGcat)) |>
        mutate(BLEACHINGcat = as.factor(BLEACHINGcat)) |>
        bind_cols() %>%
        data.frame(model.matrix(~-1 + BLEACHINGcat, data = .)) |>
        mutate(BLEACHINGany = ifelse(BLEACHINGcat0 == 1, 0, 1),
          BLEACHINGsevere = ifelse((BLEACHINGcat3 + BLEACHINGcat4 + BLEACHINGcat5) <1 ,0, 1),
          time = REPORT_YEAR - min(REPORT_YEAR),
          Zone = factor(Zone, levels = c("Northern GBR", "Central GBR", "Southern GBR")))
      save(bleaching.binary, file = "../data/modelled/bleaching.binary.RData")

      ## now a version that exludes LTMP data outside of December - March
      bleaching.binary.2 <-
        bleaching.full |>
        filter((Project == "LTMP" & Season == "Bleaching") | Project == "Aerial") |>
        filter(!is.na(BLEACHINGcat)) |>
        droplevels() |> 
        mutate(BLEACHINGcat = as.factor(BLEACHINGcat)) |>
        bind_cols() %>%
        data.frame(model.matrix(~-1 + BLEACHINGcat, data = .)) |>
        mutate(BLEACHINGany = ifelse(BLEACHINGcat0==1, 0,1),
          BLEACHINGsevere = ifelse((BLEACHINGcat3 + BLEACHINGcat4 + BLEACHINGcat5) <1, 0, 1),
          time = REPORT_YEAR - min(REPORT_YEAR),
          Zone = factor(Zone, levels = c("Northern GBR", "Central GBR", "Southern GBR")))
      save(bleaching.binary.2, file = "../data/modelled/bleaching.binary.2.RData")
    }
    ## Bleaching severe - This is the one used..
    {
      dist.table[['Bleaching']] <- list()
      bleaching.severe.glmmTMB <- glmmTMB(
        BLEACHINGsevere ~ time*Zone + (1|REEF_NAME),
        data = bleaching.binary, family = binomial())
      save(bleaching.severe.glmmTMB, file = '../data/modelled/bleaching.severe.glmmTMB.RData')

      summary(bleaching.severe.glmmTMB)
      plot(ACF(bleaching.severe.glmmTMB, resType="pearson"), alpha=0.05)

      dist.table[['Bleaching']][['Intercept']] <-
        calcFreqs(bleaching.severe.glmmTMB)[[1]] |>
        dplyr::select(-time) |>
        mutate(Disturbance = 'Bleaching', Stat = 'Intercept') |>
        dplyr::select(Disturbance, Stat, everything())
      dist.table[['Bleaching']][['Slope']] <- calcFreqs(bleaching.severe.glmmTMB)[[2]] |>
        mutate(Disturbance = 'Bleaching', Stat = 'Slope') |>
        dplyr::select(Disturbance, Stat, everything())  
      dist.table[['Bleaching']][['PercentChange']] <-
        ((calcFreqs.matrix(bleaching.severe.glmmTMB)[[2]]-1)*100) |>
        mutate(Disturbance = 'Bleaching', Stat = 'PercentChange') |>
        dplyr::select(Disturbance, Stat, everything())  
      dist.table[['Bleaching']][['InterceptProb']] <-
        calcFreqs.matrix(bleaching.severe.glmmTMB)[[1]] |>
        mutate(Disturbance = 'Bleaching', Stat = 'InterceptProb') |>
        dplyr::select(Disturbance, Stat, everything())  
      df <- modelEachReef(bleaching.severe.glmmTMB,
        form = BLEACHINGsevere ~ time
      )
      d1_2021 <- plotEffects(bleaching.severe.glmmTMB,
        firstYear = min(bleaching.binary$REPORT_YEAR),
        individualReefs = df,
        brewer_pal(palette = "Reds")(5)[5], points = FALSE, type = 2
      )
      g.severe.bleaching_2021 = d1_2021
      save(g.severe.bleaching_2021, file='../data/modelled/g.severe.bleaching_2021.RData')
    }
    ## ---- Bleaching severe version 2 (exclude LTMP not in summer) - This is the one used..
    {
      dist.table2 <- dist.table
      dist.table2[["Bleaching"]] <- list()
      bleaching.severe.2.glmmTMB <-
        glmmTMB(
          BLEACHINGsevere ~ time * Zone + (1 | REEF_NAME),
          data = bleaching.binary.2,
          family = binomial()
        )
      save(bleaching.severe.2.glmmTMB, file = "../data/modelled/bleaching.severe.2.glmmTMB.RData")

      summary(bleaching.severe.2.glmmTMB)
      plot(ACF(bleaching.severe.2.glmmTMB, resType="pearson"), alpha=0.05)

      dist.table2[["Bleaching"]][["Intercept"]] <-
        calcFreqs(bleaching.severe.2.glmmTMB)[[1]] %>%
        dplyr::select(-time) %>%
        mutate(Disturbance = "Bleaching", Stat = "Intercept") %>%
        dplyr::select(Disturbance, Stat, everything())
      dist.table2[["Bleaching"]][["Slope"]] <-
        calcFreqs(bleaching.severe.2.glmmTMB)[[2]] %>%
        mutate(Disturbance = "Bleaching", Stat = "Slope") %>%
        dplyr::select(Disturbance, Stat, everything())
      df <- modelEachReef(bleaching.severe.2.glmmTMB,
        form = BLEACHINGsevere ~ time
      )
      d1_2021 = plotEffects(bleaching.severe.2.glmmTMB,
        firstYear = min(bleaching.binary$REPORT_YEAR),
        individualReefs = df, brewer_pal(palette = "Reds")(5)[5], points = FALSE, type = 2
      )
      g.severe.bleaching_2021.2 = d1_2021
      save(g.severe.bleaching_2021, file="../data/modelled/g.severe.bleaching_2021.2.RData")
    }
    ## ----end
  }
}

## ---- Compilation Plots
{
  {
    hues <- RColorBrewer::brewer.pal(4, "Blues")
    load(file = '../data/modelled/g.severe.bleaching_2021.RData')
    load(file = '../data/modelled/g.severe.cots_2021.RData')
    load(file = '../data/modelled/g.severe.cyclones_2021.RData')

    g.severe.bleaching1 <- format_plots(g.severe.bleaching_2021, sublabel = "(a)")
    g.severe.bleaching1 |> grid.draw()

    g.severe.cots1 <- format_plots(g.AO.cots_2021, sublabel = "(b)")
    g.severe.cots1 |> grid.draw()

    g.severe.cyclones1 <- format_plots(g.severe.cyclones_2021, sublabel = "(c)")
    g.severe.cyclones1 |> grid.draw()

    gplot <- combine_panels(g.severe.bleaching1, g.severe.cots1, g.severe.cyclones1)
    gplot <- add_strip_thumbs(gplot)

    gplot <- patchwork::wrap_plots(
      textGrob(expression(Probability ~ of ~ being ~ impacted ~ phantom("(")),
        just = "centre", rot = 90, gp = gpar(fontsize = 16)
      ),
      wrap_plots(gplot, ncol = 1),
      widths = c(0.2, 10), ncol = 2
    )

    ggsave(
      file = "../outputs/figures/Disturbances_severe_compilation_newA_2021.pdf",
      gplot, width = 9, height = 7
    )
    ggsave(
      file = "../outputs/figures/Disturbances_severe_compilation_newA_2021.png",
      gplot, width = 9, height = 7, dpi = 600
    )
  }
}
