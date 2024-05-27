source("helper_functions.R")
cr_check_packages()


## Load the data
{
  
  load(file='../data/modelled/mod.northern_brms.beta.ry.disp.RData')
  load(file='../data/modelled/mod.central_brms.beta.ry.disp.RData')
  load(file='../data/modelled/mod.southern_brms.beta.ry.disp.RData')
}

## Settings
{
  ## BRMS tow, beta ry disp
  final_year <- 2022
  mceiling <- function(x, base) {
    base * ceiling(x / base)
  }
  final_year_seq <- mceiling(final_year, 5)

  ## Number of reefs
  nd <-
    manta.tow |>
    group_by(Region) |>
    summarise(
      Year = mean(range(as.numeric(as.character(Year)))),
      N = paste0("N=", length(unique(REEF_NAME)))
    ) |>
    bind_rows(manta.tow |> summarise(
      Year = mean(range(as.numeric(as.character(Year)))),
      N = paste0("N=", length(unique(REEF_NAME)))
    ) |>
      mutate(Region = "GBR"))

  model_source <- "brms.beta.ry.disp"
  include_n <- FALSE
  include_gbr <- FALSE

  hues <- RColorBrewer::brewer.pal(4, "Blues")
}

make_all_banners()

## Prepare Data
{
  dat.northern <- sym(paste0('dat.northern_',model_source))
  dat.central <- sym(paste0('dat.central_',model_source))
  dat.southern <- sym(paste0('dat.southern_',model_source))
  
  newdata <-
    dat.northern |>
    eval() |>
    mutate(Region = "Northern GBR") |>
    rbind(dat.central |> eval() |> mutate(Region = "Central GBR")) |>
    rbind(dat.southern |> eval() |> mutate(Region = "Southern GBR")) |>
    mutate(Region = factor(Region, levels = unique(Region))) |>
    rename_with(recode,
      lower.HPD = "lower", upper.HPD = "upper",
      lower.CL = "lower", upper.CL = "upper",
      conf.low = "lower", conf.high = "upper",
      mean = "response", estimate = "response"
    )
  if (!include_gbr) {
    newdata <- newdata |>
      filter(Region != "GBR") |>
      droplevels()
  }
  ## write_csv(newdata, file = paste0("../data/modelled/modelled_", model_source, ".csv"))
}

## Inital Plot
{
  g1 <-
    newdata |> 
    ggplot(aes(y = response, x = as.numeric(as.character(Year)))) +
    geom_blank(aes(y = 0.10, x = 1995)) +
    geom_blank(aes(y = 0.35, x = 1995)) +
    facet_wrap(~Region,
      nrow = 3, scales = "free_x",
      labeller = labeller(Region = setNames(paste0(
        "\n",
        levels(newdata$Region), "\n"
      ), levels(newdata$Region)))
    ) +
    geom_blank() +
    geom_pointrange(aes(ymin = lower, ymax = upper), size = rel(0.2)) +
    geom_line(aes(x = as.numeric(as.character(Year))), color = "blue") +
    scale_y_continuous(expression(Coral ~ cover ~ ("%")),
      labels = function(x) x * 100,
      expand = c(0, 0),
      limits = c(0, 0.50)) +
    scale_x_continuous("", breaks = seq(1985, final_year_seq, by = 5),
      limits = c(1985, final_year)) +
    theme_classic(base_size = 7) +
      theme(
        strip.background = element_rect(
          fill = hues[2], color = "black",
          size = 0.5
        ),
        panel.background = element_rect(color = "black"),
        axis.title.y = element_text(
          size = rel(1.5),
          margin = margin(r = 1, unit = "lines")
        ),
        axis.text.x = element_text(size = rel(1.5)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = rel(1.5)),
        panel.grid.minor = element_line(size = 0.1, color = NA),
        panel.grid.major = element_line(size = 0.1, color = "gray70"),
        panel.grid.minor.x = element_line(
          size = 0, color = "white",
          linetype = NULL
        ),
        panel.grid.major.x = element_line(
          size = 0, color = "white",
          linetype = NULL
        ),
        strip.text = element_text(
          margin = margin(t = 0.1, b = 0.1, unit = "lines"),
          ## size = 20, lineheight = 0.5, face = "bold", hjust = 0.6, vjust = 0.5
          size = rel(20/11), lineheight = 0.5, face = "bold", hjust = 0.6, vjust = 0.5
        ),
        plot.margin = unit(c(0, 0, 2, 0), "pt"),
        panel.spacing.x = unit(10, "pt")
      )
  
  if (include_n) {
    g1 <- g1 + geom_text(data = nd, aes(y = Inf, x = Year, label = N), vjust = 1.2)
  }
  g1
  ## ggsave(file = paste0("../output/figures/threePanels.Bars_Stacked_", model_source, "_", ifelse(include_n, "with_n", ""), ".pdf"), g1, width = 5, height = 9, units = "in", dpi = 300)
}

## Add Banner
{
  gT <- ggplot_gtable(ggplot_build(g1))
  facets <- grep("strip-t-1-1", gT$layout$name)
  gg <- with(
    gT$layout[facets, ],
    gtable_add_grob(gT,
      ggplotGrob(banner_thumb_northern),
      t = t, l = 5, b = b, r = 6, name = "pic_predator"
    )
  )
  facets <- grep("strip-t-1-2", gT$layout$name)
  gg <- with(
    gg$layout[facets, ],
    gtable_add_grob(gg,
      ggplotGrob(banner_thumb_central),
      t = t, l = 5, b = b, r = 6, name = "pic_predator"
    )
  )
  facets <- grep("strip-t-1-3", gT$layout$name)
  gg <- with(
    gg$layout[facets, ],
    gtable_add_grob(gg,
      ggplotGrob(banner_thumb_southern),
      t = t, l = 5, b = b, r = 6, name = "pic_predator"
    )
  )
  grid.draw(gg)
}

## Save Plot
{
  save(gg,
    file = paste0(
      "../data/spatial/threePanels.Bars_Stacked_",
      model_source, ".RData"
    )
  )

  ggsave(
    file = paste0(
      "../outputs/figures/threePanels.Bars_Stacked_",
      model_source, "_", ifelse(include_n, "with_n", ""), ".pdf"
    ),
    gg,
    width = 8.5, height = 8.5 * (9/5), units = "cm",
  )

  ggsave(
    file = paste0(
      "../outputs/figures/threePanels.Bars_Stacked_",
      model_source, "_", ifelse(include_n, "with_n", ""), ".png"
    ),
    gg,
    width = 8.5, height = 8.5 * (9/5), units = "cm", dpi = 600
  )
  ## png(
  ##   file = paste0(
  ##     "../outputs/figures/threePanels.Bars_Stacked_",
  ##     model_source, "_", ifelse(include_n, "with_n", ""), ".png"
  ##   ),
  ##   width = 5, height = 9, units = "in", res = 300
  ## )
  ## grid.draw(gg)
  ## dev.off()


  
}
