source("helper_functions.R")
cr_check_packages()

finalYear <- 2022

## Read in the data sources
{
  load(file = "../data/spatial/spatial_3Zone.RData")
  load(file = "../data/modelled/cots.sum.all_3Zone.RData")
  load(file = "../data/modelled/bleaching.sum.all_3Zone.RData")
  load(file = "../data/modelled/cyclones.sum.all_3Zone.RData")
}


hues <- RColorBrewer::brewer.pal(4, "Blues")
make_all_banners()
make_color_palettes()
## Generate the banner
{
  aus <- oz::ozRegion(sections = c(3, 11:13))
  aus <- oz::ozRegion()
  aus_df <- rbind(
    xy2df(aus$lines[[3]]),
    xy2df(aus$lines[[13]]),
    xy2df(aus$lines[[12]])[nrow(xy2df(aus$lines[[12]])):1, ],
    xy2df(aus$lines[[11]])
  )

  aus_ps <- SpatialPolygons(list(Polygons(list(Polygon(aus_df)),
    ID = "QLD"
  )))
  aus_sf <- aus_ps |>
    st_as_sf() |>
    st_set_crs(st_crs(4326))
  spatial_3Zone_sf <- spatial_3Zone |>
    st_as_sf() |>
    st_set_crs(st_crs(4326))

  banner_thumb <-
    ggplot(aus_sf) +
    geom_blank(aes(x = 190, y = -20)) +
    coord_map() +
    geom_sf(data = spatial_3Zone_sf, fill = NA, colour = "black", size = 0.2) +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = NA),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.background = element_blank(),
      panel.spacing = unit(0, "pt"),
      plot.margin = unit(c(0, 0, 0, 0), "pt")
    )

  banner_thumb_1 <-
    banner_thumb +
    geom_sf(data = spatial_3Zone_sf, fill = hues[4], colour = NA) +
    geom_sf(data = spatial_3Zone_sf, fill = NA, colour = "black") +
    geom_sf(fill = "white", colour = hues[4])
}
## Generate the individual banners
{
  banner_thumb_all <-
    banner_thumb +
    geom_sf(data = spatial_3Zone_sf, fill = hues[4], color = NA) +     
    geom_sf(data = spatial_3Zone_sf, fill = NA, color = "black", size = 0.2)  +
    geom_sf(fill = "white", color = hues[4]) +
    geom_blank(aes(x = 250, y = -20))

  banner_thumb_northern <-
    banner_thumb +
    geom_sf(data = spatial_3Zone_sf[1,], fill = hues[4], color = NA) +     
    geom_sf(data = spatial_3Zone_sf[1,], fill = NA, color = "black", size = 0.2)  +
    geom_sf(fill = "white", color = hues[4]) +
    geom_blank(aes(x = 250, y = -20))

  banner_thumb_central <-
    banner_thumb +
    geom_sf(data = spatial_3Zone_sf[2,], fill = hues[4], color = NA) +     
    geom_sf(data = spatial_3Zone_sf[2,], fill = NA, color = "black", size = 0.2)  +
    geom_sf(fill = "white", color = hues[4]) +
    geom_blank(aes(x = 250, y = -20))

  banner_thumb_southern <-
    banner_thumb +
    geom_sf(data = spatial_3Zone_sf[3,], fill = hues[4], color = NA) +     
    geom_sf(data = spatial_3Zone_sf[3,], fill = NA, color = "black", size = 0.2)  +
    geom_sf(fill = "white", color = hues[4]) +
    geom_blank(aes(x = 250, y = -20))

  save(banner_thumb_all,
    banner_thumb_northern,
    banner_thumb_central,
    banner_thumb_southern,
    file = "../data/spatial/banner_thumbs.RData"
  )
}


## COTs
{
  cots.sum.all <-
    cots.sum.all |>
    filter(Location != "Great Barrier Reef") |>
    filter(REPORT_YEAR < 2023) |> 
    droplevels()
  cots.sum.all <-
    cots.sum.all |>
    mutate(
      Zone = str_replace(Zone, "(.*) GBR", "\\1"),
      Location = factor(str_replace(Zone, "(.*) GBR", "\\1"),
        levels = c("Northern", "Central", "Southern")
      )
    )
  labs <-
    cots.sum.all |>
    pull(Location) |>
    levels()
  
  cots.dat <-
    cots.sum.all |>
    filter(
      REPORT_YEAR > 1985, REPORT_YEAR < 2025,
      COTScat %in% c("IO", "AO")
    )

  gcots <-
    cots.dat |>
    mutate(COTS.p = COTS.p / 100) |> 
    disturbance_bar_plot(
      var = COTS.p,
      var_cat = COTScat,
      xoffset = -0.3,
      bar_col = NA
      ## bar_col = "black"
    ) +
    scale_fill_manual("COTS outbreak status",
      breaks = c("IO", "AO"),
      labels = c("IO", "AO"),
      ## values = scales:::brewer_pal(palette = "Greens")(3)[-1])
      ## values = generate_tints(scales::brewer_pal(palette = "YlGnBu")(3)[1], 3)[-1]
      ## values = generate_tints(viridis_pal()(3)[3], 3)[-1]
      ## values = generate_tints(viridis_pal()(3)[2], 3)[-1]
      values = cots_palette
    ) 
  gcots
  
  ## gcots <-
  ##   cots.dat |> 
  ##   ggplot(aes(y = COTS.p, x = REPORT_YEAR - 0.3)) +
  ##   geom_bar(stat = "identity", position = "stack",
  ##     aes(fill = COTScat),
  ##     width = 0.3, show.legend = TRUE)+
  ##   gggrid::grid_panel(
  ##     grob = grobLabel,
  ##     mapping =  aes(label =  label),
  ##     data = cots.dat |>
  ##       mutate(label = case_when(
  ##         Zone == "Northern" ~ "b",
  ##         Zone == "Central" ~ "d",
  ##         Zone == "Southern" ~ "f"
  ##       ))
  ##   ) +
  ##   scale_fill_manual("COTS outbreak status", breaks = c("IO","AO"),
  ##     labels = c("IO", "AO"),
  ##     values = scales:::brewer_pal(palette = "Greens")(3)[-1]) +
  ##   facet_wrap(Location~., nrow = 3, scales = "free_x",
  ##     labeller = labeller(Location = setNames(paste0("", labs, "\n"), labs)),
  ##     strip.position = "right") +
  ##   scale_y_continuous(expression(Reefs~impacted~("%")),
  ##     expand = c(0, 0), lim = c(0, 100)) +
  ##   scale_x_continuous("", breaks = seq(1985, 2020, by = 5),
  ##     position = "bottom",
  ##     limits = c(1985, 2023.5)) +
  ##   theme_classic() +
  ##   theme(strip.background = element_rect(fill = hues[2],
  ##     color = "black", size = 0.5),
  ##     panel.border = element_rect(fill = NA, color = "black"),
  ##     axis.title.y = element_text(size = rel(1.5),
  ##       margin = margin(l = 0.5, r = 0.5, unit = "lines")),
  ##     axis.text.x = element_text(size = rel(1.5)),
  ##     axis.title.x = element_blank(),
  ##     axis.text.y = element_text(size = rel(1.5)),
  ##     panel.grid.minor = element_line(size = 0.1, color = NA),
  ##     panel.grid.major = element_line(size = 0.1, color = "gray70"),
  ##     panel.grid.minor.x = element_line(size = 0.1, color = NA,
  ##       linetype = "dashed"),
  ##     panel.grid.major.x = element_line(size = 0.1, color = "gray70",
  ##       linetype = "dashed"),
  ##     plot.margin = unit(c(0, 0, 2, 0), "pt"),
  ##     panel.spacing.x = unit(10, "pt"),
  ##     legend.position = "bottom",
  ##     strip.text = element_text(margin = margin(t = 0.5, b = 0.5,
  ##       l = 0.25, r = 0.25, unit = "lines"),
  ##       size=20,
  ##       lineheight = 0.5,
  ##       face = "bold",
  ##       hjust = 0.9,
  ##       vjust = -1)) +
  ##   guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))

  ## gcots1 <- gcots
  ## gcots

  ggcots = gcots+
    gggrid::grid_panel(
      grob = grobLabel,
      mapping =  aes(label =  label),
      data = cots.dat %>% mutate(label = case_when(
        Location == "Northern GBR" ~  "b", Location == "Central GBR" ~  "d", Location == "Southern GBR" ~ "e"))
    ) +
    facet_wrap(Location~., ncol =  1, strip.position =  "right", scales='free_x', labeller=labeller(Location=setNames(paste0("", labs.shorter, "\n"), labs))) +
    theme(
      panel.spacing.y = unit(5, 'pt'),
      plot.title =  element_blank(),
      axis.title.y=element_text(size=rel(1.5),
        margin=margin(r=1,unit='lines')),
      axis.text.x=element_text(size=rel(1.5)),
      axis.title.x=element_blank(),
      axis.text.y=element_text(size=rel(1.5)),
      plot.margin=unit(c(0,0,2,0),'pt'),
      panel.spacing.x=unit(1,'pt') 
    ) +
    guides(fill = guide_legend(title.position = 'top'))
  gcots
}

## Bleaching plot
{
  bleaching.sum.all <-
    bleaching.sum.all |>
    filter(Location !=  "Great Barrier Reef") |>
    filter(REPORT_YEAR < 2023) |> 
    droplevels()
  bleaching.dat <-
    bleaching.sum.all |>
    filter(REPORT_YEAR > 1985, BLEACHINGcat != "0") |>
    mutate(
      Zone = str_replace(Zone, "(.*) GBR", "\\1"),
      Location = factor(str_replace(Zone, "(.*) GBR", "\\1"),
        levels = c("Northern", "Central", "Southern")
      )
    )
  gbleaching <- 
    bleaching.dat |> 
    disturbance_bar_plot(
      var = BLEACHING.p,
      var_cat = BLEACHINGcat,
      xoffset = 0.3,
      bar_col = NA
      ## bar_col = "black"
    ) +
    scale_fill_manual("Bleaching severity",
      breaks = c(1, 2, 3, 4, 5),
      labels = c(1, 2, 3, 4, 5),
      ## values = generate_tints(scales::brewer_pal(palette = "YlGnBu")(3)[2], 6)[-1]
      ## values = c(scales::brewer_pal(palette = "Reds")(6)[-1])
      ## values = generate_tints(viridis_pal()(3)[1], 6)[-1]
      ## values = generate_tints(viridis_pal(option = "turbo")(3)[1], 6)[-1]
      values = bleaching_palette
    ) 
  gbleaching 

  ## gbleaching <-
  ##   bleaching.dat |>
  ##   ggplot(aes(y = 100 * BLEACHING.p, x = REPORT_YEAR + 0.3)) +
  ##   geom_bar(stat = "identity", position = "stack",
  ##     aes(fill = BLEACHINGcat), width = 0.3, show.legend = TRUE) +
  ##   scale_fill_manual("Bleaching severity", breaks = c(1, 2, 3, 4, 5),
  ##     labels = c(1, 2, 3, 4, 5),
  ##     values = c(scales::brewer_pal(palette = "Reds")(6)[-1])) +
  ##   facet_grid(Location~., scales = "fixed",
  ##     labeller = labeller(Location = setNames(paste0("", labs, "\n"), labs))) +
  ##   scale_y_continuous(expression(Reefs~impacted~("%")),
  ##     expand = c(0, 0),
  ##     lim = c(0, 100)) +
  ##   scale_x_continuous("", breaks = seq(1985, 2020, by = 5),
  ##     position  =  "bottom",
  ##     limits = c(1985, 2023.5)) +
  ##   theme_classic() +
  ##   theme(strip.background = element_rect(fill = hues[2], color = "black",
  ##     size = 0.5),
  ##     panel.border = element_rect(fill = NA, color = "black"),
  ##     axis.title.y = element_text(size = rel(1.5),
  ##       margin = margin(l = 0.5, r = 0.5, unit = "lines")),
  ##     axis.text.x = element_text(size = rel(1.0)),
  ##     axis.title.x = element_blank(),
  ##     axis.text.y = element_text(size = rel(1.0)),
  ##     panel.grid.minor = element_line(size = 0.1, color = NA),
  ##     panel.grid.major = element_line(size = 0.1, color = "gray70"),
  ##     panel.grid.minor.x = element_line(size = 0.1, color = NA,
  ##       linetype = "dashed"),
  ##     panel.grid.major.x = element_line(size = 0.1, color = "gray70",
  ##       linetype = "dashed"),
  ##     plot.margin = unit(c(0, 0, 2, 0), "pt"),
  ##     panel.spacing.x  =  unit(10, "pt"),
  ##     legend.position  =   "bottom",
  ##     strip.text = element_text(margin = margin(t = 0.5, b = 0.5,
  ##       l = 0.25, r = 0.25, unit = "lines"),
  ##       size = 20,
  ##       lineheight = 0.5,
  ##       face = "bold",
  ##       hjust = 0.9,
  ##       vjust = -1)) +
  ##   guides(fill  =  guide_legend(title.position = "top", title.hjust  =  0.5))
  ## gbleaching1 <- gbleaching
  ## gbleaching


  ggbleaching = gbleaching + theme(panel.background=element_blank(), legend.justification = 'left', legend.direction = 'horizontal') + #ggtitle('a)') +
    gggrid::grid_panel(
      grob = grobLabel,
      mapping =  aes(label =  label),
      data = bleaching.dat %>% mutate(label = case_when(
        Location == "Northern GBR" ~  "b", Location == "Central GBR" ~  "d", Location == "Southern GBR" ~ "e"))
    ) +
    facet_wrap(Location~., ncol =  1, strip.position =  "right", scales='free_x', labeller=labeller(Location=setNames(paste0("", labs.shorter, "\n"), labs))) +
    theme(
      panel.spacing.y = unit(5, 'pt'),
      plot.title =  element_blank(),
      axis.title.y=element_text(size=rel(1.5), margin=margin(r=1,unit='lines')),
      axis.text.x=element_text(size=rel(1.5)), axis.title.x=element_blank(),
      axis.text.y=element_text(size=rel(1.5)),
      plot.margin=unit(c(0,0,2,0),'pt'),
      panel.spacing.x=unit(10,'pt') 
    ) +
    guides(fill = guide_legend(title.position = 'top'))
  gbleaching
  
}

## Cyclones plot
{
  cyclones.sum.all <-
    cyclones.sum.all |>
    filter(Location != "Great Barrier Reef") |>
    filter(REPORT_YEAR < 2023) |> 
    droplevels()
  cyclones.dat <-
    cyclones.sum.all |>
    filter(REPORT_YEAR > 1985, CYCLONEcat != 0) |>
    mutate(
      Zone = str_replace(Zone, "(.*) GBR", "\\1"),
      Location = factor(str_replace(Zone, "(.*) GBR", "\\1"),
        levels = c("Northern", "Central", "Southern")
      )
    )
  gcyclones <-
    cyclones.dat |>
    mutate(CYCLONE.p = CYCLONE.p / 100) |> 
    disturbance_bar_plot(
      var = CYCLONE.p,
      var_cat = CYCLONEcat,
      xoffset = 0,
      bar_col = NA
      ## bar_col = "black"
    ) +
    scale_fill_manual("Cyclone severity",
      breaks = c(1, 2, 3),
      labels = c(1, 2, 3),
      ## values = generate_tints(scales::brewer_pal(palette = "YlGnBu")(3)[3], 4)[-1]
      ## values = c(scales::brewer_pal(palette = "Blues")(4)[-1])
      ## values = generate_tints(viridis_pal()(3)[2], 4)[-1]
      ## values = generate_tints(viridis_pal(option = "turbo")(3)[3], 4)[-1]
      values = cyclone_palette
    ) 

  ## gcyclones <-
  ##   cyclones.dat |>
  ##   ggplot(aes(y = CYCLONE.p, x = REPORT_YEAR + 0)) +
  ##   geom_bar(stat = "identity", position = "stack",
  ##     aes(fill = CYCLONEcat),
  ##     width = 0.3, show.legend = TRUE) +
  ##   scale_fill_manual("Cyclone severity", breaks = c(1, 2, 3),
  ##     labels = c(1, 2, 3),
  ##     values = c(scales::brewer_pal(palette = "Blues")(4)[-1])) +
  ##   facet_grid(Location~., scales = "fixed",
  ##     labeller = labeller(Location = setNames(paste0("", labs, "\n"), labs))) +
  ##   scale_y_continuous(expression(Reefs~impacted~("%")),expand = c(0,0),lim = c(0,100)) +
  ##   scale_x_continuous("", breaks = seq(1985, 2020, by = 5),
  ##     position  =  "bottom",
  ##     limits = c(1985, 2023.5)) +
  ##   theme_classic() +
  ##   theme(
  ##     strip.background = element_rect(fill = hues[2], color = "black",
  ##       size = 0.5),
  ##     panel.border = element_rect(fill = NA, color = "black"),
  ##     axis.title.y = element_text(size = rel(1.5),
  ##       margin = margin(l = 0.5, r = 0.5, unit = "lines")),
  ##     axis.text.x = element_text(size = rel(1.0)),
  ##     axis.title.x = element_blank(),
  ##     axis.text.y = element_text(size = rel(1.0)),
  ##     panel.grid.minor = element_line(size = 0.1, color = NA),
  ##     panel.grid.major = element_line(size = 0.1, color = "gray70"),
  ##     panel.grid.minor.x = element_line(size = 0.1, color = NA, linetype = "dashed"),
  ##     panel.grid.major.x = element_line(size = 0.1, color = "gray70", linetype = "dashed"),
  ##     plot.margin = unit(c(2, 5, 5, 0), "pt"),
  ##     panel.spacing.y  =  unit(15, "pt"),
  ##     legend.position  =  "bottom",
  ##     strip.text = element_text(margin = margin(t = 0.5, b = 0.5,
  ##       l = 0.25, r = 0.25, unit = "lines"),
  ##       size = 20,
  ##       lineheight = 0.5,
  ##       face = "bold",
  ##       hjust = 0.9,
  ##       vjust = -1)) +
  ##   guides(fill  =  guide_legend(title.position = "top", title.hjust  =  0.5))
  ## gcyclones1 <- gcyclones
  gcyclones
  
  ggcyclones = gcyclones + theme(panel.background=element_blank(), legend.justification = c(0.6, 0.5), legend.direction = 'horizontal') + ggtitle('a)')+
    facet_grid(Location~., scales='fixed', labeller=labeller(Location=setNames(paste0("", labs.shorter, "\n"), labs)))+
    theme(
      panel.spacing.y = unit(5, 'pt')) +
    guides(fill = guide_legend(title.position = 'top'))
}


## Combine each plot
{
  gbleaching <- gbleaching +
    theme(panel.background = element_blank())
  gcyclones <- gcyclones +
    theme(panel.background = element_blank())
  gcots <- ggplot_gtable(ggplot_build(gcots))
  gbleaching <- ggplot_gtable(ggplot_build(gbleaching))
  gcyclones <- ggplot_gtable(ggplot_build(gcyclones))
  panels <- grepl("panel", gbleaching$layout$name)

  pp <- c(subset(gcots$layout, grepl("panel", gcots$layout$name), se = t:r))

  # Overlap panels for second plot on those of the first plot
  gT <- gtable_add_grob(gcots, gbleaching$grobs[grepl("panel", gcots$layout$name)], 
    pp$t, pp$l, pp$b, pp$l, name='bleaching')
  gT <- gtable_add_grob(gT, gcyclones$grobs[grepl("panel", gcots$layout$name)], 
    pp$t, pp$l, pp$b, pp$l, name='cyclones')

  legCots <- gcots$grobs[[which(gcots$layout$name == "guide-box")]]
  legCots$widths[c(1,5)] <- unit(0.5, "cm")
  legBleaching <- gbleaching$grobs[[which(grepl("guide-box", gbleaching$layout$name))]]
  legBleaching$widths[c(1,5)] <- unit(0.5, "cm")
  legCyclones <- gcyclones$grobs[[which(grepl("guide-box", gcyclones$layout$name))]]
  legCyclones$widths[c(1,5)] <- unit(0.5, "cm")
  leg <- gtable:::cbind_gtable(
    gtable:::cbind_gtable(legBleaching, legCyclones, "first"),
    legCots, "first")
  gT$grobs[[which(gT$layout$name == "guide-box")]] <- leg
  grid::grid.draw(gT)  # Note: Legend does not fit

 }

{
  ## banner_thumb_v <-
  ##   ggplot(aus_sf) +
  ##   geom_sf(data = spatial_3Zone_sf, fill = NA, colour = "black", size = 0.2) +
  ##   theme_classic() +
  ##   theme(
  ##     panel.background = element_rect(fill = NA),
  ##     axis.text.y = element_blank(),
  ##     axis.text.x = element_blank(),
  ##     axis.title.y = element_blank(),
  ##     axis.title.x = element_blank(),
  ##     axis.ticks = element_blank(),
  ##     axis.line = element_blank(),
  ##     plot.background = element_blank(),
  ##     panel.spacing = unit(0, "pt"),
  ##     plot.margin = unit(c(0, 0, 0, 0), "pt")
  ##   )

  ## banner_thumb_northern_v <-
  ##   banner_thumb_v +
  ##   geom_sf(data = spatial_3Zone_sf[1, ], fill = hues[4], color = NA) +
  ##   geom_sf(data = spatial_3Zone_sf[1, ], fill = NA, color = "black") +
  ##   geom_sf(fill = "white", color = hues[4]) +
  ##   coord_sf() +
  ##   theme(plot.margin = unit(c(2, 0, 100, 0), "pt"))

  ## banner_thumb_central_v <-
  ##   banner_thumb_v +
  ##   geom_sf(data = spatial_3Zone_sf[2, ], fill = hues[4], color = NA) +
  ##   geom_sf(data = spatial_3Zone_sf[2, ], fill = NA, color = "black") +
  ##   geom_sf(fill = "white", color = hues[4]) +
  ##   coord_sf() +
  ##   theme(plot.margin = unit(c(2, 0, 100, 0), "pt"))

  ## banner_thumb_southern_v <-
  ##   banner_thumb_v +
  ##   geom_sf(data = spatial_3Zone_sf[3, ], fill = hues[4], color = NA) +
  ##   geom_sf(data = spatial_3Zone_sf[3, ], fill = NA, color = "black") +
  ##   geom_sf(fill = "white", color = hues[4]) +
  ##   coord_sf() +
  ##   theme(plot.margin = unit(c(2, 0, 100, 0), "pt"))

  facets <- grep("strip-r-1", gT$layout$name)
  c(subset(gg$layout, grepl("strip-r-1", gg$layout$name), se = t:r))
  gg <- with(gT$layout[facets,],
    gtable_add_grob(gT, ggplotGrob(banner_thumb_northern_v),t=7, l=6, b=7, r=6, name="pic_predator"))
  grid.draw(gg)

  facets <- grep("strip-r-2", gT$layout$name)
  c(subset(gg$layout, grepl("strip-r-2", gg$layout$name), se = t:r))
  gg <- with(gg$layout[facets,],
    gtable_add_grob(gg, ggplotGrob(banner_thumb_central_v), t=9, l=6, b=9, r=6, name="pic_predator"))

  facets <- grep("strip-r-3", gT$layout$name)
  c(subset(gg$layout, grepl("strip-r-3", gg$layout$name), se = t:r))
  gg_2021 <- with(gg$layout[facets,],
    gtable_add_grob(gg, ggplotGrob(banner_thumb_southern_v),t=11, l=6, b=11, r=6, name="pic_predator"))

  grid::grid.draw(gg_2021)
  gplot <- patchwork::wrap_plots(gg_2021)

  ggsave(
    ## file = "../outputs/figures/Disturbances_severe_compilation_newA_2021.pdf",
    file = "../outputs/figures/figure_2.png",
    gplot,
    width = 18, height = 18*1.05, units = "cm", dpi = 600
  )
  ggsave(
    ## file = "../outputs/figures/Disturbances_severe_compilation_newA_2021.pdf",
    file = "../outputs/figures/figure_2.pdf",
    gplot,
    width = 18, height = 18*1.05, units = "cm", dpi = 600
  )

  ## pdf(file='../outputs/figures/Disturbances_all_no_label_transposed.pdf', width = 7, height = 7*1.05)
  ## grid.draw(gg_2021)
  ## dev.off()

  ## ## resize to 8.3cm (3.27in) or 3.27x600 = 1960.63
  
  ## system("convert -density 600 ../outputs/figures/Disturbances_all_no_label_transposed.pdf -resize 1960.63x ../outputs/figures/Disturbances_all_no_label_transposed.pdf")
  

}


