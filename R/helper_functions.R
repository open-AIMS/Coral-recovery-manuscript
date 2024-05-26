###################################################################
## The following function checks to ensure that all the required ##
## packages are available on the system.                         ##
###################################################################
cr_check_packages <- function() {
  ## require(gdata) # load this first as it masks many
  require(tidyverse)
  require(viridis)
  require(RColorBrewer)
  require(gtable)
  require(grid)
  require(nlme)
  ## require(gridExtra)
  ## require(xtable)
  ## library(broom)
  ## require(rgdal)
  ## require(rgeos)
  require(sp)
  require(oz)
  ## require(maps)
  ## require(mapdata)
  ## require(ggsn)
  require(scales)
  ## require(mapping) # consider replacing this with a self contained function
  ## require(maptools)
  ## require(raster)

  ## require(INLA)
  ## require(rstanarm)
  ## require(coda)
  require(sf)
  require(patchwork)
  require(glmmTMB)
  ## require(brms)
  require(emmeans)
}


COTScategories = function(x) {
    case_when(
        x == 0 ~ 'Zero',
        x > 0 & x < 0.22 ~ 'NO',
        x >= 0.22 & x < 1 ~ 'IO',
        x >= 1 ~ 'AO'
    )
}

COE_COTScategories = function(x) {
    case_when(
        x == 0 ~ 0,
        x > 0 & x <= 0.1 ~1,
        x > 0.1 & x <= 0.3 ~2,
        x > 0.3 & x <= 0.6 ~3,
        x > 0.6 & x <= 0.9 ~4,
        x > 0.9 ~5
    )
}
COE_BLEACHINGcategories = function(x) {
    case_when(
        x < 0.025 ~ 0,
        ## x > 0 & x <= 0.1 ~1,
        x >= 0.025 & x <= 0.1 ~1,
        x > 0.1 & x <= 0.3 ~2,
        x > 0.3 & x <= 0.6 ~3,
        x > 0.6 & x <= 0.9 ~4,
        x > 0.9 ~5
    )
}

CoralTrends_calcPercent = function(x) {
    ifelse(x=='0', 0,
    ifelse(x=='1', 0.05,
    ifelse(x=='1L', 0.025,
    ifelse(x=='1U', 0.075,
    ifelse(x=='2', 0.2,
    ifelse(x=='2L', 0.15,
    ifelse(x=='2U', 0.25,
    ifelse(x=='3', 0.4,
    ifelse(x=='3L', 0.35,
    ifelse(x=='3U', 0.45,
    ifelse(x=='4', 0.625,
    ifelse(x=='4L', 0.5625,
    ifelse(x=='4U', 0.6875,
    ifelse(x=='5', 0.875,
    ifelse(x=='5L',0.8125,0.9375)))))))))))))))
}


ML_gClip <- function(shp, bb){
    if(class(bb) == "matrix") {
        if (identical(dim(bb), c(2L,2L))) {
            b_poly <- as(raster:::extent(as.vector(t(bb))), "SpatialPolygons")
        } else b_poly = bb
    } else if (class(bb) =='SpatialPolygons') {
        b_poly=bb
    } else b_poly <- as(raster:::extent(bb), "SpatialPolygons")
  gIntersection(shp, b_poly, byid = T)
}


######################################################################
## The following function generates Zone categories based on De'ath ##
## 2012's latitudinal divisions.                                    ##
##   parameters:                                                    ##
##      x:     a numeric vector of latitudes                        ##
##   returns:  a categorical vector of Zones                        ##
######################################################################
CoralTrends_calc3ZoneLocations <- function(x) {
    factor(ifelse(x<= -10.68 & x > -15.4, 'Northern',  #glenn's version is -11.8
           ifelse(x<= -15.4 & x > -20.0, 'Central',
           ifelse(x<= -20.0 & x > -23.92, 'Southern','Outside'))))
}
CoralTrends_calc3ZoneLocation <- function(dat) {
    load(paste0(DATA_PATH, 'primary/gbr_3Zone.RData'))
    dat %>%
        st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(gbr_3Zone)) %>%
        st_join(gbr_3Zone) %>%
        cbind(Longitude=st_coordinates(.)[,1],Latitude=st_coordinates(.)[,2]) %>%
        st_drop_geometry
}

## The following function converts a cover into a category representing the median of the category
median_cover_cat <- function(dat) {
    n <- length(dat)
    if (n %% 2 == 0) {
        d <- paste(unique(sort(dat)[c((n-1)/2, (n+1)/2)]), collapse='/')
    } else {
        d <- paste(unique(sort(dat)[(n+1)/2]))
    }
    factor(d)
}


################################################################
## The following function converts a list with x and y into a ##
## data.frame                                                 ##
##   parameters:                                              ##
##      xy:    a list with elements x and y                   ##
##   returns:  a data.frame with fields x and y               ##
################################################################
xy2df<-function(xy) {
  data.frame(x=xy$x,y=xy$y)
}


grobLabel <- function(data, coords) {
  grid::textGrob(data$label[which.max(coords$y)],
                     x=unit(0, "npc") + 1*grid::unit(0.5,"lines"),
                     y=unit(1, "npc") + -1*unit(0.5,"lines"),
                     gp =  gpar(fontsize =  15),
                     just=c("left","top"))
}

disturbance_bar_plot <- function(dat, var = NULL, var_cat = NULL, xoffset = 0) {
  var <- enquo(var)
  var_cat <- enquo(var_cat)
  
  dat |>
    ggplot(aes(y = 100 * {{var}}, x = REPORT_YEAR + xoffset)) +
    geom_bar(stat = "identity", position = "stack",
      aes(fill = {{var_cat}}), width = 0.3, show.legend = TRUE) +
    facet_grid(Location~., scales = "fixed",
      labeller = labeller(Location = setNames(paste0("", labs, "\n"), labs))) +
    scale_y_continuous(expression(Reefs~impacted~("%")),
      expand = c(0, 0),
      lim = c(0, 100)) +
    scale_x_continuous("", breaks = seq(1985, 2020, by = 5),
      position  =  "bottom",
      ## limits = c(1985, 2022.5)) +
      limits = c(1985, 2022.5),
      expand = c(0.02, 0)) +
    theme_classic() +
    theme(strip.background = element_rect(fill = hues[2], color = "black",
      size = 0.5),
      panel.border = element_rect(fill = NA, color = "black"),
      axis.title.y = element_text(size = rel(1.5),
        margin = margin(l = 0.5, r = 0.5, unit = "lines")),
      axis.text.x = element_text(size = rel(1.0)),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = rel(1.0)),
      panel.grid.minor = element_line(size = 0.1, color = NA),
      panel.grid.major = element_line(size = 0.1, color = "gray70"),
      panel.grid.minor.x = element_line(size = 0.1, color = NA,
        linetype = "dashed"),
      panel.grid.major.x = element_line(size = 0.1, color = "gray70",
        linetype = "dashed"),
      plot.margin = unit(c(10, 5, 2, 0), "pt"),
      panel.spacing.x  =  unit(200, "pt"),
      panel.spacing.y = unit(20, "pt"),
      legend.position  =   "bottom",
      strip.text = element_text(margin = margin(t = 0.5, b = 0.5,
        l = 0.25, r = 0.25, unit = "lines"),
        size = 20,
        lineheight = 0.5,
        face = "bold",
        hjust = 0.9,
        vjust = -1)) +
    guides(fill  =  guide_legend(title.position = "top", title.hjust  =  0.5))
}



make_base_banner <- function(spatial_3Zone) {
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
  
  return(list(
    aus_sf = aus_sf,
    banner_thumb = banner_thumb,
    spatial_3Zone_sf = spatial_3Zone_sf
  ))
}

make_banner_thumbs <- function(banner_thumb, spatial_3Zone_sf, hues) {
  banner_thumb_all <-
    banner_thumb +
    geom_sf(data = spatial_3Zone_sf, fill = hues[4], color = NA) +
    geom_sf(data = spatial_3Zone_sf, fill = NA, color = "black", size = 0.2) +
    geom_sf(fill = "white", color = hues[4]) +
    geom_blank(aes(x = 250, y = -20))

  banner_thumb_northern <-
    banner_thumb +
    geom_sf(data = spatial_3Zone_sf[1, ], fill = hues[4], color = NA) +
    geom_sf(data = spatial_3Zone_sf[1, ], fill = NA, color = "black", size = 0.2) +
    geom_sf(fill = "white", color = hues[4]) +
    geom_blank(aes(x = 250, y = -20))

  banner_thumb_central <-
    banner_thumb +
    geom_sf(data = spatial_3Zone_sf[2, ], fill = hues[4], color = NA) +
    geom_sf(data = spatial_3Zone_sf[2, ], fill = NA, color = "black", size = 0.2) +
    geom_sf(fill = "white", color = hues[4]) +
    geom_blank(aes(x = 250, y = -20))

  banner_thumb_southern <-
    banner_thumb +
    geom_sf(data = spatial_3Zone_sf[3, ], fill = hues[4], color = NA) +
    geom_sf(data = spatial_3Zone_sf[3, ], fill = NA, color = "black", size = 0.2) +
    geom_sf(fill = "white", color = hues[4]) +
    geom_blank(aes(x = 250, y = -20))

  save(banner_thumb_all,
    banner_thumb_northern,
    banner_thumb_central,
    banner_thumb_southern,
    file = "../data/spatial/banner_thumbs.RData"
  )
  return(list(
    banner_thumb_all = banner_thumb_all,
    banner_thumb_northern = banner_thumb_northern,
    banner_thumb_central = banner_thumb_central,
    banner_thumb_southern = banner_thumb_southern
  ))
}

make_banner_thumbs_v <- function(aus_sf, banner_thumb, spatial_3Zone_sf, hues) {
  banner_thumb_v <-
    ggplot(aus_sf) +
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

  banner_thumb_northern_v <-
    banner_thumb_v +
    geom_sf(data = spatial_3Zone_sf[1, ], fill = hues[4], color = NA) +
    geom_sf(data = spatial_3Zone_sf[1, ], fill = NA, color = "black") +
    geom_sf(fill = "white", color = hues[4]) +
    theme(plot.margin = unit(c(2, 0, 100, 0), "pt"))

  banner_thumb_central_v <-
    banner_thumb_v +
    geom_sf(data = spatial_3Zone_sf[2, ], fill = hues[4], color = NA) +
    geom_sf(data = spatial_3Zone_sf[2, ], fill = NA, color = "black") +
    geom_sf(fill = "white", color = hues[4]) +
    coord_sf() +
    theme(plot.margin = unit(c(2, 0, 100, 0), "pt"))

  banner_thumb_southern_v <-
    banner_thumb_v +
    geom_sf(data = spatial_3Zone_sf[3, ], fill = hues[4], color = NA) +
    geom_sf(data = spatial_3Zone_sf[3, ], fill = NA, color = "black") +
    geom_sf(fill = "white", color = hues[4]) +
    coord_sf() +
    theme(plot.margin = unit(c(2, 0, 100, 0), "pt"))

  return(list(
    banner_thumb_v = banner_thumb_v,
    banner_thumb_northern_v = banner_thumb_northern_v,
    banner_thumb_central_v = banner_thumb_central_v,
    banner_thumb_southern_v = banner_thumb_southern_v
  ))
}

make_all_banners <- function() {
  load(file = "../data/spatial/spatial_3Zone.RData")
  banner_base_list <- make_base_banner(spatial_3Zone)
  banner_list <- make_banner_thumbs(
    banner_base_list$banner_thumb,
    banner_base_list$spatial_3Zone_sf,
    hues = hues
  )
  banner_list |> list2env(env = globalenv())
  banner_list_v <- make_banner_thumbs_v(
    banner_base_list$aus_sf,
    banner_base_list$banner_thumb,
    banner_base_list$spatial_3Zone_sf,
    hues = hues
  )
  banner_list_v |> list2env(env = globalenv())
}
