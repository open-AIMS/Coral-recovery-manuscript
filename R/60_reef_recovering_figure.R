source("helper_functions.R")
cr_check_packages()


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
## More data processing
{
  prop.recovered.bench <-
     manual |>
     dplyr::filter(recovered.bench != "NA") |>
     group_by(Region, recovery.year) |>
     dplyr::count(recovered.bench) |>
     mutate(recovered.bench = factor(recovered.bench)) |>
     mutate(
       data4plot = ifelse(recovered.bench == "no", n * -1, n),
       Region = factor(Region,
         levels = c("Northern", "Central", "Southern")
       )
     ) |>
     ungroup()

   prop.recovered.prior <-
     manual |>
     dplyr::filter(recovered.prior != "NA") |>
     group_by(Region, recovery.year) |>
     dplyr::count(recovered.prior) |>
     mutate(recovered.prior = factor(recovered.prior)) |>
     mutate(data4plot = ifelse(recovered.prior == "no", n * -1, n)) |>
     ungroup()


   prop.all.bench <-
     manual |>
     dplyr::filter(recovered.bench != "NA") |>
     dplyr::count(recovered.bench)

   prop.all.prior <-
     manual |>
     dplyr::filter(recovered.bench != "NA") |>
     dplyr::count(recovered.prior)

   Years <- expand.grid(
     recovery.year = seq(from = 1985, to = 2022, by = 1),
     Region = levels(prop.recovered.bench$Region)
   ) |>
     mutate(recovery.year = factor(recovery.year))
   Years |>
     pull(recovery.year) |>
     levels()
}
## Plot a
{
    prop.recovered.prior |>
    filter(!is.na(Region)) |> 
    ggplot(aes(x = recovery.year, y = data4plot,
      fill = recovered.prior)) +
    geom_bar(stat = "identity", colour = "black",
      linewidth = 0.05) +
    scale_x_continuous("Year of recovery",
      breaks = c(1990, 2000, 2010, 2020)) +
    scale_y_continuous("Number of recoveries",
      breaks = c(-20, -15, -10, -5, 0, 5, 10, 15, 20),
      labels = c(20, 15, 10, 5, 0, 5, 10, 15, 20), limits = c(-17, 17)) +
    scale_fill_manual("", values = c("grey", "white"),
      breaks = c("no", "yes"), labels = c("Partial recovery", "Full recovery")) +
    geom_hline(yintercept = 0) +
    ggtitle("a) Recovery to prior benchmark") +
    facet_grid(~Region) +
    theme_classic() +
    theme(
      legend.position = c(0.175, 0.9),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 8),
      legend.key.width = unit(0.25, "cm"),
      legend.key.height = unit(0.25, "cm"),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 7)
    )

  prop.prior.plot
}
## Plot b
{
  prop.bench.plot <-
    prop.recovered.bench %>%
    filter(!is.na(Region)) |> 
    ggplot(
      aes(x = recovery.year, y = data4plot, fill = recovered.bench)
    ) +
    geom_bar(stat = "identity", colour = "black",
      show.legend = F, linewidth = 0.05) +
    scale_x_continuous("Year of recovery",
      breaks = c(1990, 2000, 2010, 2020)) +
    scale_y_continuous("Number of recoveries",
      breaks = c(-20, -15, -10, -5, 0, 5, 10, 15, 20), labels = c(20, 15, 10, 5, 0, 5, 10, 15, 20), limits = c(-18, 17)) +
    scale_fill_manual("", values = c("grey", "white"),
      breaks = c("no", "yes"), labels = c("Partial recovery", "Full recovery")) +
    geom_hline(yintercept = 0) +
    ggtitle("b) Recovery to historical benchmark") +
    facet_grid(~Region) +
    theme_classic() +
    theme(
      legend.position = c(0.9, 0.8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 7)
    )

  prop.bench.plot
}
## Plot c
{
  prior.recovery.hist <-
    manual |> ggplot(aes(x = rel.perc.prior.recovery)) +
    geom_histogram(fill = "grey", colour = "black",
      binwidth = 10, boundary = 10, linewidth = 0.05) +
    geom_vline(xintercept = 100, linetype = "dashed",
      colour = "red", linewidth = 0.2) +
    geom_vline(xintercept = 50, linetype = "dashed",
      colour = "blue", linewidth = 0.2) +
    ggtitle("c)") +
    scale_x_continuous("Relative percent recovery") +
    scale_y_continuous("Count") +
    facet_wrap(~Region) + # ,scales = 'free_x')+
    theme_classic() +
    theme(axis.text = element_text(size = 8))

  prior.recovery.hist
}
## Plot d
{
  bench.recovery.hist <-
    manual |>
    ggplot(aes(x = rel.perc.benchmark.recovery)) +
    geom_histogram(fill = "grey", colour = "black",
      binwidth = 10, boundary = 10, linewidth = 0.05) +
    geom_vline(xintercept = 100, linetype = "dashed",
      colour = "red", linewidth = 0.2) +
    geom_vline(xintercept = 50, linetype = "dashed",
      colour = "blue", linewidth = 0.2) +
    ggtitle("d)") +
    scale_x_continuous("Relative percent recovery") +
    scale_y_continuous("Count") +
    facet_wrap(~Region) + # ,scales = 'free_x')+
    theme_classic() +
    theme(axis.text = element_text(size = 8))

  bench.recovery.hist
}
## Put plots together
{
  layout <-
    "AABB
CCDD"

  recovery.plot.multipanel <-
    prop.prior.plot +
    prop.bench.plot +
    prior.recovery.hist +
    bench.recovery.hist +
    plot_layout(design = layout)

  recovery.plot.multipanel

  ggsave(
    filename = "../outputs/figures/figure_6.pdf",
    recovery.plot.multipanel,
    height = 20, width = 18, units = "cm"
  )

  ggsave(
    filename = "../outputs/figures/figure_6.png",
    recovery.plot.multipanel,
    height = 20, width = 18, units = "cm", dpi = 600
  )
}







