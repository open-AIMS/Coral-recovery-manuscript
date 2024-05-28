source("helper_functions.R")
cr_check_packages()

cr_paths()

## NOTE: the first script (10_get_data.R) is provided so as to
## illustrate the importation and broad processing steps taken in
## preparing the data for modelling. Access to the primary database
## and datasets is not provided within this repository. Access to the
## primary data can be requested by emailing the senior author
## directly.
source("10_get_data.R")

## Figure 2
source("20_disturbance_all_figure.R")

## Figure 3
source("30_disturbance_severe_figure.R")

## Figure 4
## NOTE: this script is provided so as to illustrate the steps taken
## to model coral cover. Access to the input data is not provided
## within this repository. Access to the primary data can be requested
## by emailing the senior author directly.
source("40_coral_trends_modells.R")
source("45_coral_trends_figure.R")

## Figure 5
source("50_disturbance_year_figure.R")

## Figure 6
source("60_reef_recovery_figure.R")

## Figure 7
source("70_recovery_rate_figure.R")

## Figure 8
source("80_drivers_pf_recovery_figure.R")
