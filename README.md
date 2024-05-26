Readme
=============

Code repository associated with analyses and Figures in Emslie et al
(2024) "Increasing disturbance frequency undermines coral recovery"


```
/root
|- data
|   |- processed
|   |- modelled
|- R
   |- 00_main.R
   |- 10_get_data.R
|- outputs
   |- figures.R

```



# Analysis steps

Step through each of the lines of `00_main.R` (first ensuring that the
`R` folder is the working directory):

1. Source all the additional helper functions required for the various
   processing analysis routines

2. Run the function `cr_check_packages()` that `requires` each of the
   necessary packages. If any are missing, you should install them
   (and rerun this step) before proceeding.

3. Run the function `cr_paths()` to ensure all the necessary paths are
   in place.

4. *Optional* - examine script `10_get_data.R` 

   NOTE: this first processing/analysis script (`10_get_data.R`) is
   provided so as to illustrate the importation and broad processing
   steps taken in preparing the data for modelling. Access to the
   primary database and datasets is not provided within this
   repository. Access to the primary data can be requested by emailing
   the senior author directly.

5. Run the `20_disturbance_all_figure.R` script in order to generate
   Figure 2.

   This script takes the following inputs:

   - `../data/spatial/spatial_3Zone.RData`
   - `../data/modelled/cots.sum.all_3Zone.RData`
   - `../data/modelled/bleaching.sum.all_3Zone.RData`
   - `../data/modelled/cyclones.sum.all_3Zone.RData`

   This script produces the following outputs:

   - `../outputs/figures/Disturbances_all_no_label_transposed.pdf`

6. Run the `30_disturbance_severe_figure.R` script in order to
   generate Figure 3.
   
   This script takes the following inputs:

   - `../data/modelled/bleaching.full_3Zone.RData`
   - `../data/modelled/cots.full_3Zone.RData`
   - `../data/modelled/cyclones.full_3Zone.RData`
   - `../data/processed/all.reefs.cyclones.RData`
   - `../data/modelled/cots.sum.all_3Zone.RData`

   This script produces the following outputs:

   - `../outputs/figures/Disturbances_severe_compilation_newA_2021.pdf`
   
7. Explore the `40_coral_trends_models.R` script in order to generate
   the models underlying Figure 4.
   
   NOTE: this script is provided so as to illustrate the steps taken
   to model coral cover. Access to the input data is not provided
   within this repository. Access to the primary data can be requested
   by emailing the senior author directly.

8. Run the `45_coral_trends_figure.R` script in order to
   generate Figure 4.

   This script takes the following inputs:

   - `../data/modelled/mod.northern_brms.beta.ry.disp.RData`
   - `../data/modelled/mod.central_brms.beta.ry.disp.RData`
   - `../data/modelled/mod.southern_brms.beta.ry.disp.RData`

   This script produces the following outputs:

   - `../outputs/figures/threePanels_brms.beta.ry.disp_.pdf`
   

