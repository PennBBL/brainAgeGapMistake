# Statistical Pitfalls in Brain Age Analyses

This repository contains all of the scripts for the manuscript, Statistical Pitfalls in Brain Age Analyses. It also contains scripts to produce a figure for a poster at OHBM, and a couple of scripts to help illustrate important concepts relevant to the paper to peers.

`equation8.R` creates Figure 1.

`realDataExample.R` creates Figure 2.

`followalong.R` creates a worksheet to demonstrate the nature of residuals and the proposed modification under a model where the brain feature does not covary with age at all.

`mbagSimulation.R` demonstrates that the reason MBAG is not a vector of residuals is because of regression on residuals.

`posterFigure.R` creates a gganatogram for Ellyn's OHBM poster.

`regressionToTheMeanSimulation.R` illustrates that regression to the mean results in under prediction of subjects above the mean on the outcome variable, on average, and over prediction of subjects below the mean on the outcome variable, on average.
