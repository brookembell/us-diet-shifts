# Repository for U.S. Diet Shifts Paper

This Github repository contains code and data for the manuscript titled, "Shifting to recommended diets in the US improves health outcomes, but has differing environmental, social, and economic impacts," currently under review.

This repository is also an R project, and the `us-diet-shifts.Rproj` file must be opened in RStudio first in order for the code to run properly.

## System Requirements

R code can be run on MacOS, Windows, or Linux. See more here: https://cran.rstudio.com/.

This code was run using R (version 4.4.0) and RStudio Desktop (version 2024.04.1) on MacOS 15.6.1. 

RStudio requires a 64-bit operating system and R version 3.6.0+.

## R Packages

The following R packages will need to be installed prior to running the code.

```
install.packages(c("tidyverse", "data.table", "foreach", "doParallel", "magrittr"))
```

Installing these packages will typically take a few minutes.

## Organization

There are three main folders:

- code: contains the R code needed to run all analyses that were conducted for the manuscript. the main script is title `run_models.RMD`. other scripts and functions that are sourced are located in the "functions" and "source_code" folders.
- data: contains all datasets needed to run the analyses. 
- figures: contains the model output.
  * CRA: contains the output from the Comparative Risk Assessment (CRA).
  * envecosoc: contains the output from the environmental, economic, and social ("envecosoc") models
  
## Running the Code

When running all of the options (e.g., all food groups, all diet patterns, all 1,000 simulations, etc.), this code takes a couple of days to run (on a high performance computer).

Therefore, for testing and/or review purposes, it is *highly* recommended to only run a subset of the options (e.g., 10 simulations, only 2 diet patterns, etc.).
