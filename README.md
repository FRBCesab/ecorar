# ecorar

R code and data to reproduce figures of Loiseau, Mouquet et al.'s 2020 article.

<hr />



## General

This repository is structured as follow:

- `data/`: contains all data required to reproduce figures
- `R/`: contains R functions developed for this project
- `man/`: contains R functions documentation
- `analyses/`: contains R scripts (one per figure) and a setup file defining parameters



## Notes

- All required packages will be installed (if necessary) and loaded.
- Figures will be stored in `figures/`



## Usage

Clone the repository and run this command in R/RStudio:

```r
source("make.R")
```

You can customize figures by editing `analyses/setup.R`.

Cheers!
