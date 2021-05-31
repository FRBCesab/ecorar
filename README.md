# ecorar

Research compendium to reproduce analyses and figures of the following article:

> Loiseau N, Mouquet N, Casajus N, Grenié M, Guéguen M, Maitner B, Mouillot D, Ostling A, Renaud J, Tucker C, Velez L, Thuiller W & Violle C (2020) Global distribution and conservation status of ecologically rare mammal and bird species. _Nature Communications_, **11**, 5071. DOI: [10.1038/s41467-020-18779-w](http://dx.doi.org/10.1038/s41467-020-18779-w).

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

Cheers!
