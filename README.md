# SyncER - Synchronicity testing for Event Records

*SyncER* is a toolbox developed to **sync**hronize age-depth models and evaluate the potential synchronous deposition of correlated **event** deposits between geological **records**.

This notebook presents the *SyncER* workflow by illustrating the different steps on a synthetic test dataset. You can adjust the executable cells to apply it to your own dataset, but it is recommendable to instead copy the different steps to your personal R script for execution.

When using *SyncER*, please cite the following paper:

> Wils, K, Ramisch, A. SyncER: Synchronicity testing for Event Records using tied age-depth models. *Scientific Reports.*

If you make use of the built-in compatibility with *rbacon* for age-depth modelling, you should also cite the original reference for this package:

> Blaauw M, Christen JA. (2011) Flexible paleoclimate age-depth models using an autoregressive gamma process. *Bayesian Analysis* 6 (3), 457-474. DOI: 10.1214/11-ba618.

If you make use of the built-in compatibility with *rplum* for age-depth modelling, you should cite the above reference for *rbacon* as well as the one below:

> Aquino-López, M. A., Blaauw, M., Christen, J.A., and Sanderson, N. K. (2018) Bayesian analysis of 210Pb dating. *Journal of Agricultural, Biological and Environmental Statistics* 23 (3), 317-333. DOI: 10.1007/s13253-018-0328-7

To use the *SyncER* package, download it from [GitHub](https://github.com/katleenwils/SyncER) and load it into your environment. The example dataset used throughout the vignette is distributed with the package, so no separate download is needed.

```{r SyncER-install}
install.packages("remotes")
remotes::install_github("katleenwils/SyncER", build_vignettes = TRUE)
library(SyncER)
```
