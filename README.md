# SyncER - Synchronicity testing for Event Records

#### Usage

*SyncER* is a toolbox developed to **sync**hronize age-depth models and evaluate the potential synchronous deposition of correlated **event** deposits between geological **records**.

This notebook presents the *SyncER* workflow by illustrating the different steps on a synthetic test dataset. You can adjust the executable cells to apply it to your own dataset, but it is recommendable to instead copy the different steps to your personal R script for execution.

When using *SyncER*, cite the following paper:

> Wils, K, Ramisch, A. (2026) SyncER: Synchronicity testing for Event Records using tied age-depth models. *Scientific Reports.*

If you make use of the built-in compatibility with *rbacon* for age-depth modelling, you should also cite the original reference for this package:

> Blaauw M, Christen JA. (2011) Flexible paleoclimate age-depth models using an autoregressive gamma process. *Bayesian Analysis* 6 (3), 457-474. DOI: 10.1214/11-ba618.

If you make use of the built-in compatibility with *rplum* for age-depth modelling, you should cite the above reference for *rbacon* as well as the one below:

> Aquino-López, M. A., Blaauw, M., Christen, J.A., and Sanderson, N. K. (2018) Bayesian analysis of 210Pb dating. *Journal of Agricultural, Biological and Environmental Statistics* 23 (3), 317-333. DOI: 10.1007/s13253-018-0328-7

To use the *SyncER* package, download it from [GitHub](https://github.com/katleenwils/SyncER) and load it into your environment. The example dataset used throughout the vignette is distributed with the package, so no separate download is needed.

```{r SyncER-install}
install.packages("remotes")
remotes::install_github("katleenwils/SyncER")
library(SyncER)
```

Once *SyncER* is available on CRAN, you can use the following commands:

```{r SyncER-install}
install.packages("SyncER")
library(SyncER)
```

#### Reporting bugs

If you run into a problem or unexpected behaviour, please open an issue on the [GitHub issue tracker](https://github.com/katleenwils/SyncER/issues) or contact the developer. Questions and suggestions are welcome there too, and will be accommodated as much as reasonably possible. A report is most useful when it includes:

- the version of *SyncER* you are running (`packageVersion("SyncER")`),
- the output of `sessionInfo()`,
- a minimal example that reproduces the problem — the bundled example dataset is ideal,
- what you expected to happen, and what happened instead.

#### Disclaimer

*SyncER* is research software and is provided **as is**, without warranty of any kind, either expressed or implied. The entire risk as to the quality and performance of the software is with you. Interpretations drawn from *SyncER* outputs remain the responsibility of the user. Synchronicity scores are conditional on the quality of the underlying age-depth models, and synchronized models gain **precision** without necessarily gaining **accuracy** — relative synchronicity is only ever evaluated with respect to the other records included.

SyncER © 2026 Katleen Wils. Licensed under the GNU General Public License v3.0 (GPL-3); see the LICENSE.md file for the full text.
