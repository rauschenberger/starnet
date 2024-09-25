
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/rauschenberger/starnet?svg=true)](https://ci.appveyor.com/project/rauschenberger/starnet)
[![R-CMD-check](https://github.com/rauschenberger/starnet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rauschenberger/starnet/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/rauschenberger/starnet/graph/badge.svg)](https://app.codecov.io/gh/rauschenberger/starnet)

# Stacked Elastic Net Regression

Implements stacked elastic net regression (combining multiple mixing parameters rather than selecting one mixing parameter).

## Installation

Install the current release from
[CRAN](https://CRAN.R-project.org/package=starnet):

``` r
install.packages("starnet")
```

or the latest development version from
[GitHub](https://github.com/rauschenberger/starnet):

``` r
#install.packages("remotes")
remotes::install_github("rauschenberger/starnet")
```

## Reference

Armin Rauschenberger
[![AR](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0001-6498-4801),
Enrico Glaab
[![EG](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-3977-7469),
and Mark A. van de Wiel
[![MvdW](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-4780-8472)
(2021).
"Predictive and interpretable models via the stacked elastic net".
*Bioinformatics*
37(14):2012-2016.
[doi: 10.1093/bioinformatics/btaa535](https://doi.org/10.1093/bioinformatics/btaa535).

[![CRAN version](https://www.r-pkg.org/badges/version/starnet)](https://CRAN.R-project.org/package=starnet)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/starnet)](https://CRAN.R-project.org/package=starnet)
[![Total CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/starnet)](https://CRAN.R-project.org/package=starnet)

## Disclaimer

The R package `starnet` implements stacked elastic net regression ([Rauschenberger et al., 2021](https://doi.org/10.1093/bioinformatics/btaa535)).

Copyright &copy; 2019 Armin Rauschenberger, University of Luxembourg, Luxembourg Centre for Systems Biomedicine (LCSB), Biomedical Data Science (BDS)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
