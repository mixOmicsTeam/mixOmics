
![](http://mixomics.org/wp-content/uploads/2019/07/MixOmics-Logo-1.png)

This repository contains the `R` package [now hosted on
Bioconductor](http://bioconductor.org/packages/release/bioc/html/mixOmics.html)
and our current `GitHub` version.

## Installation

**(Mac Users Only:)** Ensure you have installed
[XQuartz](https://www.xquartz.org/) first.

#### Latest Bioconductor Release

Make sure you have the latest R version and the latest `BiocManager`
package installed following [these
instructions](https://www.bioconductor.org/install/) (if you use legacy
R versions (\<=3.5.0) refer to the instructions at the end of the
mentioned page), you can then install `mixOmics` using the following
code:

``` r
## install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
## install mixOmics
BiocManager::install('mixOmics')
```

#### Latest `GitHub` Version

Install the [devtools](https://github.com/r-lib/devtools) package in R,
then load it and install the latest stable version of `mixOmics` from
`GitHub` (as bug-free as it can be):

``` r
## install devtools if not installed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
## install mixOmics
devtools::install_github("mixOmicsTeam/mixOmics")
```

You can also install the development version:

``` r
devtools::install_github("mixOmicsTeam/mixOmics", ref="devel")
```

Check after installation that the following code does not throw any
error (especially Mac users - refer to [installation
instructions](#installation)) and that the welcome message confirms you
have installed latest version as in the latest package [DESCRIPTION
file](https://github.com/mixOmicsTeam/mixOmics/blob/master/DESCRIPTION#L4):

``` r
library(mixOmics) 
#> Loaded mixOmics ?.?.?
```

## Contribution

We welcome community contributions concordant with [our code of
conduct](https://github.com/mixOmicsTeam/mixOmics/blob/master/CODE_OF_CONDUCT.md).
We strongly recommend adhering to [Bioconductor’s coding
guide](https://bioconductor.org/developers/how-to/coding-style/) for
software consistenncy.

### Bug reports and pull requests

To report a bug (or offer a solution for a bug\!):
<https://github.com/mixOmicsTeam/mixOmics/issues>. We fully welcome and
appreciate well-formatted and detailed pull requests. Preferrably with
tests on our datasets.

### Discussion forum

We wish to make our discussions transparent so please direct your
questions to our discussion forum
<https://mixomics-users.discourse.group>. This forum is aimed to host
discussions on choices of multivariate analyses, bug report as well as
comments and suggestions to improve the package. We hope to create an
active community of users, data analysts, developers and R programmers
alike\! Thank you\!

## About the `mixOmics` team

`mixOmics` is collaborative project between Australia (Melbourne),
France (Toulouse), and Canada (Vancouver). The core team includes
Kim-Anh Lê Cao - <https://lecao-lab.science.unimelb.edu.au> (University
of Melbourne), Florian Rohart - <http://florian.rohart.free.fr>
(Toulouse) and Sébastien Déjean -
<https://perso.math.univ-toulouse.fr/dejean/>. We also have key
contributors, past (Benoît Gautier, François Bartolo) and present (Al
Abadi, University of Melbourne) and several collaborators including
Amrit Singh (University of British Columbia), Olivier Chapleur (IRSTEA,
Paris), Antoine Bodein (Universite de Laval) - **it could be you too, if
you wish to be involved\!**.

The project started at the *Institut de Mathématiques de Toulouse* in
France, and has been fully implemented in Australia, at the *University
of Queensland*, Brisbane (2009 – 2016) and at the *University of
Melbourne*, Australia (from 2017). We focus on the development of
computational and statistical methods for biological data integration
and their implementation in `mixOmics`.

## Why this toolkit?

`mixOmics` offers a wide range of novel multivariate methods for the
exploration and integration of biological datasets with a particular
focus on variable selection. Single ‘omics analysis does not provide
enough information to give a deep understanding of a biological system,
but we can obtain a more holistic view of a system by combining multiple
‘omics analyses. Our `mixOmics` R package proposes a whole range of
multivariate methods that we developed and validated on many biological
studies to gain more insight into ‘omics biological studies.

## Want to know more?

www.mixOmics.org (tutorials and resources)

Our latest bookdown vignette:
<https://mixomicsteam.github.io/Bookdown/>.

## Different types of methods

We have developed 17 novel multivariate methods (the package includes 19
methods in total). The names are full of acronyms, but are represented
in this diagram. *PLS* stands for *Projection to Latent Structures*
(also called Partial Least Squares, but not our prefered nomenclature),
*CCA* for *Canonical Correlation Analysis*.

That’s it\! Ready\! Set\! Go\!

Thank you for using
`mixOmics`\!

![](http://mixomics.org/wp-content/uploads/2012/04/framework-mixOmics-June2016.jpg)

## What’s New

#### November 2019

  - Parallel computing improved for `tune` and `perf` functions, and
    even more on Unix-like systems.

  - Fixed margin error problem with `plotLoadings`. See the example in
    [this issue](https://github.com/mixOmicsTeam/mixOmics/issues/45).

  - `cim` bug which overestimated correlations for single component now
    fixed.

  - `perf.sgccda` now supports calculation of average combined AUC.

#### September 2019

  - You can now customise `auroc` plots in version 6.8.3. See example
    [here](https://github.com/mixOmicsTeam/mixOmics/issues/35).

#### August 2019

  - The infamous `plot.perf` bug has been fixed in version 6.8.2. See
    example [here](https://github.com/mixOmicsTeam/mixOmics/issues/27).
