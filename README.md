
<!-- PLEASE DO NOT EDIT ./README.md BY HAND, EDIT ./inst/README.Rmd AND RENDER TO CREATE ./README.md -->

[![R-CMD-check](https://github.com/mixOmicsteam/mixOmics/workflows/R-CMD-check/badge.svg?event=pull_request)](https://github.com/mixOmicsteam/mixOmics/actions)
[![bioc release](https://img.shields.io/badge/bioc%20release-6.20.0-green.svg)](https://www.bioconductor.org/packages/mixOmics)
[![BioC status](http://www.bioconductor.org/shields/build/release/bioc/mixOmics.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/mixOmics)
[![codecov](https://codecov.io/github/mixOmicsTeam/mixOmics/graph/badge.svg?token=PzhCBOuU8o)](https://codecov.io/github/mixOmicsTeam/mixOmics)
[![download](http://www.bioconductor.org/shields/downloads/release/mixOmics.svg)](https://bioconductor.org/packages/stats/bioc/mixOmics)
[![last commit](https://img.shields.io/github/last-commit/mixOmicsTeam/mixOmics.svg)](https://github.com/mixOmicsTeam/mixOmics/commits/master)
[![license](https://img.shields.io/badge/license-GPL%20(%3E=%202)-lightgrey.svg)](https://choosealicense.com/)
[![dependencies](http://bioconductor.org/shields/dependencies/release/mixOmics.svg)](http://bioconductor.org/packages/release/bioc/html/mixOmics.html#since)

![](http://mixomics.org/wp-content/uploads/2019/07/MixOmics-Logo-1.png)

This repository contains the `R` package [now hosted on
Bioconductor](http://bioconductor.org/packages/release/bioc/html/mixOmics.html)
and our stable and development `GitHub` versions.

## Installation

(**macOS users only:** Ensure you have installed
[XQuartz](https://www.xquartz.org/) first.)

Make sure you have the latest R version and the latest `BiocManager`
package installed following [these
instructions](https://www.bioconductor.org/install/) (if you use legacy
R versions (\<=3.5.0) refer to the instructions at the end of the
mentioned page).

``` r
## install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

Ensure the following returns `TRUE`, or follow the guidelines provided
by the output.

``` r
BiocManager::valid()
```

For installation in R, see options a) and b). For Docker containers, see
c).

#### a) Latest `Bioconductor` Release

You can then install `mixOmics` using the following code:

``` r
## install mixOmics
BiocManager::install('mixOmics')
```

#### b) `GitHub` Versions

##### Stable version

Install the latest stable version (see below for latest
[development](https://github.com/ajabadi/mixOmics#development-version)
version) of `mixOmics` from `GitHub` (as bug-free as it can be):

``` r
BiocManager::install("mixOmicsTeam/mixOmics") 
```

Check after installation that the following code does not throw any
error (especially Mac users - refer to [installation
instructions](#installation)) and that the welcome message confirms you
have installed [the latest
version](https://github.com/mixOmicsTeam/mixOmics/blob/master/DESCRIPTION#L4):

``` r
library(mixOmics) 
#> Loaded mixOmics ?.?.?
```

##### Development version

You can also install the [development
version](https://github.com/mixOmicsTeam/mixOmics/blob/devel/DESCRIPTION#L4)
for new features yet to be widely tested (see [What’s
New](/https://github.com/ajabadi/mixOmics#whats-new)):

``` r
BiocManager::install("mixOmicsTeam/mixOmics@devel")
```

#### c) `Docker` container of the stable GitHub version

<details>
<summary>
Click to expand
</summary>

**Note: this requires root privileges**

1)  Install Docker following instructions at
    <https://docs.docker.com/docker-for-mac/install/>

**if your OS is not compatible with the latest version** download an
older version of Docker from the following link:

-   MacOS: <https://docs.docker.com/docker-for-mac/release-notes/>
-   Windows: <https://docs.docker.com/docker-for-windows/release-notes/>

Then open your system’s command line interface (e.g. Terminal for MacOS
and Command Promot for Windows) for the following steps.

**MacOS users only:** you will need to launch Docker Desktop to activate
your root privileges before running any docker commands from the command
line.

2)  Pull mixOmics container

``` bash
docker pull mixomicsteam/mixomics
```

3)  Ensure it is installed

The following command lists the running images:

``` bash
docker images
```

This lists the installed images. The output should be something similar
to the following:

    $ docker images 
      > REPOSITORY                       TAG       IMAGE ID       CREATED         SIZE
      > mixomicsteam/mixomics            latest    e755393ac247   2 weeks ago     4.38GB

4)  Active the container

Running the following command activates the container. You must change
`your_password` to a custom password of your own. You can also customise
ports (8787:8787) if desired/necessary. see
<https://docs.docker.com/config/containers/container-networking/> for
details.

``` bash
docker run -e PASSWORD=your_password --rm -p 8787:8787 mixomicsteam/mixomics
```

5)  Run

In your web browser, go to `http://localhost:8787/` (change port if
necessary) and login with the following credentials:

*username*: rstudio  
*password*: (your_password set in step 4)

6)  Inspect/stop

The following command lists the running containers:

``` bash
sudo docker ps
```

The output should be something similar to the following:

``` bash
$ sudo docker ps
  > CONTAINER ID   IMAGE                   COMMAND   CREATED         STATUS         PORTS                    NAMES
  > f14b0bc28326   mixomicsteam/mixomics   "/init"   7 minutes ago   Up 7 minutes   0.0.0.0:8787->8787/tcp   compassionate_mestorf
```

The listed image ID can then be used to stop the container (here
`f14b0bc28326`)

``` bash
docker stop f14b0bc28326
```

</details>

## Contribution

We welcome community contributions concordant with [our code of
conduct](https://github.com/mixOmicsTeam/mixOmics/blob/master/CODE_OF_CONDUCT.md).
We strongly recommend adhering to [Bioconductor’s coding
guide](https://bioconductor.org/developers/how-to/coding-style/) for
software consistency if you wish to contribute to `mixOmics` R codes.

### Bug reports and pull requests

To report a bug (or offer a solution for a bug!) visit:
<https://github.com/mixOmicsTeam/mixOmics/issues>. We fully welcome and
appreciate well-formatted and detailed pull requests. Preferably with
tests on our datasets.

<details>
<summary>
Set up development environment
</summary>

-   Install the latest version of R
-   Install RStudio
-   Clone this repo, checkout master branch, pull origin and then run:

``` r
install.packages("renv", Ncpus=4)
install.packages("devtools", Ncpus=4)

# restore the renv environment
renv::restore()

# or to initialise renv
# renv::init(bioconductor = TRUE)

# update the renv environment if needed
# renv::snapshot()

# test installation
devtools::install()
devtools::test()

# complete package check (takes a while)
devtools::check()
```

</details>

### Discussion forum

We wish to make our discussions transparent so please direct your
analysis questions to our discussion forum
<https://mixomics-users.discourse.group>. This forum is aimed to host
discussions on choices of multivariate analyses, as well as comments and
suggestions to improve the package. We hope to create an active
community of users, data analysts, developers and R programmers alike!
Thank you!

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
you wish to be involved!**.

The project started at the *Institut de Mathématiques de Toulouse* in
France, and has been fully implemented in Australia, at the *University
of Queensland*, Brisbane (2009 – 2016) and at the *University of
Melbourne*, Australia (from 2017). We focus on the development of
computational and statistical methods for biological data integration
and their implementation in `mixOmics`.

## Why this toolkit?

`mixOmics` offers a wide range of novel multivariate methods for the
exploration and integration of biological datasets with a particular
focus on variable selection. Single ’omics analysis does not provide
enough information to give a deep understanding of a biological system,
but we can obtain a more holistic view of a system by combining multiple
’omics analyses. Our `mixOmics` R package proposes a whole range of
multivariate methods that we developed and validated on many biological
studies to gain more insight into ’omics biological studies.

## Want to know more?

www.mixOmics.org (tutorials and resources)

Our latest bookdown vignette:
<https://mixomicsteam.github.io/Bookdown/>.

## Different types of methods

We have developed 17 novel multivariate methods (the package includes 19
methods in total). The names are full of acronyms, but are represented
in this diagram. *PLS* stands for *Projection to Latent Structures*
(also called Partial Least Squares, but not our preferred nomenclature),
*CCA* for *Canonical Correlation Analysis*.

That’s it! Ready! Set! Go!

Thank you for using `mixOmics`!

![](http://mixomics.org/wp-content/uploads/2012/04/framework-mixOmics-June2016.jpg)

## What’s New

#### March 2022

-   bug fix implemented for [Issue
    \#196](https://github.com/mixOmicsTeam/mixOmics/issues/196).
    `perf()` can now handle features with a `(s)pls` which have near
    zero variance.
-   bug fix implemented for [Issue
    \#192](https://github.com/mixOmicsTeam/mixOmics/issues/192).
    `predict()` can now handle when the testing and training data have
    their columns in different orders.
-   bug fix implemented for [Issue
    \#178](https://github.com/mixOmicsTeam/mixOmics/issues/178). If the
    `indY` parameter is used in `block.spls()`, `circosPlot()` can now
    properly identify the
    ![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
    dataframe.
-   bug fix implemented for [Issue
    \#172](https://github.com/mixOmicsTeam/mixOmics/issues/172).
    `perf()` now returns values for the `choice.ncomp` component when
    `nrepeat`
    ![\< 3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%3C%203 "< 3")
    whereas before it would just return `NA`s.
-   bug fix implemented for [Issue
    \#171](https://github.com/mixOmicsTeam/mixOmics/issues/171). `cim()`
    now can take `pca` objects as input.
-   bug fix implemented for [Issue
    \#161](https://github.com/mixOmicsTeam/mixOmics/issues/161).
    `tune.spca()` can now handle `NA` values appropriately.
-   bug fix implemented for [Issue
    \#150](https://github.com/mixOmicsTeam/mixOmics/issues/150).
    Provided users with a specific error message for when `plotArrow()`
    is run on a `(mint).(s)plsda` object.
-   bug fix implemented for [Issue
    \#122](https://github.com/mixOmicsTeam/mixOmics/issues/122).
    Provided users with a specific error message for when a `splsda`
    object that has only one sample associated with a given class is
    passed to `perf()`.
-   bug fix implemented for [Issue
    \#120](https://github.com/mixOmicsTeam/mixOmics/issues/120).
    `plotLoadings()` now returns the loading values for features from
    **all** dataframes rather than just the last one when operating on a
    `(mint).(block).(s)plsda` object.
-   bug fix implemented for [Issue
    \#43](https://github.com/mixOmicsTeam/mixOmics/issues/43).
    Homogenised the way in which `tune.mint.splsda()` and
    `perf.mint.splsda()` calculate balanced error rate (BER) as there
    was disparity between them. Also made the global BER a weighted
    average of BERs across each study.
-   enhancement implemented for [Issue
    \#30/#34](https://github.com/mixOmicsTeam/mixOmics/issues/34). The
    parameter `verbose.call` was added to most of the methods. This
    parameter allows users to access the specific values input into the
    call of a function from its output.
-   bug fix implemented for [Issue
    \#24](https://github.com/mixOmicsTeam/mixOmics/issues/24).
    `background.predict()` can now operate on `mint.splsda` objects and
    can be used as part of `plotIndiv()`.

#### July 2021

-   new function `plotMarkers` to visualise the selected features in
    block analyses (see
    <https://github.com/mixOmicsTeam/mixOmics/issues/134>)
-   `tune.spls` now able to tune the selected variables on both `X` and
    `Y`. See `?tune.spls`
-   new function `impute.nipals` to impute missing values using the
    nipals algorithm
-   new function `tune.spca` to tune the number of selected variables
    for pca components
-   `circosPlot` now has methods for `block.spls` objects. It can now
    handle similar feature names across blocks. It is also much more
    customisable. See advanced arguments in `?circosPlot`
-   new `biplot` function for `pca` and `pls` objects. See
    `?mixOmics::biplot`
-   `plotDiablo` now takes `col.per.group` (see \#119)

#### April 2020

-   weighted consensus plots for DIABLO objects now consider
    per-component weights

#### March 2020

-   `plotIndiv` now supports (weighted) consensus plots for block
    analyses. See the example in [this
    issue](https://github.com/mixOmicsTeam/mixOmics/issues/57)
-   `plotIndiv(..., ind.names=FALSE)` [warning
    issue](https://github.com/mixOmicsTeam/mixOmics/issues/59) now fixed

#### January 2020

-   `perf.block.splsda` now supports calculation of combined AUC
-   `block.splsda` bug which could drop some classes with
    `near.zero.variance=TRUE` now fixed
