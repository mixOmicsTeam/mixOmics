## version 6.18.0

### new features / enhancements / changes

* new function `plotMarkers` to visualise the selected features in block analyses (see #134)
* `auroc` title now fixed (#135)
* `cimDiablo` takes `trim` argument to customise outlier filtering (#136)
* `plotIndiv.pca` default shape set to `16`
* `circosPlot` & `network` now support blocks with similar feature names 
* `circosPlot` now has methods for block.spls objects
* `circosPlot` now takes new formal and advanced args for customisation. See `?circosPlot`.

### bug fixes

* `plotVar` legend colour mismatch bug fixed
* `plotDiablo` error `undefined variable (Y)` fixed
* `nipals` initialisation bug with high-variance high-NA rate column fixed 
* `cim` bug with high NA rate data fixed using imputation by the column mean

## version 6.16.0

### new features / enhancements / changes

* `biplot` default colours now match `plotIndiv`
* `biplot` now takes `cex` argument
* `nipals` now takes `center` and `scale` arguments
* `nipals` now only outputs `p`, `t` and `eig`
* new function `impute.nipals` to impute missing values using NIPALS algorithm
* `nipals` function checks for orthogonality of components when NA present
* `plot.tune` legends improved
* `plot.tune` now uses colour-blind friendly colours
* new `tune.spls` function to perform variable selection on both X and Y
* `tune.spls` now chooses optimal keepX even if nrepeat < 3
* `tune.spls` now takes `validation` arg
* `tune.spca` is now much faster
* `circosPlot` links can now have adjustable width (#118)
* `plotDiablo` now takes `col.per.group` (#119)
* terminology change: `consensus` renamed to `average` in Diablo context

### bug fixes

* `plotVar` bug fixed

## version 6.14.0

### new features / enhancements / changes

* `circosPlot`: The radial location of feature names can now be cutomised using `var.adj`
* added `plot` and `print` methods for `nipals` ouput (#87)
* all Discriminat Analyses now run solely on `mode=regression` (#79)
* `cim` argument change: `threshold` replaced by `cutoff`
* `nipals` and `pca` with missing values allow skipping reconstitution of the input matrix
* `tune.block.splsda` now allows random number seed also for parallel processing (#72)
* New `biplot` methods for the pca family (#90)

### bug fixes

* `plotIndiv`: Legend bug which misspecified the groups resolved
* `plotIndiv`: Legends now ordered as inputted, and not alphabetically
* `plot` method issue for `spca` resolved 
* `plotLoadings.spca` bug with `var.names` now fixed (#81)
* `ipca` deprecation warning fixed

-------------------------------------------------------------------------------
## version 6.12.0

### new features / enhancements

* `plotLoadings`'s infamous *figure margins too large* error now handled and informative condition thrown
* `circosPlot`'s `lines` argument default to `FALSE` now
* `circosPlot`'s inconsistentcy of blocks with identical `X` names fixed
*  consensus and weighted consensus plots now supported for `plotIndiv` with relevant block analyses
* `plotLoadings`'s feature name trimming can be customised
* `block.splsda` bug which could drop some Y factors with `near.zero.variance=TRUE` fixed
* `perf.block.splsda` now supports calculation of combined and per-block AUC
* model improvement significance can be custmoised in all `perf` and `tune` functions
* `perf.block.splsda` is now much faster and supports FORK clusters
* `tune.(s)pls(da)`, `perf.(s)plsda` now support FORK clusters

### bug fixes

* `circosPlot`'s faded `lines` bug when many `NA`s present fixed
* `tune.block.splsda()` bug when using fixed `test.keepX` over two or more blocks fixed
* `circosPlot` and `plotLoadings` bug caused by features with `NA`s fixed
* `plotIndiv(..., ind.names = FALSE)` warning/bug fixed
* `tune.block.splsda` bug on Windows parallelisation fixed
* `perf` and `tune` functions' issue  when choosing the optimum component resolved
* added option to suppress `auroc` from printing all the AUCs
-------------------------------------------------------------------------------
## version 6.10.0

### new features / enhancements

* parallel processing on `tune.block.splsda` improved
* `tune.block.splsda` now supports more distances
* You can now customise `auroc` plots. Refer to documentation for more info

### bug fixes

* single factor multilevel error in `pls` fixed
* fixed over-estimated correlation of `cim` for `mixo_(s)pls` objects with single component 
* margin error in `cim` now handled properly
* fixed `plotLoadings` error for very long variable names
* `predict` function bug for single sample prediction fixed
* `plotLoadings` bug for long variable names fixed
* Fixed `tune.spls` and `pef.plsda` bugs when using `cpus` argument for parallel 
processing
* `perf.plot` bug in extracting names fixed
* Few fixes for `tune.splsda` with AUC

### minor improvements

* missing values in `plotIndiv`'s `group` argument no more throws error
* `mixOmics::predict` function documentation now more accessible
* names of `linnerud` datasets fixed
* `plot.perf` now respects `ylim` arguments for custom y range
* package startup message with direct liks to useful resources
* `mixOmics` function documentation disambiguated with instruction on how to get
package help
* Updated onLoad message with discussion forum info, bug reports, and more
* Dropped legacy `comp.tol` argument from `pca`
* `plot.perf` now respects `ylim` arguments for custom y range
* Added Code of Conduct

-------------------------------------------------------------------------------
## version 6.8.0

* NOW HOSTED ON BIOCONDUCTOR

-------------------------------------------------------------------------------

## version < 6.8.0

* Refer to `./inst/legacy/NEWS-old` on GitHub repo