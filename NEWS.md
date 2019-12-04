
## version 6.12.0

### new features / enhancements

* `perf.block.splsda` now supports calculation of combined and per-block AUC
* improvement significance can be custmoised in all `perf` and `tune` functions
* `perf.block.splsda` is now much faster and supports FORK clusters
* added option to suppress `auroc` from printing all the AUCs
* `tune.(s)pls(da)`, `perf.(s)plsda` now support FORK clusters

### bug fixes

* `tune.block.splsda` bug on Windows parallelisation fixed
* `perf` and `tune` functions' issue  when choosing the optimum component resolved
-------------------------------------------------------------------------------
## version 6.10.2

### new features / enhancements

* parallel processing on `tune.block.splsda` improved
* `tune.block.splsda` now supports more distances

### bug fixes

* single factor multilevel error in `pls` fixed

-------------------------------------------------------------------------------

## version 6.10.1

### bug fixes

* fixed over-estimated correlation of `cim` for `mixo_(s)pls` objects with single component 
* margin error in `cim` now handled properly
* fixed `plotLoadings` error for very long variable names

-------------------------------------------------------------------------------

## version 6.8.6

### bug fixes

* `predict` function bug for single sample prediction fixed
* `plotLoadings` bug for long variable names fixed

### minor improvements

* missing values in `plotIndiv`'s `group` argument no more throws error
* `mixOmics::predict` function documentation now more accessible

-------------------------------------------------------------------------------

## version 6.8.5

### bug fixes

* names of `linnerud` datasets fixed.


-------------------------------------------------------------------------------

## version 6.8.4

### minor improvements

* package startup message with direct liks to useful resources
* `mixOmics` function documentation disambiguated with instruction on how to get
package help.

-------------------------------------------------------------------------------

## version 6.8.3

### new features / enhancements

* You can now customise `auroc` plots. Refer to documentation for more info.

### bug fixes

* Fixed `tune.spls` and `pef.plsda` bugs when using `cpus` argument for parallel 
processing

### minor improvements

* `auroc` help files now updated with latest changes

-------------------------------------------------------------------------------

## version 6.8.2

### minor improvements

* Updated onLoad message with discussion forum info, bug reports, and more
* Dropped legacy `comp.tol` argument from `pca`

-------------------------------------------------------------------------------

## version 6.8.1

### bug fixes

* `perf.plot` bug in extracting names fixed
* Few fixes for tune.splsda with AUC

### minor improvements

* `plot.perf` now respects `ylim` arguments for custom y range
* Added code of conduct
* Updated DESCRIPTION with bug reports and biocViews
* Updated README

-------------------------------------------------------------------------------

## version 6.8.0

* NOW HOSTED ON BIOCONDUCTOR

-------------------------------------------------------------------------------

## version < 6.8.0

* Refer to `./inst/legacy/NEWS-old` on GitHub repo.