############# build all documentation for the package and update namespace
# file.remove(list.files('man/', pattern = '.Rd', full.names = TRUE))
roxygen2::roxygenise()
############# load the pkg to access functions
devtools::load_all()
fun = 'spca'

# fun = NULL
test_fun <- function(fun=fun){
  if (!is.null(fun)) {
    ## run file tests
    devtools::test_file(file = sprintf('tests/testthat/test-%s.R', fun))
    ## run examples
    source(sprintf('examples/%s-example.R', fun))
  }
}
test_fun(fun = fun)

# devtools::test_file()
############# run ./tests tests which then redirects to testthat.R. testhat/helper-*.R files are run before tests
devtools::test()
############# check without examples
devtools::check(args = "--no-examples")
############# test examples
devtools::test_file(file = 'tests/testthat/test-pca.R')
## path to Rd files - will let you know if there are any errors/warnings
# testthat::test_example("man/parent_base_ext.Rd")
############# manual tests
## test files that cannot be included in standard testthat pipeline
# invisible(lapply(list.files("tests/manual", full.names = TRUE), source))
############# check
# devtools::check()
############# build it
## by default, it is made in parent directoy of the package, we manually set it.
## testhat has changed wd
# setwd('/Users/alabadi/Documents/_Projects/_Personal/someWrappers')

  # binary_dir <- "../__binary"
  # if(!dir.exists(binary_dir)){
  #   dir.create(binary_dir)
  # }
  # pkg <- devtools::build(path = binary_dir)
  # ############# install it in my library
  # install.packages(pkg, repos = NULL)

# CMD + SHift + B
