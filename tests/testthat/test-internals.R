context("internals")

test_that(".get.ind.colors works as expected", code = {
    foo <- .get.ind.colors(group = NULL, col = 'black', n_ind = 10)
    expect_true(length(foo) == 10)
    expect_true(identical(unique(foo), 'black'))
    foo <- .get.ind.colors(group = factor(mtcars$am[1:10]), col = NULL, n_ind = 10, 
                           col.per.group = c('1'='red', '0'='blue'))
    expect_true(length(foo) == 10)
    expect_true(all(names(foo) %in% c('red', 'blue')))
    expect_true(all(foo %in% c('0', '1')))
})

test_that(".are.colors works as expected", code = {
    expect_error(.are.colors(c('red', 'foo')))
})

test_that(".get.colors works as expected", code = {
    expect_error(.get.colors(col = c('red', 'foo'), n_col = 10))
    expect_true(length(.get.colors(col = c('red'), n_col = 10)) == 10)
})

test_that(".get.character.vector works as expected", code = {
    expect_identical(.get.character.vector(arg = TRUE, vec = letters[1:10]), letters[1:10])
    expect_identical(.get.character.vector(arg = LETTERS[1:10], vec = letters[1:10]), LETTERS[1:10])
    expect_null(.get.character.vector(arg = FALSE, vec = letters[1:10]))
    expect_error(.get.character.vector(arg = letters, vec = letters[1:10]))
})

test_that(".check_test.keepX works as expected for single and list X", code = {
    expect_equal(.check_test.keepX(test.keepX = c(3, 8, 1), X = mtcars), c(3, 8, 1))
    expect_condition(.check_test.keepX(test.keepX = c(1, 2.2, 3), X = mtcars))
    expect_condition(.check_test.keepX(test.keepX = c(1, '2', 3), X = mtcars))
    expect_condition(.check_test.keepX(test.keepX = seq(3, ncol(mtcars)+1, 1), X = mtcars))
    
    expect_equal(.check_test.keepX(
        test.keepX = list('a'=c(1,4), 'b'=c(3,5,8)), 
        X =          list('a'=mtcars, 'b'=mtcars)
    ),
    list('a'=c(1,4), 'b'=c(3,5,8)))
    
    expect_condition(.check_test.keepX(
        test.keepX = list('a'=c(1,4), 'b'=c(3,5,8)),
        X =          list('a'=mtcars, 'b'=mtcars, 'c'=mtcars  )
    ))
    expect_condition(.check_test.keepX(
        test.keepX = list('a'=c(1,4), 'b'=c(3,5,8), 'c'=c(4,8)),
        X =          list('a'=mtcars, 'b'=mtcars)
    ))
    expect_condition(.check_test.keepX(
        test.keepX = list('a'=c(1,4), 'b'=c(3,ncol(mtcars)+1)), 
        X =          list('a'=mtcars, 'b'=mtcars               )
    ))
    
})
test_that(".check_ncomp works as expected for single and list X", code = {
    expect_true(.check_ncomp(ncomp=2, X = mtcars) == 2)
    expect_true(.check_ncomp(ncomp=2, X = list(a=mtcars, b=mtcars)) == 2)
    expect_condition(.check_ncomp(ncomp=min(dim(mtcars))+1, X = mtcars))
    expect_condition(.check_ncomp(ncomp=4, X = list(a=mtcars, b=mtcars[,1:3])))
    
})
