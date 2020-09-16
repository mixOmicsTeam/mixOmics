test_that(".get.pch works as expected", code = {
    expect_identical(
        .get.pch(pch = NULL, pch.levels = NULL),
        list(
            pch = structure(1L, .Label = "16", class = "factor"),
            pch.levels = c(`16` = 16),
            pch.legend = FALSE
        )
    )
    
    expect_identical(
        .get.pch(pch = mtcars$cyl[1:5]),
        list(
            pch = structure(
                c(2L, 2L, 1L, 2L, 3L),
                .Label = c("4", "6", "8"),
                class = "factor"
            ),
            pch.levels = c(`4` = 16, `6` = 15, `8` = 17),
            pch.legend = TRUE
        )
    )
    
    
    expect_identical(
        .get.pch(
            pch = mtcars$cyl[1:5],
            pch.levels = c(`4` = 21,  `8` = 17, `6` = 11)
        ),
        list(
            pch = structure(
                c(2L, 2L, 1L, 2L, 3L),
                .Label = c("4", "6", "8"),
                class = "factor"
            ),
            pch.levels = c(`4` = 21, `8` = 17, `6` = 11),
            pch.legend = TRUE
        )
    )
    expect_identical(
        .get.pch(pch = 21),
        list(pch = structure(1L, .Label = "21", class = "factor"), pch.levels = c(`21` = 21), pch.legend = FALSE)
    )
    
})
