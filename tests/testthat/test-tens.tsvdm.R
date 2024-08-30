context("tsvdm")
# bltodo: add tests for permutations of different configurations and input sizes
# like the `tred` library

test_that(
  "reconstructing original matrix from tsvdm",
  code = {
    test_tensor <- array(1:24, dim = c(2, 4, 3))
    tsvdm_decomposition <- tsvdm(test_tensor)
    u <- tsvdm_decomposition$u
    s <- tsvdm_decomposition$s
    v <- tsvdm_decomposition$v
    vt <- ft(v)

    expect_equal(
      test_tensor,
      m_product(u, s, vt)
    )
  }
)

test_that(
  "reconstructing original matrix from tsvdm without transform",
  code = {
    test_tensor <- array(1:24, dim = c(2, 4, 3))
    tsvdm_decomposition <- tsvdm(test_tensor, transform = FALSE)
    u <- tsvdm_decomposition$u
    s <- tsvdm_decomposition$s
    v <- tsvdm_decomposition$v
    vt <- ft(v)

    expect_equal(
      test_tensor,
      facewise_product(u, s, vt)
    )
  }
)
