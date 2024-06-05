test_that("vse function test", {
    set.seed(123)
    x <- rnorm(1024)
    x.vse <- vse(x)

    # Add conditions related to your function's expected output
    # vse returns 0.5±0.01
    expect_type(x.vse, "double")
    expect_equal(x.vse, 0.5, tolerance = 0.01)
})
