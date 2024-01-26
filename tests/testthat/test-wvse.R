test_that("wvse function test", {
    set.seed(123)
    x <- rnorm(1024)
    vse <- wvse(x)

    # Add conditions related to your function's expected output
    expect_type(vse, "double")
    expect_equal(vse, 0.5, tolerance = 0.01)
})
