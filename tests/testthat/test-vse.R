test_that("vse function test", {
    set.seed(123)
    x <- rnorm(1024)
    vse <- vse(x)

    # Add conditions related to your function's expected output
    # vse returns a list with vse = 0.5±0.01
    expect_type(vse, "list")
    expect_equal(vse$vse, 0.5, tolerance = 0.01)
})
