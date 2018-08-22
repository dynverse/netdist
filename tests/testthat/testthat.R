context("Check nettools consistency.")


set.seed(0)
a <- matrix(rnorm(1000),ncol=100)
set.seed(1)
b <- matrix(rnorm(1000),ncol=100)

test_that("Test empty-full distance:", {
    e <- matrix(0, ncol=100, nrow=100)
    f <- e + 1
    diag(f) <- 0
    dd <- netdist(e, f, d="HIM", components=FALSE, n.cores=1)
    expect_equal(as.numeric(dd), 1)
})

test_that("Test distance:", {
  e <- matrix(runif(100*100), ncol=100, nrow=100)
  f <- matrix(runif(100*100), ncol=100, nrow=100)
  d1 <- netdist(e,f, d="HIM", components=FALSE, n.cores=1)
})