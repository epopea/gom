library(gom)
data <- data.frame(x1 = round(stats::runif(n = 100, 1, 2), 0),
                   x2 = round(stats::runif(n = 100, 1, 3), 0),
                   Id = 1:100)

gomml <- gom_ml(data.object = data, case.id = "Id", initial.lambda = "random", MC_iter = 300)

test_that("Length of gom-ml object equals the number of profiles", {
  expect_equal(length(gomml), 2)
})


