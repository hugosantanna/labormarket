make_lmbias_fixture <- function() {
  set.seed(101)
  counts <- matrix(
    c(
      12, 2, 1, 1,
      2, 12, 2, 1,
      1, 2, 12, 2,
      1, 1, 2, 12,
      6, 4, 2, 1
    ),
    nrow = 5,
    byrow = TRUE
  )
  cells <- which(counts > 0, arr.ind = TRUE)
  cell_id <- rep(seq_len(nrow(cells)), counts[cells])
  df <- data.frame(
    f1 = factor(cells[cell_id, 1]),
    f2 = factor(cells[cell_id, 2])
  )
  df$x <- rnorm(nrow(df))
  fe1 <- c(-2, -1, 0, 1, 2)
  fe2 <- c(-1.5, -0.5, 0.5, 1.5)
  df$y <- fe1[as.integer(df$f1)] + fe2[as.integer(df$f2)] +
    0.4 * df$x + rnorm(nrow(df), sd = 0.2)
  df$xbeta <- 0
  df$cluster_id <- rep(seq_len(ceiling(nrow(df) / 5)), each = 5,
                       length.out = nrow(df))

  list(
    data = df,
    model = fixest::feols(y ~ x | f1 + f2, data = df)
  )
}

test_that("lmbias errors on non-fixest model", {
  expect_error(lmbias(lm(y ~ x, data = data.frame(y = 1:10, x = 1:10)),
                       data.frame(y = 1:10, x = 1:10)),
               "fixest object")
})

test_that("lmbias errors on bad vcov", {
  skip_if_not_installed("fixest")
  # Create a minimal fixest model for validation testing
  df <- data.frame(y = rnorm(20), x = rnorm(20), f1 = rep(1:4, 5), f2 = rep(1:5, 4))
  df$xbeta <- df$x
  mod <- fixest::feols(y ~ x | f1 + f2, data = df)
  expect_error(lmbias(mod, df, vcov = "HC3"), "'arg' should be one of")
})

test_that("lmbias no longer requires the unused xvarname column", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("lfe")
  fixture <- make_lmbias_fixture()
  fixture$data$xbeta <- NULL
  expect_no_error(lmbias(fixture$model, fixture$data, R = 2, vcov = "HC0"))
})

test_that("lmbias errors on model with fewer than 2 FEs", {
  skip_if_not_installed("fixest")
  df <- data.frame(y = rnorm(20), x = rnorm(20), f1 = rep(1:4, 5))
  df$xbeta <- df$x
  mod <- fixest::feols(y ~ x | f1, data = df)
  expect_error(lmbias(mod, df), "exactly two fixed effects")
})

test_that("HC1 uses the connected two-way FE rank", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("lfe")
  fixture <- make_lmbias_fixture()

  set.seed(202)
  hc0 <- lmbias(fixture$model, fixture$data, vcov = "HC0", R = 6)
  set.seed(202)
  hc1 <- lmbias(fixture$model, fixture$data, vcov = "HC1", R = 6)

  n <- nrow(fixture$data)
  n_fe1 <- nlevels(fixture$data$f1)
  n_fe2 <- nlevels(fixture$data$f2)
  n_slopes <- length(fixture$model$coefficients)
  model_rank <- n_fe1 + n_fe2 - 1 + n_slopes
  expected_hc1_factor <- n / (n - model_rank)

  expect_equal(
    unname(hc1["Bias", 1:3]),
    unname(hc0["Bias", 1:3]) * expected_hc1_factor,
    tolerance = 1e-8
  )
})

test_that("corrected correlation is computed from corrected moments", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("lfe")
  fixture <- make_lmbias_fixture()

  set.seed(303)
  result <- lmbias(fixture$model, fixture$data, R = 10)
  corrected <- result["Corrected", ]
  expected <- corrected[3] / sqrt(corrected[1] * corrected[2])

  expect_gt(corrected[1], 0)
  expect_gt(corrected[2], 0)
  expect_equal(unname(corrected[4]), unname(expected), tolerance = 1e-12)
})

test_that("cluster accepts a column name or an aligned vector", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("lfe")
  fixture <- make_lmbias_fixture()

  set.seed(404)
  by_name <- lmbias(fixture$model, fixture$data, R = 4,
                    cluster = "cluster_id")
  set.seed(404)
  by_vector <- lmbias(fixture$model, fixture$data, R = 4,
                      cluster = fixture$data$cluster_id)

  expect_equal(by_name, by_vector)
})

test_that("lmbias returns 3x4 matrix with correct names", {
  skip_on_cran()
  skip_if_not_installed("fixest")
  skip_if_not_installed("lfe")

  set.seed(42)
  lm_obj <- simlabormarket(nk = 2, nl = 2, nt = 5, ni = 2000, lambda = 0.5)
  panel <- lm_obj@panel
  panel$fid <- as.factor(panel$fid)
  panel$i   <- as.factor(panel$i)

  # Restrict to the largest connected set
  cf <- lfe::compfactor(list(panel$i, panel$fid))
  largest <- names(which.max(table(cf)))
  panel <- panel[cf == largest, ]
  panel$i   <- droplevels(panel$i)
  panel$fid <- droplevels(panel$fid)

  mod <- fixest::feols(lw ~ 1 | i + fid, data = panel)
  panel$xbeta <- 0  # intercept-only model

  res <- lmbias(mod, panel, R = 10)

  expect_true(is.matrix(res))
  expect_equal(dim(res), c(3, 4))
  expect_equal(rownames(res), c("Original", "Bias", "Corrected"))
  expect_equal(colnames(res)[3], "Covariance")
  expect_equal(colnames(res)[4], "Correlation")
})
