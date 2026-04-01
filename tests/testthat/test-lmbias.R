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

test_that("lmbias errors on malformed xvarname input", {
  skip_if_not_installed("fixest")
  df <- data.frame(y = rnorm(20), x = rnorm(20), f1 = rep(1:4, 5), f2 = rep(1:5, 4))
  mod <- fixest::feols(y ~ x | f1 + f2, data = df)
  expect_error(lmbias(mod, df, xvarname = c("x", "z")), "single character string")
})

test_that("lmbias errors on model with fewer than 2 FEs", {
  skip_if_not_installed("fixest")
  df <- data.frame(y = rnorm(20), x = rnorm(20), f1 = rep(1:4, 5))
  df$xbeta <- df$x
  mod <- fixest::feols(y ~ x | f1, data = df)
  expect_error(lmbias(mod, df), "at least two fixed effects")
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

test_that("lmbias accepts cluster vectors", {
  skip_on_cran()
  skip_if_not_installed("fixest")
  skip_if_not_installed("lfe")

  set.seed(42)
  lm_obj <- simlabormarket(nk = 2, nl = 2, nt = 5, ni = 1500, lambda = 0.5)
  panel <- lm_obj@panel
  panel$fid <- as.factor(panel$fid)
  panel$i   <- as.factor(panel$i)

  cf <- lfe::compfactor(list(panel$i, panel$fid))
  largest <- names(which.max(table(cf)))
  panel <- panel[cf == largest, ]
  panel$i   <- droplevels(panel$i)
  panel$fid <- droplevels(panel$fid)

  mod <- fixest::feols(lw ~ 1 | i + fid, data = panel)

  res <- lmbias(mod, panel, R = 5, cluster = panel$fid)
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(3, 4))
})
