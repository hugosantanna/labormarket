test_that("simlabormarket returns LaborMarket S4 object", {
  lm_obj <- simlabormarket(nk = 3, nl = 3, nt = 2, ni = 50, lambda = 0.3)
  expect_s4_class(lm_obj, "LaborMarket")
})

test_that("panel has expected columns and dimensions", {
  ni <- 80
  nt <- 3
  lm_obj <- simlabormarket(nk = 2, nl = 2, nt = nt, ni = ni, lambda = 0.2)
  panel <- lm_obj@panel

  expect_true(is.data.frame(panel))
  expect_equal(nrow(panel), ni * nt)
  expect_true(all(c("i", "k", "t", "fid", "spell", "gender", "psi", "alpha", "lw") %in% names(panel)))
})

test_that("mincer = TRUE populates age, educ, experience", {
  lm_obj <- simlabormarket(nk = 2, nl = 2, nt = 2, ni = 50, lambda = 0.3, mincer = TRUE)
  panel <- lm_obj@panel

  expect_true(all(c("age", "educ", "experience") %in% names(panel)))
  expect_true(all(panel$age > 0))
  expect_true(all(panel$educ > 0))
})

test_that("gender is binary (0/1)", {
  lm_obj <- simlabormarket(nk = 2, nl = 2, nt = 2, ni = 100, lambda = 0.2)
  expect_true(all(lm_obj@panel$gender %in% c(0, 1)))
})

test_that("worker and firm effects are fixed across the panel", {
  lm_obj <- simlabormarket(nk = 3, nl = 3, nt = 4, ni = 100, lambda = 0.3)
  panel <- lm_obj@panel

  worker_alpha_counts <- panel[, .(n_alpha = data.table::uniqueN(alpha)), by = i]
  firm_psi_counts <- panel[, .(n_psi = data.table::uniqueN(psi)), by = fid]

  expect_true(all(worker_alpha_counts$n_alpha == 1))
  expect_true(all(firm_psi_counts$n_psi == 1))
})

test_that("input validation errors on bad params", {
  expect_error(simlabormarket(nk = 1), "nk must be an integer >= 2")
  expect_error(simlabormarket(nl = 0), "nl must be an integer >= 2")
  expect_error(simlabormarket(nt = 1), "nt must be an integer >= 2")
  expect_error(simlabormarket(ni = -1), "ni must be a positive integer")
  expect_error(simlabormarket(ratiog = 0), "ratiog must be a number strictly between 0 and 1")
  expect_error(simlabormarket(ratiog = 1), "ratiog must be a number strictly between 0 and 1")
  expect_error(simlabormarket(lambda = 0), "lambda must be a number strictly between 0 and 1")
  expect_error(simlabormarket(lambda = 1.5), "lambda must be a number strictly between 0 and 1")
  expect_error(simlabormarket(mincer = "yes"), "mincer must be TRUE or FALSE")
})
