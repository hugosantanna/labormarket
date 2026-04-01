#' Corrects the Limited Mobility Bias Using Bootstrap
#'
#' The `lmbias` function applies a bootstrapping technique to correct the limited mobility bias
#' (often referred to as small sample bias) in fixed-effects models. This implementation is based
#' on the method proposed by Azkarate-Askasua and Zerecero (2022). It uses a Rademacher random
#' vector and the fixest package's model specifications for its calculations.
#'
#' @param femodel A model object of class `fixest` (from the `fixest` package) that contains the
#'        fixed effects estimates to be corrected. Must have at least two fixed effects.
#' @param data A data.frame or data.table object that contains the dataset used for the fixed effects model.
#' @param vcov A character string specifying the type of heteroskedasticity correction to be
#'        used. Must be one of `"HC0"` or `"HC1"` (default).
#' @param R An integer that sets the number of bootstrap replications. Default is 1000.
#' @param xvarname A character string giving the name of the column in `data` that contains
#'        the linear predictor (X * beta). Default is `"xbeta"`.
#' @param parallel A character string that determines the type of parallel operation to be used
#'        (if any). Possible values are `"no"`, `"multicore"`, and `"snow"`. Default is `"no"`.
#' @param cluster An optional vector or factor that provides the data's cluster IDs if cluster
#'        bootstrapping is desired.
#' @param ncpus An integer specifying the number of CPUs to be used for parallel processing.
#'        This argument is relevant only if `parallel` is set to `"multicore"` or `"snow"`.
#'        Default uses the option set in `getOption("boot.ncpus", 1L)`.
#' @param cl An optional parallel or snow cluster for use with `parallel = "snow"`. If not
#'        supplied and `parallel = "snow"`, a cluster on the local machine is created.
#'
#' @return A 3x4 matrix with rows `"Original"`, `"Bias"`, and `"Corrected"`, and columns for
#'         variance of each fixed effect, covariance, and correlation.
#'
#' @references
#' Azkarate-Askasua, M. and Zerecero, M. (2022). "The Macroeconomics of Firms' Financing."
#' Working Paper. \url{https://github.com/mazkarate/bias_correction}
#'
#' @seealso [simlabormarket()] for simulating labor market data suitable for this function.
#'
#' @import Matrix
#' @importFrom boot boot
#' @importFrom lfe compfactor
#' @import fixest
#' @export lmbias
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' # Simulate a labor market with high mobility for connectivity
#' lm_sim <- simlabormarket(nk = 2, nl = 2, nt = 5, ni = 2000, lambda = 0.5)
#' panel <- lm_sim@panel
#' panel$fid <- as.factor(panel$fid)
#' panel$i   <- as.factor(panel$i)
#'
#' # Restrict to the largest connected set
#' cf <- lfe::compfactor(list(panel$i, panel$fid))
#' largest <- names(which.max(table(cf)))
#' panel <- panel[cf == largest, ]
#' panel$i   <- droplevels(panel$i)
#' panel$fid <- droplevels(panel$fid)
#'
#' # Fit a two-way fixed effects model
#' model <- fixest::feols(lw ~ 1 | i + fid, data = panel)
#'
#' # Add the linear predictor column (match model obs)
#' panel <- panel[!is.na(predict(model, newdata = panel)), ]
#' panel$xbeta <- 0
#'
#' # Correct limited mobility bias (small R for demo)
#' lmbias(model, panel, R = 50)
#' }

lmbias <- function(femodel, data, vcov = "HC1", R = 1000, xvarname = "xbeta",
                    parallel = c("no", "multicore", "snow"), cluster =  NULL,
                    ncpus = getOption("boot.ncpus", 1L), cl = NULL) {

# Checking if the variables are correct
  if (!inherits(femodel, "fixest")) {
    stop("The model must be a fixest object (e.g., from feols())")
  }

  if(!("data.frame" %in% class(data))) {
    stop("The data must be a data.frame or data.table object")
  }

  vcov <- match.arg(vcov, c("HC0", "HC1"))
  parallel <- match.arg(parallel)

  if (!is.numeric(R) || length(R) != 1 || R < 1) {
    stop("R must be a positive integer (number of bootstrap replications)")
  }
  R <- as.integer(R)

  if (!is.character(xvarname) || length(xvarname) != 1) {
    stop("xvarname must be a single character string")
  }

  if (length(femodel$fixef_vars) < 2) {
    stop("The model must have at least two fixed effects")
  }

  # fixed effects names
  fe1 <- femodel$fixef_vars[1]
  fe2 <- femodel$fixef_vars[2]

  # Subset data to observations used by the model (handles singleton removal)
  nobs_model <- femodel$nobs
  if (nrow(data) != nobs_model) {
    obs_removed <- femodel$obs_selection$obsRemoved
    if (!is.null(obs_removed)) {
      data <- data[obs_removed, ]  # negative indices in fixest = rows to drop
    }
    if (nrow(data) != nobs_model) {
      stop("Data row count (", nrow(data), ") does not match model nobs (", nobs_model, "). ",
           "Please pass the same data used to fit the model.")
    }
  }

  cluster_var <- NULL
  if (!is.null(cluster)) {
    if (length(cluster) == 1 && is.character(cluster)) {
      if (!(cluster %in% names(data))) {
        stop("cluster must be either a column name in data or a vector with one value per observation")
      }
      cluster_vec <- data[[cluster]]
    } else if (length(cluster) == nrow(data)) {
      cluster_vec <- cluster
    } else {
      stop("cluster must be either a column name in data or a vector with one value per observation")
    }
    cluster_var <- ".__lmbias_cluster__"
    data[[cluster_var]] <- as.factor(cluster_vec)
  }

  data[[fe1]] <- factor(data[[fe1]])
  data[[fe2]] <- factor(data[[fe2]])

  # fixed effects (originals)
  beta_fe1_original <- fixef(femodel)[[1]]
  beta_fe2_original <- fixef(femodel)[[2]]

  # #obs
  nobs <- femodel$nobs

  nfe1 <- as.integer(femodel$fixef_sizes[1])
  nfe2 <- as.integer(femodel$fixef_sizes[2])
  ncovar <- length(femodel$coefficients)
  k  <- nfe1 + nfe2 + ncovar

  # converting to factors
  data[[fe1]] <- as.factor(data[[fe1]])
  data[[fe2]] <- as.factor(data[[fe2]])

  # checking the connected set
  if(nlevels(compfactor(list(factor(data[[fe1]]), factor(data[[fe2]])))) != 1) {
    stop("The fixed effects must be connected in a single graph. \n
         Please check the data and the model with lfe::compfactor.")
  }

  # left hand side
  lhs <- as.character(formula(femodel))[3]

    # bootstrap function
  X1 <- sparse.model.matrix(
    formula(paste0("~",
      paste(femodel$fixef_vars[1]), "-1"
    )),
    data = data
  )
  X2 <- sparse.model.matrix(
    formula(paste0("~",
      paste(femodel$fixef_vars[2]), "-1"
    )),
    data = data
  )

  # what type of vcov
  if( vcov == "HC0") {
    phi <- sqrt(resid(femodel)^2)
  } else if (vcov == "HC1") {
    phi <- sqrt(resid(femodel)^2 * nobs / (nobs - k))
  }

  # Clustering method
  calculate_bias <- function(data = data, indices, phi, nobs, lhs, nfe1,
                              nfe2, ncovar, X1, X2, cluster_var) {

    # generating the random -1, 1 vetor (rademacher)
    if(is.null(cluster_var)) {

      r  <- sample(c(-1, 1), nobs, replace = TRUE)
      data$r <- r * phi

    } else {

      groups <- unique(data[[cluster_var]])
      g <- length(groups)
      r  <- sample(c(-1, 1), g, replace = TRUE)
      data$r <- r[match(data[[cluster_var]], groups)] * phi

    }


    new_fml <- formula(paste0("r ~ ", lhs))
    model <- feols(new_fml, data = data)
    fe <- fixef(model)
    beta_fe1 <- fe[[1]]
    beta_fe2 <- fe[[2]]


    # creating the function: beta * A * beta
    # I will exploit the fact we only need a scalar in the end to bypass
    #    the gigantic dense matrix in the middle

    # Generating the vector of ones
    ones <- rep(1, nobs)

    # Moment 1 --------------------------------------------------------------
    m1 <- t(beta_fe1) %*% crossprod(X1) %*% beta_fe1
    m1 <- (1/(nobs - 1)) * (m1 - (1/nobs) * (t(beta_fe1) %*% t(X1) %*% ones %*% t(ones) %*% X1 %*% beta_fe1))


    # Moment 2 --------------------------------------------------------------
    m2 <- (t(beta_fe2) %*% crossprod(X2) %*% beta_fe2)
    m2 <- (1/(nobs - 1)) * (m2 - (1/nobs) * (t(beta_fe2) %*% t(X2) %*% ones %*% t(ones) %*% X2 %*% beta_fe2))

    # Moment 3 --------------------------------------------------------------
    m3_1 <- t(beta_fe2) %*% crossprod(X2, X1) %*% beta_fe1
    m3_1 <- (0.5/(nobs - 1)) * (m3_1 - (1/nobs) * (t(beta_fe2) %*% t(X2) %*% ones %*% t(ones) %*% X1 %*% beta_fe1))

    m3_2 <- t(beta_fe1) %*% crossprod(X1, X2) %*% beta_fe2
    m3_2 <- (0.5/(nobs - 1)) * (m3_2 - (1/nobs) * (t(beta_fe1) %*% t(X1) %*% ones %*% t(ones) %*% X2 %*% beta_fe2))

    m3 <- m3_1 + m3_2

    return(c(as.numeric(m1), as.numeric(m2), as.numeric(m3)))
  }

  # Starting the bootstrapping
  boot_mbias <- boot(data, statistic = calculate_bias, phi = phi,
                      nobs = nobs, lhs = lhs, R = R, parallel = parallel,
                      ncpus = ncpus, cl = cl, X1 = X1, X2 = X2, cluster_var = cluster_var)

  boot_res <- colMeans(boot_mbias$t)
  results <- cbind(
    var(beta_fe1_original[femodel$fixef_id[[1]]]),
    var(beta_fe2_original[femodel$fixef_id[[2]]]),
    cov(beta_fe1_original[femodel$fixef_id[[1]]], beta_fe2_original[femodel$fixef_id[[2]]]),
    cor(beta_fe1_original[femodel$fixef_id[[1]]], beta_fe2_original[femodel$fixef_id[[2]]])
  )

  results <- rbind(results, cbind(t(boot_res), (boot_res[3]/(sqrt(boot_res[1]) * sqrt(boot_res[2])))))
  results <- rbind(results, results[1,] - results[2,])

  colnames(results) <- c(
    paste0("Var of ", fe1),
    paste0("Var of ", fe2),
    paste0("Covariance"),
    paste0("Correlation")
  )
  rownames(results) <- c("Original", "Bias", "Corrected")

  return(results)
}
