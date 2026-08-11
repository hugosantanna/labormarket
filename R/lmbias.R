#' Corrects the Limited Mobility Bias Using Bootstrap
#'
#' The `lmbias` function applies a bootstrapping technique to correct the limited mobility bias
#' (often referred to as small sample bias) in fixed-effects models. This implementation is based
#' on the method proposed by Azkarate-Askasua and Zerecero (2022). It uses a Rademacher random
#' vector and the fixest package's model specifications for its calculations.
#'
#' @param femodel A model object of class `fixest` (from the `fixest` package) that contains the
#'        fixed effects estimates to be corrected. The model must have exactly two additive
#'        fixed effects whose worker--firm graph is connected.
#' @param data A data.frame or data.table object that contains the dataset used for the fixed effects model.
#' @param vcov A character string specifying the type of heteroskedasticity correction to be
#'        used. Must be one of `"HC0"` or `"HC1"` (default).
#' @param R An integer that sets the number of bootstrap replications. Default is 1000.
#' @param xvarname Deprecated and retained for backward compatibility. The wild bootstrap
#'        refits the original model specification and does not require a precomputed linear
#'        predictor.
#' @param parallel A character string that determines the type of parallel operation to be used
#'        (if any). Possible values are `"no"`, `"multicore"`, and `"snow"`. Default is `"no"`.
#' @param cluster An optional character string naming a cluster-ID column in `data`, or a
#'        vector/factor with one cluster ID per row of `data`. One Rademacher draw is used per
#'        cluster in each replication.
#' @param ncpus An integer specifying the number of CPUs to be used for parallel processing.
#'        This argument is relevant only if `parallel` is set to `"multicore"` or `"snow"`.
#'        Default uses the option set in `getOption("boot.ncpus", 1L)`.
#' @param cl An optional parallel or snow cluster for use with `parallel = "snow"`. If not
#'        supplied and `parallel = "snow"`, a cluster on the local machine is created.
#'
#' @return A 3x4 matrix with rows `"Original"`, `"Bias"`, and `"Corrected"`, and columns for
#'         variance of each fixed effect, covariance, and correlation. The corrected
#'         correlation is `NA` if either corrected variance is nonpositive.
#'
#' @references
#' Azkarate-Askasua, M. and Zerecero, M. (2022). "The Macroeconomics of Firms' Financing."
#' Working Paper. \url{https://github.com/mazkarate/bias_correction}
#'
#' @seealso [simlabormarket()] for simulating labor market data suitable for this function.
#'
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
#' # Match the model observations
#' panel <- panel[!is.na(predict(model, newdata = panel)), ]
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

  if (!is.data.frame(data)) {
    stop("The data must be a data.frame or data.table object")
  }

  vcov <- match.arg(vcov, c("HC0", "HC1"))
  parallel <- match.arg(parallel)

  if (!is.numeric(R) || length(R) != 1L || !is.finite(R) ||
      R < 1 || R != floor(R)) {
    stop("R must be a positive integer (number of bootstrap replications)")
  }
  R <- as.integer(R)

  # Kept in the signature so existing calls continue to work. Earlier versions
  # constructed a matrix from xvarname but never used it in the estimator.
  if (!missing(xvarname) && (!is.character(xvarname) || length(xvarname) != 1L)) {
    stop("xvarname must be a single character string")
  }

  model_method <- if (!is.null(femodel$method_type)) {
    femodel$method_type
  } else {
    femodel$method
  }
  if (!is.null(model_method) && !identical(model_method, "feols")) {
    stop("The model must be estimated with fixest::feols()")
  }

  if (length(femodel$fixef_vars) != 2L) {
    stop("The model must have exactly two fixed effects")
  }

  if (!is.null(femodel$slope_flag) && any(femodel$slope_flag != 0L)) {
    stop("Varying-slope fixed effects are not supported")
  }

  # fixed effects names
  fe1 <- femodel$fixef_vars[1]
  fe2 <- femodel$fixef_vars[2]

  # Resolve clusters before applying the model's observation selection so that
  # a vector aligned with the original data is subset in the same way.
  nobs_model <- femodel$nobs
  input_nobs <- nrow(data)
  cluster_values <- NULL
  if (!is.null(cluster)) {
    if (is.character(cluster) && length(cluster) == 1L && cluster %in% names(data)) {
      cluster_values <- data[[cluster]]
    } else if ((is.atomic(cluster) || is.factor(cluster)) &&
               length(cluster) %in% c(input_nobs, nobs_model)) {
      cluster_values <- cluster
    } else {
      stop("cluster must be a column name in data or a vector with one value per observation")
    }
  }

  # Subset data to observations used by the model (handles NA and singleton removal).
  if (nrow(data) != nobs_model) {
    obs_removed <- femodel$obs_selection$obsRemoved
    if (!is.null(obs_removed)) {
      kept_rows <- seq_len(input_nobs)[obs_removed]
      data <- data[obs_removed, ]  # negative indices in fixest = rows to drop
      if (!is.null(cluster_values) && length(cluster_values) == input_nobs) {
        cluster_values <- cluster_values[kept_rows]
      }
    }
    if (nrow(data) != nobs_model) {
      stop("Data row count (", nrow(data), ") does not match model nobs (", nobs_model, "). ",
           "Please pass the same data used to fit the model.")
    }
  }

  if (!is.null(cluster_values)) {
    if (length(cluster_values) != nobs_model) {
      stop("cluster must have one value per observation used by the model")
    }
    if (anyNA(cluster_values)) {
      stop("cluster IDs cannot contain missing values")
    }
  }

  data[[fe1]] <- factor(data[[fe1]])
  data[[fe2]] <- factor(data[[fe2]])

  # A connected two-way FE design has rank n_fe1 + n_fe2 - 1. The
  # normalization removes one, rather than one per FE dimension.
  if (nlevels(compfactor(list(data[[fe1]], data[[fe2]]))) != 1L) {
    stop("The fixed effects must be connected in a single graph. \n
         Please check the data and the model with lfe::compfactor.")
  }

  original_fe <- fixef(femodel)
  beta_fe1_original <- original_fe[[1]]
  beta_fe2_original <- original_fe[[2]]

  nobs <- femodel$nobs
  nfe1 <- length(beta_fe1_original)
  nfe2 <- length(beta_fe2_original)
  ncovar <- length(femodel$coefficients)
  k <- nfe1 + nfe2 - 1L + ncovar

  if (vcov == "HC1" && nobs <= k) {
    stop("HC correction is undefined because the model has no residual degrees of freedom")
  }

  model_rhs <- as.character(formula(femodel))[3L]

  # what type of vcov
  if (vcov == "HC0") {
    phi <- abs(resid(femodel))
  } else {
    phi <- abs(resid(femodel)) * sqrt(nobs / (nobs - k))
  }

  calculate_bias <- function(data, phi, nobs, model_rhs, cluster_values) {
    # The wild bootstrap keeps the design fixed. Ordinary bootstrap indices
    # were deliberately ignored in earlier versions; sim = "parametric" below
    # now expresses that design directly and avoids generating unused indices.
    if (is.null(cluster_values)) {
      r <- sample(c(-1, 1), nobs, replace = TRUE)
    } else {
      groups <- unique(cluster_values)
      g <- length(groups)
      group_r <- sample(c(-1, 1), g, replace = TRUE)
      r <- group_r[match(cluster_values, groups)]
    }

    data$.__lmbias_response__ <- r * phi
    new_fml <- formula(paste0(".__lmbias_response__ ~ ", model_rhs))
    model <- feols(new_fml, data = data, notes = FALSE)
    fe <- fixef(model)

    fe1_by_observation <- fe[[1]][model$fixef_id[[1]]]
    fe2_by_observation <- fe[[2]][model$fixef_id[[2]]]

    c(
      var(fe1_by_observation),
      var(fe2_by_observation),
      cov(fe1_by_observation, fe2_by_observation)
    )
  }

  # Parametric mode is used as a parallel replication engine: calculate_bias
  # generates the Rademacher outcome while retaining the original design.
  boot_mbias <- boot(data, statistic = calculate_bias, phi = phi,
                     nobs = nobs, model_rhs = model_rhs, R = R,
                     sim = "parametric", parallel = parallel,
                     cluster_values = cluster_values, ncpus = ncpus, cl = cl)

  boot_res <- unname(colMeans(boot_mbias$t))
  original_moments <- c(
    var(beta_fe1_original[femodel$fixef_id[[1]]]),
    var(beta_fe2_original[femodel$fixef_id[[2]]]),
    cov(beta_fe1_original[femodel$fixef_id[[1]]],
        beta_fe2_original[femodel$fixef_id[[2]]])
  )
  corrected_moments <- original_moments - boot_res

  correlation_from_moments <- function(moments) {
    if (any(!is.finite(moments)) || moments[1] <= 0 || moments[2] <= 0) {
      return(NA_real_)
    }
    moments[3] / sqrt(moments[1] * moments[2])
  }

  corrected_correlation <- correlation_from_moments(corrected_moments)
  if (is.na(corrected_correlation) &&
      (corrected_moments[1] <= 0 || corrected_moments[2] <= 0)) {
    warning("Corrected correlation is undefined because a corrected variance is nonpositive")
  }

  results <- rbind(
    Original = c(original_moments, correlation_from_moments(original_moments)),
    Bias = c(boot_res, correlation_from_moments(boot_res)),
    Corrected = c(corrected_moments, corrected_correlation)
  )

  colnames(results) <- c(
    paste0("Var of ", fe1),
    paste0("Var of ", fe2),
    paste0("Covariance"),
    paste0("Correlation")
  )
  results
}
