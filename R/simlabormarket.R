#' Simulate a Labor Market with Gender Bias
#'
#' This function creates a simulated labor market with parameters that control
#' the number of firm types, the fraction of females in the labor market,
#' percentage of movers, etc. The simulated panel includes worker and firm
#' fixed effects together with optionally Mincer-style human capital variables
#' (age, education, experience).
#'
#' @param nk Integer >= 2. Number of firm types or clusters.
#' @param ratiog Numeric in (0, 1). Fraction of females in the labor market.
#' @param lambda Numeric in (0, 1). Probability of a worker moving firms each period.
#' @param nl Integer >= 2. Number of worker types or clusters.
#' @param nt Integer >= 2. Number of time periods.
#' @param ni Positive integer. Number of individuals.
#' @param mincer Logical. If `TRUE`, wages are generated using a Mincer equation
#'        (returns on education, experience, and experience squared) together with
#'        worker and firm random effects. Default is `FALSE`.
#'
#' @return A [LaborMarket-class] S4 object containing the simulated panel and parameters.
#'
#' @seealso [lmbias()] for correcting limited mobility bias using the simulated data,
#'   [LaborMarket-class] for the S4 class definition.
#'
#' @importFrom reshape2 melt
#' @importFrom data.table data.table setnames setkey
#' @importFrom crayon blue green red yellow
#' @import methods
#' @importFrom stats cor cov dbeta dexp dnorm formula qnorm rbeta resid rnorm runif var
#' @export simlabormarket
#'
#' @examples
#' # Simulate a small labor market
#' lm_obj <- simlabormarket(nk = 3, nl = 3, nt = 3, ni = 100, lambda = 0.3)
#'
#' # With Mincer wage equation
#' lm_mincer <- simlabormarket(nk = 3, nl = 3, nt = 3, ni = 100,
#'                              lambda = 0.3, mincer = TRUE)
#'

simlabormarket <- function(nk = 6, ratiog = 0.45, lambda = 0.05, nl = 10, nt = 4, ni = 100000, mincer = FALSE) {

  # Input validation
  if (!is.numeric(nk) || length(nk) != 1 || nk < 2 || nk != as.integer(nk)) {
    stop("nk must be an integer >= 2")
  }
  if (!is.numeric(nl) || length(nl) != 1 || nl < 2 || nl != as.integer(nl)) {
    stop("nl must be an integer >= 2")
  }
  if (!is.numeric(nt) || length(nt) != 1 || nt < 2 || nt != as.integer(nt)) {
    stop("nt must be an integer >= 2")
  }
  if (!is.numeric(ni) || length(ni) != 1 || ni < 1 || ni != as.integer(ni)) {
    stop("ni must be a positive integer")
  }
  if (!is.numeric(ratiog) || length(ratiog) != 1 || ratiog <= 0 || ratiog >= 1) {
    stop("ratiog must be a number strictly between 0 and 1")
  }
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0 || lambda >= 1) {
    stop("lambda must be a number strictly between 0 and 1")
  }
  if (!is.logical(mincer) || length(mincer) != 1) {
    stop("mincer must be TRUE or FALSE")
  }

  # Labor Market Inner Parameters
  #-------------------------------------------------------------------------------------------------

  # alpha specific sd (mixture model)
  alpha_sd = rnorm(nl)
  # firm fixed effect standard deviation
  psi_sd = 1
  # Drawing the random effects averages for firms
  psi_mean = qnorm(1:nk/(nk+1)) * psi_sd
  # Drawing the match effects seeds for individuals
  alpha_mean = qnorm(1:nl/(nl+1))

  # sorting effect
  csort = 0.5
  # network effect
  cnetw = 0.2
  # probability of moving standard deviation
  csig = 0.5
	# Starting point for making firms
	fsize = 10
  # standard deviation for wages
  w_sigma = 0.5
  # maximum years of education
  neduc = 24


  # Matrix Building
  #-------------------------------------------------------------------------------------------------

  # imposing a sorting bias
  sort_gap = min(psi_mean)/10

  # worker gender (male = 1/female = 2)
  ng = 2

  # Generating the transition matrix based on the sorting effect and network effect
  # we expect closer to zero values to have higher moving probability

  # transition matrix
  G = array(0,c(nl, nk, nk, ng))
  for (g in 1:ng) for (l in 1:nl) for (k in 1:nk){

    if(g == 2) {
      G[l, k, , g] = dnorm(psi_mean  - cnetw * psi_mean[k] - csort * alpha_mean[l] - sort_gap, sd = csig)
    } else {
      G[l, k, ,g] = dnorm(psi_mean - cnetw * psi_mean[k] - csort * alpha_mean[l], sd = csig)
    }
    # normalize to get transition matrix
    G[l, k, , g] = G[l, k, , g]/sum(G[l, k, , g])
  }

  # Steady State Matrix
  #-------------------------------------------------------------------------------------------------

  # we then solve for the stationary distribution over psis for each alpha value
  H = array(1/nk,c(nl,nk, ng)) # make a initial matrix for L x K
  for (g in 1:2) for (l in 1:nl) for (i in 1:100) {
    H[l, ,g] = t(G[l, , ,g]) %*% H[l, ,g]
  }

  # creating the education matrix
  #-------------------------------------------------------------------------------------------------

  # generating the matrix
  educ_prob = matrix(0, nrow = nl, ncol = neduc)

  # Parameters for the beta distribution for the extreme types
  alpha_low = 2  # Lower alpha for lower education (type 1)
  beta_high = 1.5    # Higher beta for lower education (type 1)
  alpha_high = 2.5   # Higher alpha for higher education (type 10)
  beta_low = 1   # Lower beta for higher education (type 10)

  # generating quantiles
  quantiles <- seq(1/(2*neduc), 1 - 1/(2*neduc), length.out = neduc)

  # Interpolating alpha and beta parameters for the types
  alphas = seq(alpha_low, alpha_high, length.out = nl)
  betas = seq(beta_high, beta_low, length.out = nl)


  # populating the matrix
  for (l in 1:nl) {

    # Generate probabilities for each education year for the individual type
    probs = dbeta(quantiles, alphas[l], betas[l])
    educ_prob[l,] = probs
  
  }

  # normalizing the matrix rowwise
  educ_prob = educ_prob/rowSums(educ_prob)

  # we simulate a panel
  network    = array(0,c(ni,nt))
  spellcount = array(0,c(ni,nt))
  age        = array(0,c(ni,nt))
  experience = array(0,c(ni,nt))

  # age where we start education
  start_age = 6

  # time invariante worker matrix (worker type, gender, education)
  A = array(0,c(ni, 3))

  for (i in 1:ni) {

    # we draw the gender type
    g = ifelse(runif(1) >= ratiog, 1, 2)
    # we draw the worker type
    l = sample.int(nl,1)

    if (mincer == TRUE) {
      # Get the probability weights for education
      prob_educ = educ_prob[l, ]
      # we draw the education based on worker type
      e = sample(x = 1:ncol(educ_prob), size = 1, prob = prob_educ)
      
      # We build the matrix A
      A[i,] = c(l, g, e)
    } else {
      A[i,] = c(l, g, 0)
    }

    # at time 1, we draw from H, set other initial paramemters
    # network matrix stands for the "history" of worker jobs

    # network matrix -----------------------------------------------------------
    network[i,1] = sample.int(nk,1,prob = H[l, , g])

    if (mincer) {
      # age matrix ---------------------------------------------------------------
      # Generate age pyramid that follows an exponential distribution
      min_age = max(18, e + start_age)
      p_age = dexp(seq_along(min_age:65), rate = 0.05)
      p_age = p_age / sum(p_age)
      age[i, 1] = sample(min_age:65, 1, prob = p_age)

      # experience matrix ---------------------------------------------------------
      experience[i, 1] = age[i, 1] - e - start_age
    }

    # updating spellcount and network matrix ------------------------------------
    for (t in 2:nt) {
      if (runif(1) < lambda) {
        network[i, t] = sample.int(nk, 1, prob = G[l, network[i, t - 1], , g])
        spellcount[i, t] = spellcount[i, t - 1] + 1
        if (mincer) experience[i, t] = 0
      } else {
        network[i, t]    = network[i, t - 1]
        spellcount[i, t] = spellcount[i, t - 1]
        if (mincer) experience[i, t] = experience[i, t - 1] + 1
      }
      if (mincer) age[i, t] = age[i, t - 1] + 1
    }

  }

  # Creating the dataset
  #-------------------------------------------------------------------------------------------------

	data  = data.table(melt(network, c('i', 't')))
	setnames(data, "value", "k")
	data[, spell := data.table(melt(spellcount, c('i', 't')))$value]
  if (mincer) {
    data[, age := data.table(melt(age, c('i', 't')))$value]
    data[, experience := data.table(melt(experience, c('i', 't')))$value]
  }

  # worker, education level and gender
  data[, l := A[i,1], i]
  data[, gender := A[i,2], i]
  data[, educ := A[i,3], i]

  # generating firm ids
  #-------------------------------------------------------------------------------------------------

  last_max_fid = 0
  f_class_count = ni/(fsize*nk*nt) # number of firm classes

  # creating the spell counter dataset (keep track whenever an individual moves from a firm to another)
  dspell = data[,list(len=.N),list(i,spell,k)]

  # Calculate and assign unique fids for each k
  dspell[, fid := {
    current_max_fid = pmax(1, sum(len) / f_class_count)
    fid = sample((last_max_fid + 1):(last_max_fid + current_max_fid), .N, replace=TRUE)
    last_max_fid <<- max(fid)
    .(fid)
  }, by=k]

  setkey(data,i,spell)
  setkey(dspell,i,spell)
    
  data[, fid := dspell[data, fid]]

  # generating fixed worker and firm effects
  #-------------------------------------------------------------------------------------------------
  
  # bargaining gap
  barg_gap = 0.9

  # worker fixed effects are assigned once per individual and then carried
  # across the full panel. A bargaining gap is introduced through the mean.
  worker_effects = unique(data[, .(i, l, gender)])
  worker_effects[, alpha_mean_i := alpha_mean[l]]
  worker_effects[gender == 2, alpha_mean_i := alpha_mean_i * barg_gap]
  worker_effects[, alpha := rnorm(.N, alpha_mean_i, sd = abs(alpha_sd[l]))]

  # firm fixed effects are assigned once per firm identifier and inherited by
  # all observations linked to that firm across periods.
  firm_effects = unique(data[, .(fid, k)])
  firm_effects[, psi := rnorm(.N, psi_mean[k], sd = 0.5)]

  data[worker_effects[, .(i, alpha)], alpha := i.alpha, on = "i"]
  data[firm_effects[, .(fid, psi)], psi := i.psi, on = "fid"]

  # as a final step, I shift the entire effects to have a positive minimum
  # For alpha
  min_alpha = min(data$alpha)
  if(min_alpha < 0) {
      data[, alpha := alpha - min_alpha + 0.01]
  }

  # For psi
  min_psi = min(data$psi)
  if(min_psi < 0) {
      data[, psi := psi - min_psi + 0.01]
  }
  if(mincer == TRUE) {
    data = data[, .(i, l, k, t, fid, spell, gender, age, educ, experience, psi, alpha)]
  } else {
    data = data[, .(i,l,k,t,fid,spell, gender, psi, alpha)]
  }

  # female is 1, male is 0 (making the var binary)
  data$gender = data$gender - 1

  # Final wrangling for data creation
  #-------------------------------------------------------------------------------------------------

  # year shock dummies, just an uniform distribution. For simplicity, my model is static for now
  shocks = runif(nt, -0.01, 0.05)
  data[, time_fe := shocks[t]]

  if (mincer == TRUE) {
    data[, lw  := 0.07 * educ + 0.02 * experience - 0.0005 * experience^2 + 
                  alpha + psi + time_fe + w_sigma * rnorm(.N)]
  } else {
    data[, lw  := alpha + psi + w_sigma * rnorm(.N)]
  }
  

  # Setting education, age and experience to be zeroes if mincer == FALSE
  if (mincer == FALSE) {
    data[, c("educ", "age", "experience") := 0]
  }

  # Creating the labor market simulator type of class
  #-------------------------------------------------------------------------------------------------


  # Define the constructor
  create_labor_m = function(panel, init.params, tranM, steadyM, alpha_mean, psi_mean, psi_sd, alpha_sd, csort, cnetw, csig, fsize, w_sigma, neduc, sort_gap, shocks) {
    new(
      "LaborMarket", 
      panel = panel,
      init.params = init.params,
      tranM = tranM,
      steadyM = steadyM,
      alpha_mean = alpha_mean,
      psi_mean = psi_mean,
      psi_sd = psi_sd,
      alpha_sd = alpha_sd,
      csort = csort,
      cnetw = cnetw,
      csig = csig,
      fsize = fsize,
      w_sigma = w_sigma,
      neduc = neduc,
      sort_gap = sort_gap,
      shocks = shocks
    )
  }


  # create the instance of the class
  lmarket = create_labor_m(data, list(nk = nk, nl = nl, nt = nt, ni = ni, ratiog = ratiog, lambda = lambda), G, H, alpha_mean, psi_mean, psi_sd, alpha_sd, csort, cnetw, csig, fsize, w_sigma, neduc, sort_gap, shocks)

  # output the class
  print(lmarket)
  invisible(lmarket)
}
