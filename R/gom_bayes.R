#' Bayesian Grade of Membership Mixture Model
#'
#' This function takes a data object and creates a joint posterior distribution of pure type probabilities, which
#' characterize a small set of extreme profiles, along with the posterior distribution of the grade of membership
#' estimates and the posterior distribution of the Dirichlet distribution hyperparameters.
#'
#' @param data A data frame with the categorical variables to be used be the model.
#' @param ntypes An integer indicating the number of pure type probabilities to be estimated.
#' @param alpha An array with the pure type probabilities. If specified, the model will not estimate
#' alpha and will only use the array provided by user.
#' @param burnin Number of iterations for the Markov chain achieve a stationary distribution.
#' @param ngibbs Number of iterations after the Markov chain achieved a stationary distribution.
#' @param omega The tuning parameter for the Metropolis–Hastings step.
#' @param eta The tuning parameter of the conditional Dirichlet distribution of Xi.
#' @param tau The shape parameter of the prior Gamma distribution for alpha.
#' @param beta The inverse scale parameter of the prior Gamma distribution for alpha.
#' @param gomscores Prefix for the gamma column names.
#'
#' @return A list with the posterior distributions of gamma, lambda, Xi, and a latent variable z.
#'
#' @importFrom  rlang .data
#' @export
gom_bayes <- function(data, ntypes = 2, alpha = "", burnin = 1000, ngibbs = 1000,
                      omega = 50, eta = 10, tau = 2, beta = 2, gomscores = "g"){
  data_dummy <- fastDummies::dummy_cols(data, select_columns = names(data), remove_selected_columns = TRUE)

  # Are known Dirichlet parameters supplied?

  ugom_gibbs_iter <- function(postg, ugom_X, ugom_alpha, ugom_L, ugom_g, gmeans, zmeans, Lambda, omega, eta, tau, beta){

    n = nrow(ugom_X)
    K = ncol(ugom_g)
    nvars = ncol(ugom_X)


    ugom_z = matrix(stats::runif(n*nvars), nrow = n, ncol = nvars)
    pr = matrix(0, nrow = 1, ncol = K)
    for(i in 1:n){
      for(j in 1:nvars){
        for(k in 1:K){
          if(ugom_X[i,j] != 0){
            pr[1, k] = ugom_g[i, k] * ugom_L[j,k]
          } else{
            pr[1, k] = ugom_g[i, k] * (1-ugom_L[j,k])
          }
        }
        pr = pr / sum(pr)
        z = 1
        sumpr = 0
        for(k in 1:K){
          sumpr = sumpr + pr[1, k]
          if(sumpr > ugom_z[i,j]){
            break
          }
          z = z + 1
        }
        ugom_z[i, j] = z
      }
    }

    for(j in 1:nvars){
      for(k in 1:K){
        a = 1 + sum((ugom_X[, j] != 0) * (ugom_z[, j] == k))
        b = 1 + sum((ugom_X[, j] == 0) * (ugom_z[, j] == k))
        ugom_L[j, k] = stats::qbeta(stats::runif(1), shape1 = a, shape2 = b) ### Lembrar de baixar pacote com a função pinvbeta
      }
    }
    newalpha = matrix(0, nrow = 1, ncol = K)
    for(i in 1:n){
      for(k in 1:K){
        newalpha[1, k] = ugom_alpha[k] + sum(ugom_z[i,] == k)
      }
      ugom_g[i, ] <- randomdirichlet(newalpha) ### Buscar função de distribuição dirichlet
    }

    #// Update alpha

    if(exists("dknown") == F) {
      ugom_alpha <- ugom_update_alpha(ugom_alpha, ugom_X, ugom_g, omega, eta, tau, beta)
    }

    #// Store L back as Lambda in Stata

    #st_matrix(st_local("Lambda"), ugom_L)

    #// Store a0 and xi back in Stata

    a0 = sum(ugom_alpha)
    e = ugom_alpha / a0

    #st_numscalar(st_local("a0"), a0)
    #st_matrix(st_local("Xi"), e)

    Xi <- e

    #// Process GoM-score variables, if need be

    #if(postg) {
    #  glist = tokens(st_local("gmeans"))
    #  for(k=1; k<=cols(glist); k++) {
    #    gvar = st_data(., glist[k], st_local("touse"))
    #    gvar = gvar + ugom_g[.,k]
    #    st_store(., glist[k], st_local("touse"), gvar)
    #  }
    #}

    if(postg){
      glist = gmeans
      llist = Lambda
      zlist = zmeans
      for(k in 1:ncol(glist)){
        gvar = glist[, k]
        gvar = gvar + ugom_g[, k]
        gmeans[ ,k] = gvar
        lvar = llist[, k]
        lvar = lvar + ugom_L[, k]
        Lambda[, k] = lvar
      }
      for(j in 1:nvars){
        zvar = zlist[, j]
        zvar = zvar + ugom_z[, j]
        zmeans[, j] = zvar
      }
      return(list("Lambda" = Lambda, "a0" = a0, "Xi" = Xi, "gmeans" = gmeans, "ugom_X" = ugom_X, "ugom_alpha" = ugom_alpha, "ugom_L" = ugom_L, "ugom_g" = ugom_g, "zmeans" = zmeans, "ugom_z" = ugom_z))
    }
    return(list("a0" = a0, "Xi" = Xi, "ugom_X" = ugom_X, "ugom_alpha" = ugom_alpha, "ugom_L" = ugom_L, "ugom_g" = ugom_g, "ugom_z" = ugom_z))
  }

  ugom_update_alpha <- function(ugom_alpha, ugom_X, ugom_g, omega, eta, tau, beta){
    a0 = sum(ugom_alpha)
    xi = ugom_alpha / a0

    K = length(ugom_alpha)
    N = nrow(ugom_X)

    # Get omega and eta, etc. from Stata

    #omega = strtoreal(st_local("omega"))
    #eta = strtoreal(st_local("eta"))
    #tau = strtoreal(st_local("tau"))
    #beta = strtoreal(st_local("beta"))

    # Draw candidate point for a0

    a0star = stats::qgamma(stats::runif(1), omega) * a0/omega

    # Calculate proposal ratio for a0

    ra0 = 0
    sgamma = lgamma(a0star) - lgamma(a0)
    for(k in 1:K) {
      ra0 = ra0 + xi[k]*sum(log(ugom_g[,k]))
      sgamma = sgamma + lgamma(xi[k]*a0) - lgamma(xi[k]*a0star)
    }
    ra0 = -(beta - ra0) * (a0star - a0) + N*sgamma
    ra0 = ra0 - omega*(a0/a0star - a0star/a0) +
      (tau-1)*log(a0star/a0) + (2*omega-1)*log(a0/a0star)
    ra0 = exp(ra0)

    # Update a0 if necessary

    change = 0
    if(stats::runif(1) < ra0) {
      a0 = a0star
      change = 1
    }

    # Draw candidate vector for xi

    xistar = randomdirichlet(eta*K*xi)

    # Calculate proposal ratio for xi

    re0 = 0
    sgamma = 0
    sg2 = 0
    for(k in 1:K) {
      re0 = re0 + (xistar[k]-xi[k]) * sum(log(ugom_g[,k]))
      sgamma = sgamma + lgamma(xi[k]*a0) - lgamma(xistar[k]*a0)
      sg2 = sg2 + lgamma(eta*K*xi[k]) - lgamma(eta*K*xistar[k]) +
        (xistar[k]-1)*log(xi[k]) - (xi[k]-1)*log(xistar[k])
    }
    re0 = a0*re0 + N*sgamma + sg2
    re0 = exp(re0)

    # Update if necessary

    if(stats::runif(1) < re0) {
      xi = xistar
      change = 1
    }

    if(change) {			# ugom_alpha requires updating
      if(is.null(xi*a0)){
        ugom_alpha= ugom_alpha
      } else{
        ugom_alpha = xi*a0
      }
    }
    return(ugom_alpha)
  }

  setup <- function(data, ntypes, alpha){
    ugom_X <- as.matrix(data)
    ugom_g <- matrix(1/ntypes, nrow = nrow(data), ncol = ntypes)
    if(exists("dknown")){
      ugom_alpha <- alpha
    } else{
      ugom_alpha <- matrix(0.25, nrow = 1, ncol = ntypes)
    }
    ugom_L <- matrix(stats::runif(ntypes*ncol(ugom_X)), nrow = ncol(ugom_X), ncol = ntypes)
    #ugom_L <- matrix(c(.9472316166,   .0522233748,
    #                   .9743182755,   .9457483679,
    #                   .1856478315,   .9487333737,
    #                   .8825376215,   .9440776079,
    #                   .0894258515,   .7505444902,
    #                   .9484983174,   .1121626508), nrow = 6, ncol = 2, byrow = T)

    return(list("ugom_X" = ugom_X, "ugom_g" = ugom_g, "ugom_alpha" = ugom_alpha, "ugom_L" = ugom_L))

  }

  randomdirichlet <- function(alpha){
    d <- vector(length = length(alpha))
    for(j in 1:length(alpha)){
      d[j] = stats::qgamma(shape = alpha[j], p=stats::runif(1))
      if(is.null(d[j])){
        d[j] = 1e-6
      }
    }
    return(d / sum(d))
  }

  if(alpha != ""){
    if(is.matrix(alpha)){
      if(nrow(alpha) != 1) {
        stop("`alpha` should only have 1 row")
      }
      if(ncol(alpha) != ntypes) {
        stop("`alpha` should have column dimension equal to the number of pure types")
      }
    }
    dknown = "dknown"
  }

  #Process GoM Scores
  tryCatch(gmeans <- as.data.frame(matrix(0, nrow = nrow(data_dummy), ncol = ntypes)))
  if(gomscores != ""){
    nsc <- length(gomscores)
    if(nsc < 1){
      stop("gomscores(): one stub name required")
    }
    for(i in 1:ntypes){
      #gmeans[[paste0(gomscores, "_", i)]] <- rep(0, nrow(data))
      names(gmeans)[i] <- paste0(gomscores, "_", i)
    }
  }

  #// Process output file name and options

  #      _prefix_saving `saving'
  #    local saving `"`s(filename)'"'
  #    local replace `"`s(replace)'"'
  #    local double `"`s(double)'"'
  #    local every `"`s(every)'"'
  #
  #    // Initialize Lambda matrix
  #
  #    tempname Lambda Xi
  #    local nvars : word count `varlist'
  #    matrix `Lambda' = J(`nvars',`ntypes',0)
  #    matrix `Xi'' = J(1, `ntypes', 0)
  nvars <- ncol(data_dummy)
  Lambda <- matrix(0, nrow = nvars, ncol = ntypes)
  zmeans <- matrix(0, nrow = nrow(data_dummy), ncol = nvars)
  Xi <- matrix(0, nrow = 1, ncol = ntypes)

  #// Form expression list for posting
  #
  #    for(k in 1:ntypes){
  #      for(j in 1:nvars){
  #
  #      }
  #    }
  #      forvalues k = 1/`ntypes' {
  #        forvalues j = 1/`nvars' {
  #	 	local explist `explist' (`Lambda'[`j',`k'])
  #      }
  #  }
  #  tempname a0
  #  local explist `explist' (scalar(`a0'))
  #forvalues k = 1/`ntypes' {
  #  local explist `explist' (`Xi'[1,`k'])
  #}

  #  // Initialize mata matrices

  temp <- setup(data_dummy, ntypes, alpha)
  ugom_X <- temp$ugom_X
  ugom_g <- temp$ugom_g
  ugom_alpha <- temp$ugom_alpha
  ugom_L <- temp$ugom_L
  a0 <- vector(length = burnin+ngibbs)
  alphas <- matrix(0, nrow = burnin+ngibbs, ncol = ntypes)
  g_dist <- list()
  z_dist <- list()
  lambda_dist <- list()
  cat("MCMC Iteration: \n")
  pb = utils::txtProgressBar(min = 1, max = burnin+ngibbs, initial = 1, style = 3, width = 60)
  for(i in 1:(burnin+ngibbs)){

    if(i <= burnin){
      temp <- ugom_gibbs_iter(postg = 0, ugom_X = ugom_X, ugom_alpha = ugom_alpha, ugom_L = ugom_L, ugom_g = ugom_g, gmeans = gmeans, zmeans = zmeans, Lambda = Lambda, omega = omega, eta = eta, tau = tau, beta = beta)
      Xi <- temp$Xi
      a0[i] <- temp$a0
      alphas[i,] <- temp$ugom_alpha
      ugom_X <- temp$ugom_X
      ugom_g <- temp$ugom_g
      g_dist[[i]] <- ugom_g
      ugom_alpha <- temp$ugom_alpha
      ugom_L <- temp$ugom_L
      lambda_dist[[i]] <- ugom_L
      z_dist[[i]] <- temp$ugom_z
    } else{
      temp <- ugom_gibbs_iter(postg = 1, ugom_X = ugom_X, ugom_alpha = ugom_alpha, ugom_L = ugom_L, ugom_g = ugom_g, gmeans = gmeans, zmeans = zmeans, Lambda = Lambda, omega = omega, eta = eta, tau = tau, beta = beta)
      Lambda <- temp$Lambda
      Xi <- temp$Xi
      a0[i] <- temp$a0
      alphas[i,] <- temp$ugom_alpha
      ugom_alpha <- temp$ugom_alpha
      gmeans <- temp$gmeans
      zmeans <- temp$zmeans
      ugom_X <- temp$ugom_X
      ugom_g <- temp$ugom_g
      g_dist[[i]] <- ugom_g
      ugom_L <- temp$ugom_L
      z_dist[[i]] <- temp$ugom_z
      lambda_dist[[i]] <- ugom_L
    }
    utils::setTxtProgressBar(pb,i)
  }


  #	  // Normalize GoM score variables
  #
  #	  local w : word count `gmeans'
  #	forvalues i = 1/`w' {
  #	  local gvar : word `i' of `gmeans'
  #		qui replace `gvar' = `gvar'/`ngibbs' if `touse'
  #	}
  #	  end

  #	  mata:

  lmeans <- Lambda/ngibbs
  n <- sapply(data, dplyr::n_distinct)
  lmeans <- as.data.frame(lmeans)
  lmeans$groups = rep(names(n), times = n)
  lmeans <- lmeans %>%
    dplyr::group_by(.data$groups) %>%
    dplyr::transmute(K1 = .data$V1/sum(.data$V1),
                     K2 = .data$V2/sum(.data$V2)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$groups) %>%
    round(4)
  names <- vector()
  for(i in 1:length(n)){
    names <- c(names, paste0(names(n[i]), "_", 1:n[i]))
  }
  lmeans$names <- names

  lmeans <- tibble::column_to_rownames(lmeans, "names")
  lmeans$prop <- sapply(data_dummy, mean)
  lmeans$prop <- round(lmeans$prop, 4)

  lmeans <- lmeans %>% dplyr::mutate(lmfr1 = round(.data$K1/.data$prop,4),
                                     lmfr2 = round(.data$K2/.data$prop, 4))


  return(list("gmeans" = gmeans/ngibbs, "lmeans" = lmeans, "zmeans" = zmeans/ngibbs, "a0" = a0, "Xi" = Xi, "alphas" = alphas, "Lambda_dist" = lambda_dist, "Gamma_dist" = g_dist, "z_dist" = z_dist))
}
#'
#' @examples





