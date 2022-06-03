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
#' @return An object of class \emph{gom_bayes} with the posterior distributions of gamma, lambda, Xi, and alpha.
#'
#' @importFrom  rlang .data
#' @export
#' @examples
#'
#' data <- data.frame(x1 = round(stats::runif(n = 50, 1, 2), 0),
#'                    x2 = round(stats::runif(n = 50, 1, 3), 0),
#'                    x3 = round(stats::runif(n = 50, 1, 4), 0))
#'
#' gom::gom_bayes(data, ntypes = 2, ngibbs = 250, burnin = 250)
#'
gom_bayes <- function(data, ntypes = 2, alpha = "", burnin = 250, ngibbs = 250,
                      omega = 50, eta = 10, tau = 2, beta = 2, gomscores = "g"){

  data <- dplyr::mutate(data, dplyr::across(.fns = as.factor))

  data_dummy <- fastDummies::dummy_cols(data, select_columns = names(data), remove_selected_columns = TRUE)

  # Are known Dirichlet parameters supplied?

  ugom_gibbs_iter <- function(postg, ugom_X, ugom_alpha, ugom_L, ugom_g, gmeans, Lambda, omega, eta, tau, beta){

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
        ugom_L[j, k] = stats::qbeta(stats::runif(1), shape1 = a, shape2 = b)
      }
      ugom_L[j, ] <- ugom_L[j, ]/sum(ugom_L[j, ])
    }

    newalpha = matrix(0, nrow = 1, ncol = K)
    for(i in 1:n){
      for(k in 1:K){
        newalpha[1, k] = ugom_alpha[k] + sum(ugom_z[i,] == k)
      }
      ugom_g[i, ] <- randomdirichlet(newalpha) # Buscar função de distribuição dirichlet
    }

    # Update alpha

    if(exists("dknown") == F) {
      ugom_alpha <- ugom_update_alpha(ugom_alpha, ugom_X, ugom_g, omega, eta, tau, beta)
    }

    a0 = sum(ugom_alpha)
    Xi <- ugom_alpha / a0

    if(postg){
      for(k in 1:ncol(gmeans)){
        gvar = gmeans[, k]
        gvar = gvar + ugom_g[, k]
        gmeans[ ,k] = gvar
        lvar = Lambda[, k]
        lvar = lvar + ugom_L[, k]
        Lambda[, k] = lvar
      }

      return(list("Lambda" = Lambda, "a0" = a0, "Xi" = Xi, "gmeans" = gmeans, "ugom_alpha" = ugom_alpha, "ugom_L" = ugom_L, "ugom_g" = ugom_g))
    }
    return(list("a0" = a0, "Xi" = Xi, "ugom_alpha" = ugom_alpha, "ugom_L" = ugom_L, "ugom_g" = ugom_g))
  }

  ugom_update_alpha <- function(ugom_alpha, ugom_X, ugom_g, omega, eta, tau, beta){
    a0 = sum(ugom_alpha)
    xi = ugom_alpha / a0

    K = length(ugom_alpha)
    N = nrow(ugom_X)

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

    if(is.nan(re0) && stats::runif(1) < re0) {
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

  ugom_X <- as.matrix(data_dummy)

  ugom_g <- matrix(1/ntypes, nrow = nrow(data_dummy), ncol = ntypes)
  if(exists("dknown")){
    ugom_alpha <- alpha
  } else{
    ugom_alpha <- matrix(0.25, nrow = 1, ncol = ntypes)
  }

  ugom_L <- matrix(stats::runif(ntypes*ncol(ugom_X)), nrow = ncol(ugom_X), ncol = ntypes)

  randomdirichlet <- function(alpha){
    d <- vector(length = length(alpha))
    for(j in 1:length(alpha)){
      d[j] = stats::qgamma(shape = alpha[j], p=stats::runif(1))
      if(is.nan(d[j])){
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
      names(gmeans)[i] <- paste0(gomscores, "_", i)
    }
  }


  nvars <- ncol(data_dummy)
  Lambda <- matrix(0, nrow = nvars, ncol = ntypes)
  Xi <- matrix(0, nrow = 1, ncol = ntypes)

  a0 <- vector(length = burnin+ngibbs)
  xis <- matrix(0, nrow = burnin+ngibbs, ncol = ntypes)
  alphas <- matrix(0, nrow = burnin+ngibbs, ncol = ntypes)
  g_dist <- list()
  lambda_dist <- list()
  cat("MCMC Iteration: \n")
  pb = utils::txtProgressBar(min = 0, max = burnin+ngibbs, initial = 0, style = 3, width = 60)
  for(i in 1:(burnin+ngibbs)){

    if(i <= burnin){
      temp <- ugom_gibbs_iter(postg = 0, ugom_X = ugom_X, ugom_alpha = ugom_alpha, ugom_L = ugom_L, ugom_g = ugom_g, gmeans = gmeans, Lambda = Lambda, omega = omega, eta = eta, tau = tau, beta = beta)
      Xi <- temp$Xi
      a0[i] <- temp$a0
      xis[i, ] <- temp$Xi
      alphas[i,] <- temp$ugom_alpha
      ugom_g <- temp$ugom_g
      g_dist[[i]] <- ugom_g
      ugom_alpha <- temp$ugom_alpha
      ugom_L <- temp$ugom_L
      lambda_dist[[i]] <- ugom_L
    } else{
      temp <- ugom_gibbs_iter(postg = 1, ugom_X = ugom_X, ugom_alpha = ugom_alpha, ugom_L = ugom_L, ugom_g = ugom_g, gmeans = gmeans, Lambda = Lambda, omega = omega, eta = eta, tau = tau, beta = beta)
      Lambda <- temp$Lambda
      Xi <- temp$Xi
      a0[i] <- temp$a0
      xis[i, ] <- temp$Xi
      alphas[i,] <- temp$ugom_alpha
      ugom_alpha <- temp$ugom_alpha
      gmeans <- temp$gmeans
      ugom_g <- temp$ugom_g
      g_dist[[i]] <- ugom_g
      ugom_L <- temp$ugom_L
      lambda_dist[[i]] <- ugom_L
    }
    utils::setTxtProgressBar(pb,i)
  }
  cat("\n")

  lmeans <- Lambda/ngibbs
  n <- sapply(data, dplyr::n_distinct)
  lmeans <- as.data.frame(lmeans)
  lmeans$groups = rep(names(n), times = n)
  lmeans <- lmeans %>%
    dplyr::group_by(.data$groups) %>%
    dplyr::transmute(n = dplyr::n(),
                     K1 = .data$V1/sum(.data$V1),
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
  lmeans <- lmeans %>% dplyr::select(.data$prop, n, .data$K1, .data$K2, .data$lmfr1, .data$lmfr2)


  lambdas <- as.data.frame(t(sapply(lambda_dist, `[`, 1:nvars, 1:ntypes)))
  lambdas <- round(lambdas, 5)
  names(lambdas) <- paste0("k", rep(1:ntypes, each = nvars),
                           "j_", rep(rep(1:ncol(data), unlist(lapply(as.data.frame(sapply(data, levels)), length))), times = ntypes),
                           "_l", rep(unlist(lapply(sapply(data, levels), sort)), times = ntypes))

  gammas <- as.data.frame(t(sapply(g_dist, `[`, 1:nrow(ugom_X), 1:ntypes)))
  gammas <- round(gammas, 5)
  names(gammas) <- paste0("i",rep(1:nrow(ugom_X), times = ntypes),"_k", rep(1:ntypes, each = nrow(ugom_X)))


  output <- list("gmeans" = gmeans/ngibbs, "lmeans" = lmeans, "a0" = a0, "Xi" = xis, "alphas" = alphas, "Lambda_dist" = lambdas, "Gamma_dist" = gammas)

  class(output) <- "gom_bayes"
  return(output)
}

