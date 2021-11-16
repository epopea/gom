#' @title Print for \code{gom_bayes} objects
#'
#' @description Print method for an objects of class \emph{gom_bayes}.
#'
#' @param x An object with class \emph{gom_bayes}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Prints a descriptive table with the posterior lambdas organized by variables and their categories.
#'
#' @return No return value.
#'
#' @export
#' @examples
#'
#' data <- data.frame(x1 = round(stats::runif(n = 50, 1, 2), 0),
#'                    x2 = round(stats::runif(n = 50, 1, 3), 0),
#'                    x3 = round(stats::runif(n = 50, 1, 4), 0))
#'
#' model <- gom::gom_bayes(data, ntypes = 2, ngibbs = 250, burnin = 250)
#'
#' print(model)
#'
print.gom_bayes <- function(x, ...){
  cat("Bayesian Grade of Membership Mixture Model \n\n")
  cat("Descriptive table with the marginal distribution and lambda estimates \n\n")
  print(x$lmeans)
  cat("\n")
}
