#' @title Print for \code{gom_ml} objects
#'
#' @description Print method for an objects of class \emph{gom_ml}.
#'
#' @param x An object with class \emph{gom_ml}.
#' @param k Number indicating the pure type result to be printed.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Prints the descriptive table inside the gom_ml object along with the likelihood and
#' AIC for the desired pure type set. If k is NULL (default), it chooses the model with the smallest AIC.
#'
#' @return No return value.
#'
#' @export
#' @examples
#'
#' data <- data.frame(x1 = round(stats::runif(n = 50, 1, 2), 0),
#'                    x2 = round(stats::runif(n = 50, 1, 2), 0),
#'                    x3 = round(stats::runif(n = 50, 1, 2), 0),
#'                    Id = 1:50)
#'
#' model <- gom::gom_ml(data.object = data, case.id = "Id", initial.lambda = "random", MC_iter = 150)
#'
#' print(model)
#'
print.gom_ml <- function(x, k = NULL, ...){
  if(!is.null(k) && is.numeric(k)){
    cat("Grade of Membership Mixture Model by Maximum Likelihood \n\n")

    cat("Descriptive table with the marginal distribution and lambda estimates \n\n")
    writeLines(x[[k]]$Table)
    cat("\n")
    cat("Likelihood: \t", x[[k]]$Likelihood, "\n")
    cat("AIC: \t\t", x[[k]]$AIC, "\n")
  } else{
    profile <- names(which.min(sapply(x, `[[`, 4)))
    cat("Grade of Membership Mixture Model by Maximum Likelihood \n\n")

    cat(paste("Model results for the smallest AIC:\n\n"))
    cat("Descriptive table with the marginal distribution and lambda estimates \n\n")
    writeLines(x[[profile]]$Table)
    cat("\n")
    cat("Likelihood: \t", x[[profile]]$Likelihood, "\n")
    cat("AIC: \t\t", x[[profile]]$AIC, "\n")
  }
}
