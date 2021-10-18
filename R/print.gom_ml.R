#' @title Print for \code{gom_ml} objects
#'
#' @description Print objects with class \emph{gom_ml}.
#'
#' @param x An object with class \emph{gom_ml}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
#' @examples
#' \dontrun{
#' data <- data.frame(x1 = round(stats::runif(n = 500, 1, 2), 0),
#'                    x2 = round(stats::runif(n = 500, 1, 3), 0),
#'                    x3 = round(stats::runif(n = 500, 1, 4), 0),
#'                    x4 = round(stats::runif(n = 500, 1, 5), 0),
#'                    Id = 1:500)
#'
#' model <- gom_ml(data.object = data, case.id = "Id", initial.lambda = "random")
#'
#' print(model)
#' }
print.gom_ml <- function(x, ...){
  profile <- names(which.min(sapply(x, `[[`, 4)))
  cat("Grade of Membership Mixture Model by Maximum Likelihood \n\n")

  cat(paste("Results for model with least AIC:\n\n"))
  cat("Lambda Descriptive Table \n\n")
  writeLines(x[[profile]]$Table)
  cat("\n")
  cat("Likelihood: \t", x[[profile]]$Likelihood, "\n")
  cat("AIC: \t\t", x[[profile]]$AIC, "\n")
}
