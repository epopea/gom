#' Grade of Membership Mixture Model
#'
#' This function takes a data object and creates a fixed, small set of pure types, each characterized
#' by the probabilities of answering a category of each question, along with a set of scores (or grades
#' of membership) that measure each subject’s degree of similarity to each pure type.
#'
#'
#' @param data.object A data frame of data or an object that can be coerced to a data frame.
#' @param initial.K An integer indicating the initial number of pure type probabilities to be estimated.
#' @param final.K An integer indicating the last number of pure type probabilities to be estimated.
#' @param gamma.algorithm A character that specifies the algorithm to be used to estimate the gamma values.
#' @param initial.gamma A character that specifies how the initial gamma values should be specified. If "gamma.object", then initial.gamma.object is needed.
#' @param initial.gamma.object A data frame with the initial gamma values for each observation.
#' @param gamma.fit A logical value indicating if the gamma values need to be estimated.
#' @param lambda.algorithm A character that specifies the algorithm to be used to estimate the lambda values.
#' @param initial.lambda A character that specifies how the initial lambda values should be specified. If "lambda.object", then initial.lambda.object is needed.
#' @param initial.lambda.object An array with the initial lambda values for each category of each variable.
#' @param lambda.fit A logical value indicating if the lambda values need to be estimated.
#' @param case.id A character with the name of the variable with the ID of each observation.
#' @param case.weight A character with the name of the variable with the weight of each observation.
#' @param internal.var A character vector with the name of the variables to be used by the model.
#' @param order.K A logical value indicating if the lambdas should be organized with 99 percent confidence interval.
#' @param omega.fit A logical value specifying if the model should generate missing patterns of the data before compiling.
#' @param dec.char A character with the decimal symbol to be used in the output files.
#' @param MC_iter The number of iterations to be used within the Monte Carlo Estimation.
#'
#' @return gom_ml saves two files, one containing a table with unique data used by the model with the grades of membership, and another
#' with the initial and final pure type probabilities and some summary statistics. The gom_ml function also returns a list with information about
#' the results for each number of pure types specified by the user, each of them with following components:
#'
#' \describe{
#'   \item{Gik}{A data frame with the grades of membership of each ordered unique pattern of the given data.}
#'   \item{Pkjl}{An array with dimensions equal to the number of pure types, variables and categories, respectively, with the pure type probabilities.}
#'   \item{Likelihood}{The maximum log-likelihood achieved by the model.}
#'   \item{AIC}{The Akaike Information Criterion.}
#'   \item{Table}{A table with the lambda results organized by variables and their categories.}
#' }
#' @export
gom_ml <- function (data.object = NULL,
                     initial.K = 2, final.K = initial.K,
                     gamma.algorithm = c("gradient.1992", "woodbury.1974"),
                     initial.gamma = c("equal.values", "random", "pure1", "gamma.object"),
                     initial.gamma.object = NULL,
                     gamma.fit = TRUE,
                     lambda.algorithm = c("gradient.1992", "woodbury.1974"),
                     initial.lambda = c("random", "pure1", "equal.values", "lambda.object"),
                     initial.lambda.object = NULL,
                     lambda.fit = TRUE,
                     case.id = NA,
                     case.weight = NA,
                     internal.var = NULL,
                     order.K = TRUE,
                     omega.fit = FALSE,
                     dec.char = ".",
                     MC_iter = 1000) {

  data.object <- as.data.frame(data.object)

  GoM <- '
    using namespace std;
    using namespace Rcpp;

    const int MAXITER_MODEL = 500;
    const int MAXITER_PARAMETERS = 25;

    const double ZERO = 1.0E-20;
    const double REALBIG = 1.0E+30;
    const double CTOL = 1.0E-07;

    IntegerVector   baselevel       (baselevel_);
    IntegerVector   ljlevels        (ljlevels_);
    NumericMatrix   cell            (cell_);
    CharacterVector gammaalgorithm  (gammaalgorithm_);
    NumericMatrix   G               (FG_);
    LogicalVector   gammafit        (gammafit_);
    CharacterVector lambdaalgorithm (lambdaalgorithm_);
    NumericVector   P               (FP_);
    LogicalVector   lambdafit       (lambdafit_);
    IntegerMatrix   FITP            (FITP_);

    int I = (cell.nrow() - 1);
    int J = (cell.ncol() - 2);
    int K = (G.ncol() - 1);

    int i, k, k_k, j, l, iter, subiter;
    double ploglik, difflik;
    double sumlik, p_ijl, part, partlik, newpartlik, sumG, sumP, startlik, curlik;

    double old_G[(K + 1)];
    double new_G[(K + 1)];

    NumericVector old_P(clone(FP_));
    NumericVector new_P(clone(FP_));

    vector<double> loglik(2);
    char buffer[255];

    // ## WOODBURY VARIABLES ####################

    double numer, denom;

    // ##########################################

    // ## GRADIENT VARIABLES ####################

    const int HALFSTEPS = MAXITER_PARAMETERS;

    const double ITOL = 0.0005;

    const double MAXSTEP = 1.0;
    const double MINSTEP = ZERO;

    int l_l, converged, halfstep, somefree;

    double f0, f1, norm, g_ij, bestlik, stepsize;

    int    freeG[(K + 1)];
    double dL_dG[(K + 1)];

    NumericVector dL_dP(clone(FP_));
    NumericVector freeP(clone(FP_));

    // ##########################################

    loglik[0] = 0.0;
    ploglik = loglik[0];
    // ## BEGIN loglikelihood FUNCTION ## //
    sumlik = 0.0;
    for (i = 1; i <= I; i++) {
      partlik = 0.0;
      for (j = 1; j <= J; j++) {
        l = cell(i, j) + 1 - baselevel(j);
        p_ijl = 0.0;
        for (k = 1; k <= K; k++) {
          p_ijl += (G(i, k) * P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)));
        }
        if (p_ijl < ZERO) {
          p_ijl = ZERO;
        }
        if (p_ijl < REALBIG) {
          partlik += log(p_ijl);
        }
      }
      sumlik += (cell(i, (J + 1)) * partlik);
    }
    // ## END loglikelihood FUNCTION ## //
    loglik[0] = sumlik;
    sprintf(buffer, "%-10.4f", loglik[0]);
    Rcout << "Primal Log-Likelihood is:\t" << buffer << endl << endl;
    loglik[1] = loglik[0];
    for (iter = 0; iter < MAXITER_MODEL; iter++) {
     difflik = (loglik[1] - ploglik);
     ploglik = loglik[1];
     if (iter) {
      if (fabs(difflik / loglik[1]) < CTOL) {
       break;
      }
     }
     if (gammafit[0] && gammaalgorithm[0] == "woodbury.1974") {
      // ## BEGIN Fit_G_Woodbury_1974 FUNCTION ## //
      subiter = 0;
      for (i = 1; i <= I; i++) {
       for (k = 1; k <= K; k++) {
        old_G[k] = G(i, k);
       }
       // ## BEGIN partlikelihood FUNCTION ## //
       part = 0.0;
       for (j = 1; j <= J; j++) {
        l = cell(i, j) + 1 - baselevel[j];
        p_ijl = 0.0;
        for (k = 1; k <= K; k++) {
         p_ijl += (G(i, k) * P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)));
        }
        if (p_ijl < ZERO) {
         p_ijl = ZERO;
        }
        if (p_ijl < REALBIG) {
         part += log(p_ijl);
        }
       }
       // ## END partlikelihood FUNCTION ## //
       partlik = part;
       for (subiter = 0; subiter < MAXITER_PARAMETERS; subiter++) {
        for (k = 1; k <= K; k++) {
         new_G[k] = 0.0;
         for (j = 1; j <= J; j++) {
          l = cell(i, j) + 1 - baselevel[j];
          p_ijl = 0.0;
          for (k_k = 1; k_k <= K; k_k++) {
           p_ijl += (G(i, k_k) * P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k_k)));
          }
          if (p_ijl > ZERO) {
           new_G[k] += ((G(i, k) * P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k))) / p_ijl);
          }
         }
        }
        for (k = 1; k <= K; k++) {
         G(i, k) = (new_G[k] / J);
        }
        // ## BEGIN rescale_G FUNCTION ## //
        sumG = 0.0;
        for (k = 1; k <= K; k++) {
         if (G(i, k) < ZERO) {
          G(i, k) = 0.0;
         }
         sumG += G(i, k);
        }
        if (sumG < ZERO) {
         sumG = (double)K;
          for (k = 1; k <= K; k++) {
           G(i, k) = 1.0;
          }
        }
        for (k = 1; k <= K; k++) {
         G(i, k) = (G(i, k) / sumG);
        }
        // ## END rescale_G FUNCTION ## //
        // ## BEGIN partlikelihood FUNCTION ## //
        part = 0.0;
        for (j = 1; j <= J; j++) {
         l = cell(i, j) + 1 - baselevel[j];
         p_ijl = 0.0;
         for (k = 1; k <= K; k++) {
          p_ijl += (G(i, k) * P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)));
         }
         if (p_ijl < ZERO) {
          p_ijl = ZERO;
         }
         if (p_ijl < REALBIG) {
          part += log(p_ijl);
         }
        }
        // ## END partlikelihood FUNCTION ## //
        newpartlik = part;
        if (newpartlik < partlik) {
         for (k = 1; k <= K; k++) {
          G(i, k) = old_G[k];
         }
         break;
        } else {
         if (partlik > ZERO) {
          if ((fabs(newpartlik - partlik) / partlik) < CTOL) {
           break;
          }
         }
        }
       }
      }
      sprintf(buffer, "%03d.%03d%s%015.5f%s%015.7f", iter, subiter, "\t", loglik[1], "\t", fabs((loglik[1] - ploglik) / ploglik));
      Rcout << "Fit G (Woodbury 1974):\t" << buffer << endl << endl;
      // ## END Fit_G_Woodbury_1974 FUNCTION ## //
     } else if (gammafit[0] && gammaalgorithm[0] == "gradient.1992") {
      // ## BEGIN Fit_G_Gradient_1992 FUNCTION ## //
      subiter = 0;
      startlik = loglik[1];
      for (i = 1; i <= I; i++) {
       converged = 0;
       for (subiter = 0; converged == 0; subiter++) { //##//
        // ## BEGIN cellgradient_G FUNCTION ## //
        norm = 0.0;
        for (k = 1; k <= K; k++) {
         dL_dG[k] =- (double)J;
         for (j = 1; j <= J; j++) {
          l = cell(i, j) + 1 - baselevel(j);
          p_ijl = 0.0;
          for (k_k = 1; k_k <= K; k_k++) {
           p_ijl = p_ijl + (G(i, k_k) * P(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k_k));
          }
          if (p_ijl > ZERO) {
           dL_dG[k] += (P(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k) / p_ijl);
          }
         }
         if ((G(i, k) <= ZERO && dL_dG[k] < 0.0) || ((1.0 - G(i, k)) <= ZERO && dL_dG[k] > 0.0)) {
          freeG[k] = 0;
         } else {
          freeG[k] = 1;
          norm += (dL_dG[k] * dL_dG[k]);
         }
        }
        // ## END cellgradient_G FUNCTION ## //
        if (norm <= ZERO) {
         converged = 1;
         break;
        }
        somefree = 0;
        for (k = 1; k <= K; k++) {
         old_G[k] = G(i, k);
         new_G[k] = G(i, k);
         if (freeG[k]) {
          somefree++;
         }
        }
        if (!somefree) {
         break;
        }
        // ## BEGIN partlikelihood FUNCTION ## //
        part = 0.0;
        for (j = 1; j <= J; j++) {
         l = cell(i, j) + 1 - baselevel[j];
         p_ijl = 0.0;
         for (k = 1; k <= K; k++) {
          p_ijl += (G(i, k) * P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)));
         }
         if (p_ijl < ZERO) {
          p_ijl = ZERO;
         }
         if (p_ijl < REALBIG) {
          part += log(p_ijl);
         }
        }
        // ## END partlikelihood FUNCTION ## //
        f0 = part;
        stepsize = MAXSTEP;
        bestlik = f0;
        for (halfstep = 0; halfstep < HALFSTEPS && stepsize > MINSTEP; halfstep++) {
         for (k = 1;k <= K; k++) {
          G(i, k) = old_G[k] + stepsize * dL_dG[k];
         }
         // ## BEGIN rescale_G FUNCTION ## //
         sumG = 0.0;
         for (k = 1; k <= K; k++) {
          if (G(i, k) < ZERO) {
           G(i, k) = 0.0;
          }
          sumG += G(i, k);
         }
         if (sumG < ZERO) {
          sumG = (double)K;
          for (k = 1; k <= K; k++) {
           G(i, k) = 1.0;
          }
         }
         for (k = 1; k <= K; k++) {
          G(i, k) = (G(i, k) / sumG);
         }
         // ## END rescale_G FUNCTION ## //
         // ## BEGIN partlikelihood FUNCTION ## //
         part = 0.0;
         for (j = 1; j <= J; j++) {
          l = cell(i, j) + 1 - baselevel[j];
          p_ijl = 0.0;
          for (k = 1; k <= K; k++) {
           p_ijl += (G(i, k) * P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)));
          }
          if (p_ijl < ZERO) {
           p_ijl = ZERO;
          }
          if (p_ijl < REALBIG) {
           part += log(p_ijl);
          }
         }
         // ## END partlikelihood FUNCTION ## //
         f1 = part;
         if (f1 > bestlik) {
          bestlik = f1;
          for (k = 1; k <= K; k++) {
           new_G[k] = G(i, k);
          }
         } else {
          if (bestlik > (f0 + ZERO) && bestlik > (f1 + ZERO)) {
           break;
          }
         }
         stepsize = (stepsize / 2.0);
        }
        f1 = bestlik;
        for (k = 1;k <= K; k++) {
         G(i, k) = new_G[k];
        }
        if (fabs(f0) > ZERO && fabs((f1 - f0) / f0) < CTOL) {
         converged = 1;
        }
        if (fabs(f0) > ZERO && fabs((f1-f0)/f0) < ITOL && subiter > MAXITER_PARAMETERS) {
         break;
        }
       }
      }
      sprintf(buffer, "%03d.%03d%s%015.5f%s%015.7f", iter, subiter, "\t", loglik[1], "\t", fabs((loglik[1] - ploglik) / ploglik));
      Rcout << "Fit G (Gradient 1992):\t" << buffer << endl << endl;
      // ## END Fit_G_Gradient_1992 FUNCTION ## //
     }
     // ## BEGIN loglikelihood FUNCTION ## //
     sumlik = 0.0;
     for (i = 1; i <= I; i++) {
      partlik = 0.0;
      for (j = 1; j <= J; j++) {
       l = cell(i, j) + 1 - baselevel(j);
       p_ijl = 0.0;
       for (k = 1; k <= K; k++) {
        p_ijl += (G(i, k) * P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)));
       }
       if (p_ijl < ZERO) {
        p_ijl = ZERO;
       }
       if (p_ijl < REALBIG) {
        partlik += log(p_ijl);
       }
      }
      sumlik += (cell(i, (J + 1)) * partlik);
     }
     // ## END loglikelihood FUNCTION ## //
     loglik[1] = sumlik;
     if (lambdafit[0] && lambdaalgorithm[0] == "woodbury.1974") {
     // ## BEGIN Fit_P_Woodbury_1974 FUNCTION ## //
     startlik = loglik[1];
     for (subiter = 0; subiter < MAXITER_PARAMETERS; subiter++) {
      for (k = 1; k <= K; k++) {
       for (j = 1; j <= J; j++) {
        for (l = 1; l <= ljlevels(j); l++) {
         old_P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) = P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k));
        }
       }
      }
      for (k = 1; k <= K; k++) {
       for (j = 1; j <= J; j++) {
        if (FITP(k, j)) {
         for (l = 1; l <= ljlevels(j); l++) {
          numer = 0.0;
          denom = 0.0;
          for (i = 1; i <= I; i++) {
           p_ijl = 0.0;
           for (k_k = 1; k_k <= K; k_k++) {
            p_ijl += (G(i, k_k) * P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k_k)));
           }
           if (p_ijl > ZERO) {
            g_ij = ((G(i, k) * P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k))) / p_ijl);
           } else {
            g_ij = 0.0;
           }
           denom += g_ij;
           if (l == (cell(i, j) + 1 - baselevel(j))) {
            numer += g_ij;
           }
          }
          if (denom > ZERO) {
           new_P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) = (numer / denom);
          } else {
           new_P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) = 0.0;
          }
         }
        }
       }
      }
      for (k = 1; k <= K; k++) {
       for (j = 1; j <= J; j++) {
        if (FITP(k, j)) {
         for (l = 1; l <= ljlevels(j); l++) {
          P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) = new_P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k));
         }
         // ## BEGIN rescale_P FUNCTION ## //
         sumP = 0.0;
         for (l = 1; l <= ljlevels(j); l++) {
          if (P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) < ZERO) {
           P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) = 0.0;
          }
          sumP += P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k));
         }
         if (sumP < ZERO) {
          sumP = (double)ljlevels(j);
          for (l = 1; l <= ljlevels(j); l++) {
           P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) = 1.0;
          }
         }
         for (l = 1; l <= ljlevels(j); l++) {
          P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) = (P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) / sumP);
         }
         // ## END rescale_P FUNCTION ## //
        }
       }
      }
      // ## BEGIN loglikelihood FUNCTION ## //
      sumlik = 0.0;
      for (i = 1; i <= I; i++) {
       partlik = 0.0;
       for (j = 1; j <= J; j++) {
        l = cell(i, j) + 1 - baselevel(j);
        p_ijl = 0.0;
        for (k = 1; k <= K; k++) {
         p_ijl += (G(i, k) * P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)));
        }
        if (p_ijl < ZERO) {
         p_ijl = ZERO;
        }
        if (p_ijl < REALBIG) {
         partlik += log(p_ijl);
        }
       }
       sumlik += (cell(i, (J + 1)) * partlik);
      }
      // ## END loglikelihood FUNCTION ## //
      curlik = sumlik;
      sprintf(buffer, "%03d.%03d%s%015.5f%s%015.7f", iter, subiter, "\t", startlik, "\t", fabs((startlik - ploglik) / ploglik));
      Rcout << "Fit P (Woodbury 1974):\t" << buffer << endl;
      if (curlik < startlik) {
       for (k = 1; k <= K; k++) {
        for (j = 1; j <= J; j++) {
         for (l = 1; l <= ljlevels(j); l++) {
          P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) = old_P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k));
         }
        }
       }
       loglik[1] = startlik;
       break;
      } else {
       if ((startlik > ZERO) && (fabs(curlik - startlik) / startlik) < CTOL) {
        break;
       }
      }
      startlik = curlik;
     }
     Rcout << endl;
     // ## END Fit_P_Woodbury_1974 FUNCTION ## //
     } else if (lambdafit[0] && lambdaalgorithm[0] == "gradient.1992") {
      // ## BEGIN Fit_P_Gradient_1992 FUNCTION ## //
      converged = 0;
      for (subiter = 0; converged == 0; subiter++) { //##//
       // ## BEGIN gradient_P FUNCTION ## //
       norm = 0.0;
       for (k = 1; k <= K; k++) {
        for (j = 1; j <= J; j++) {
         for (l = 1; l <= ljlevels(j); l++) {
          dL_dP(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k) = 0.0;
          for (i = 1; i <= I; i++) {
           l_l = cell(i, j) + 1 - baselevel(j);
           p_ijl = 0.0;
           for (k_k = 1; k_k <= K; k_k++) {
            p_ijl += (G(i, k_k) * P(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k_k));
           }
           if (l == l_l) {
            if (p_ijl > ZERO) {
             dL_dP(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k) += (cell(i, (J + 1)) * G(i, k) * ((1.0 / p_ijl) - 1.0));
            }
           } else {
            dL_dP(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k) -= (cell(i, (J + 1)) * G(i, k));
           }
          }
          if (((P(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k) <= ZERO) && (dL_dP(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k) < 0.0)) || (((1.0 - P(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) <= ZERO) && (dL_dP(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k) > 0.0))) {
           freeP(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k) = 0;
          } else {
           freeP(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k) = 1;
           norm += dL_dP(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k) * dL_dP(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k);
          }
         }
        }
       }
       // ## END gradient_P FUNCTION ## //
       if (norm <= ZERO) {
        converged = 1;
        break;
       }
       somefree = 0;
       for (k = 1; k <= K; k++) {
        for (j = 1; j <= J; j++) {
         for (l = 1; l <= ljlevels(j); l++) {
          old_P(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k) = P(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k);
          new_P(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k) = P(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k);
          if (freeP(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) {
           somefree++;
          }
         }
        }
       }
       if (!somefree) {
        break;
       }
       stepsize = MAXSTEP;
       // ## BEGIN loglikelihood FUNCTION ## //
       sumlik = 0.0;
       for (i = 1; i <= I; i++) {
        partlik = 0.0;
        for (j = 1; j <= J; j++) {
         l = cell(i, j) + 1 - baselevel(j);
         p_ijl = 0.0;
         for (k = 1; k <= K; k++) {
          p_ijl += (G(i, k) * P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)));
         }
         if (p_ijl < ZERO) {
          p_ijl = ZERO;
         }
         if (p_ijl < REALBIG) {
          partlik += log(p_ijl);
         }
        }
        sumlik += (cell(i, (J + 1)) * partlik);
       }
       // ## END loglikelihood FUNCTION ## //
       f0 = sumlik;
       bestlik = f0;
       for (halfstep = 0; halfstep < HALFSTEPS && converged == 0; halfstep++) {
        for (k = 1; k <= K; k++) {
         for (j = 1; j <= J; j++) {
          if (FITP(k, j)) {
           for (l = 1; l <= ljlevels(j); l++) {
            P(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k) = old_P(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k) + stepsize * dL_dP(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k);
           }
           // ## BEGIN rescale_P FUNCTION ## //
           sumP = 0.0;
           for (l = 1; l <= ljlevels(j); l++) {
            if (P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) < ZERO) {
             P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) = 0.0;
            }
            sumP += P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k));
           }
           if (sumP < ZERO) {
            sumP = (double)ljlevels(j);
            for (l = 1; l <= ljlevels(j); l++) {
             P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) = 1.0;
            }
           }
           for (l = 1; l <= ljlevels(j); l++) {
            P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) = (P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)) / sumP);
           }
           // ## END rescale_P FUNCTION ## //
          }
         }
        }
        // ## BEGIN loglikelihood FUNCTION ## //
        sumlik = 0.0;
        for (i = 1; i <= I; i++) {
         partlik = 0.0;
         for (j = 1; j <= J; j++) {
          l = cell(i, j) + 1 - baselevel(j);
          p_ijl = 0.0;
          for (k = 1; k <= K; k++) {
           p_ijl += (G(i, k) * P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)));
          }
          if (p_ijl < ZERO) {
           p_ijl = ZERO;
          }
          if (p_ijl < REALBIG) {
           partlik += log(p_ijl);
          }
         }
         sumlik += (cell(i, (J + 1)) * partlik);
        }
        // ## END loglikelihood FUNCTION ## //
        f1 = sumlik;
        if (f1 > bestlik) {
         bestlik = f1;
         for (k = 1; k <= K; k++) {
          for (j = 1; j <= J; j++) {
           for (l = 1; l <= ljlevels(j); l++) {
            new_P(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k) = P(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k);
           }
          }
         }
        } else {
         if ((bestlik > (f0 + ZERO)) && (bestlik > (f1 + ZERO))) {
          break;
         }
        }
        stepsize = (stepsize / 2.0);
       }
       f1 = bestlik;
       loglik[1] = bestlik;
       for (k = 1; k <= K; k++) {
        for (j = 1; j <= J; j++) {
         for (l = 1; l <= ljlevels(j); l++) {
          P(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k) = new_P(((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k);
         }
        }
       }
       sprintf(buffer, "%03d.%03d%s%015.5f%s%015.7f", iter, subiter, "\t", loglik[1], "\t", fabs((loglik[1] - ploglik) / ploglik));
       Rcout << "Fit P (Gradient 1992):\t" << buffer << endl;
       if ((fabs(f0) > ZERO) && (fabs((f1 - f0) / f0) < CTOL)) {
        converged = 1;
       }
       if ((fabs(f0) > ZERO) && (fabs((f1 - f0) / f0) < ITOL) && (subiter >= MAXITER_PARAMETERS)) {
        break;
       }
      }
      Rcout << endl;
      // ## END Fit_P_Gradient_1992 FUNCTION ## //
     }
    }
    // ## BEGIN loglikelihood FUNCTION ## //
    sumlik = 0.0;
    for (i = 1; i <= I; i++) {
     partlik = 0.0;
     for (j = 1; j <= J; j++) {
      l = cell(i, j) + 1 - baselevel(j);
      p_ijl = 0.0;
      for (k = 1; k <= K; k++) {
       p_ijl += (G(i, k) * P((((l * ((J + 1) * (K + 1))) + ((K + 1) * j)) + k)));
      }
      if (p_ijl < ZERO) {
       p_ijl = ZERO;
      }
      if (p_ijl < REALBIG) {
       partlik += log(p_ijl);
      }
     }
     sumlik += (cell(i, (J + 1)) * partlik);
    }
    // ## END loglikelihood FUNCTION ## //
    loglik[1] = sumlik;
    sprintf(buffer, "%-10.4f", loglik[1]);
    Rcout << "Latter Log-Likelihood is:\t" << buffer << endl;
    return(wrap(loglik));
  '

  #####################BEGIN verify.parameters function#


  GoM_Model <- inline::cxxfunction (methods::signature(
    baselevel_ = "numeric",
    ljlevels_ = "numeric",
    cell_ = "matrix",
    gammaalgorithm_ = "character",
    FG_ = "matrix",
    gammafit_ = "bool",
    lambdaalgorithm_ = "character",
    FP_ = "numeric",
    lambdafit_ = "bool",
    FITP_ = "matrix"),
    body = GoM,
    includes = "#include <cstdio>",
    plugin = "Rcpp")


  verify.parameters <- function (case.id,
                                 case.weight,
                                 data.object,
                                 internal.var,
                                 initial.gamma, initial.lambda,
                                 gamma.algorithm, lambda.algorithm,
                                 initial.gamma.object, initial.lambda.object) {
    if (!(gamma.algorithm %in% c("woodbury.1974", "gradient.1992"))) {
      stop("The gamma.algorithm information is wrong ...")
    }
    if (!(initial.gamma %in% c("random", "pure1", "equal.values", "gamma.object"))) {
      stop("The initial.gamma information is wrong ...")
    }
    if (initial.gamma == c("gamma.object") & missing(initial.gamma.object)) {
      stop("The initial GoM scores object is missing ...")
    }
    if (!(lambda.algorithm %in% c("woodbury.1974", "gradient.1992"))) {
      stop("The lambda.algorithm information is wrong ...")
    }
    if (!(initial.lambda %in% c("random", "pure1", "equal.values", "lambda.object"))) {
      stop("The initial.lambda information is wrong ...")
    }
    if (initial.lambda == c("lambda.object") & missing(initial.lambda.object)) {
      stop("The initial pure type probabilities object is missing ...")
    }
    if (is.na(case.id) | !(case.id %in% names(data.object))) {
      stop("The case.id is missing ...")
    }
    if (!is.na(case.weight) && !(case.weight %in% names(data.object))) {
      stop("The case.weight is missing ...")
    }
    if (missing(internal.var)) {
      stop("The internal.var is missing ...")
    } else if (length(internal.var[(c(internal.var) %in% names(data.object) == FALSE)]) > 0) {
      stop("The internal.var information is wrong ...")
    }
  }
  #####################END verify.parameters function#

  #####################BEGIN summary.parameters function#
  summary.parameters <- function (case.id,
                                  case.weight,
                                  data.object,
                                  internal.var,
                                  initial.K, final.K,
                                  initial.gamma, initial.lambda,
                                  gamma.algorithm, lambda.algorithm,
                                  gamma.fit, lambda.fit, order.K, omega.fit) {
    cat(paste("\n*GoMRcpp Summary***********************\n", sep = "", collapse = NULL))
    cat(paste("Input data object---------------------: ", "From input data object", "\n", sep = "", collapse = NULL))
    cat(paste("Number of profiles in initial model---: ", initial.K, "\n", sep = "", collapse = NULL))
    cat(paste("Number of profiles in final model-----: ", final.K, "\n", sep = "", collapse = NULL))
    if(gamma.fit == TRUE) {
      cat(paste("GoM scores algorithm------------------: ", gamma.algorithm, "\n", sep = "", collapse = NULL))
    } else {
      cat(paste("GoM scores algorithm------------------: None\n", sep = "", collapse = NULL))
    }
    if (initial.gamma == "gamma.object") {
      cat(paste("Initial GoM scores--------------------: From input gamma object\n", sep = "", collapse = NULL))
    } else {
      cat(paste("Initial GoM scores--------------------: ", initial.gamma, "\n", sep = "", collapse = NULL))
    }
    cat(paste("Fit GoM scores------------------------: ", gamma.fit, "\n", sep = "", collapse = NULL))
    if (lambda.fit == TRUE) {
      cat(paste("Pure type probabilities algorithm-----: ", lambda.algorithm, "\n", sep = "", collapse = NULL))
    } else {
      cat(paste("Pure type probabilities algorithm-----: None\n", sep = "", collapse = NULL))
    }
    if (initial.lambda == "lambda.object") {
      cat(paste("Initial pure type probabilities-------: From input lambda object\n", sep = "", collapse = NULL))
    } else {
      cat(paste("Initial pure type probabilities-------: ", initial.lambda, "\n", sep = "", collapse = NULL))
    }
    cat(paste("Fit pure type probabilities-----------: ", lambda.fit, "\n", sep = "", collapse = NULL))
    cat(paste("Sort profiles-------------------------: ", order.K, "\n", sep = "", collapse = NULL))
    cat(paste("Records in data object----------------: ", nrow(data.object), "\n", sep = "", collapse = NULL))
    if (omega.fit == TRUE) {
      cat(paste("Unique data patterns------------------: All patterns", "\n", sep = "", collapse = NULL))
    } else {
      cat(paste("Unique data patterns------------------: From input data object", "\n", sep = "", collapse = NULL))
    }
    cat(paste("Case label----------------------------: ", case.id, "\n", sep = "", collapse = NULL))
    if (!is.na(case.weight)) {
      cat(paste("Case weight---------------------------: ", case.weight, "\n", sep = "", collapse = NULL))
    } else {
      cat(paste("Case weight---------------------------: None\n", sep = "", collapse = NULL))
    }
    cat("Internal variables--------------------:", c(internal.var), "\n")
    if (length(names(data.object)[(names(data.object) %in% c(c(case.id), c(case.weight), c(internal.var)) == FALSE)] > 0)) {
      external.var <- names(data.object)[(names(data.object) %in% c(c(case.id), c(case.weight), c(internal.var)) == FALSE)]
    } else {
      external.var <- c("------")
    }
    cat("External variables--------------------:", external.var, "\n")
  }
  #####################END summary.parameters function#

  #####################BEGIN cell.data function#
  cell.data <- function (data.object, case.weight, internal.var, omega.fit) {
    cell <- data.frame(lapply(data.object[, internal.var], as.factor))
    cell <- data.frame(lapply(cell[, internal.var], as.numeric))
    cell <- data.frame(lapply(cell[, internal.var], as.factor))
    if (!is.na(case.weight)) {
      cell[c(case.weight)] <- data.frame(sapply(data.object[, case.weight], as.numeric))
    }
    if (max(mapply(nlevels, cell[, c(internal.var)])) == 2) {
      cat(paste("\n*Note: All internal variables are dichotomous.\n", sep = "", collapse = NULL))
    }
    cell <- cell[do.call(order, cell[c(internal.var)]), ]
    cell$Patterns <- do.call(paste, c(as.list(cell[c(internal.var)]), sep=""))
    if (omega.fit == TRUE) {
      cellomega <- as.data.frame(stats::ftable(cell[c(internal.var)]))
      cellomega$Freq <- as.numeric(cellomega$Freq + 1)
      cellomega <- cellomega[do.call(order, cellomega[c(internal.var)]), ]
      cellomega$Patterns <- do.call(paste, c(as.list(cellomega[c(internal.var)]), sep=""))
      if (is.na(case.weight)) {
        cell <- cellomega[, c(c("Patterns"), c(internal.var), c("Freq"))]
      } else if (!is.na(case.weight)) {
        cellStringFreq <- stats::aggregate(cell[c(case.weight)], list(cell$Patterns), FUN = sum, simplify = TRUE)
        names(cellStringFreq) <- c("Patterns", "Freqomega")
        cellomega <- merge(cellomega, cellStringFreq, by = c("Patterns"), all.x = TRUE)
        cellomega$Freqomega <- as.numeric(cellomega$Freqomega + 1)
        cellomega$Freqomega[is.na(cellomega$Freqomega)] <- 1
        cell <- cellomega[, c(c("Patterns"), c(internal.var), c("Freqomega"))]
        names(cell) <- c("Patterns", c(internal.var), "Freq")
      }
      cat(paste("\n*Note: ", nrow(cell), " unique data patterns has found for combinations of the all patterns.\n\n", sep = "", collapse = NULL))
    }
    if (omega.fit == FALSE) {
      if (is.na(case.weight)) {
        cell$FreqCell <- sequence(rle(cell$Patterns)$lengths)
        cellStringFreq <- stats::aggregate(cell$FreqCell, list(cell$Patterns), FUN = max, simplify = TRUE)
        names(cellStringFreq) <- c("Patterns", "Freq")
      } else if (!is.na(case.weight)) {
        cellStringFreq <- stats::aggregate(cell[c(case.weight)], list(cell$Patterns), FUN = sum, simplify = TRUE)
        names(cellStringFreq) <- c("Patterns", "Freq")
      }
      cell <- merge(cell, cellStringFreq, by = c("Patterns"))
      cell <- cell[, c(c("Patterns"), c(internal.var), c("Freq"))]
      cell <- unique(cell)
      cat(paste("\n*Note: ", nrow(cell), " unique data patterns has found in data object.\n\n", sep = "", collapse = NULL))
    }
    cell <- rbind(rep(NA, (length(c(internal.var)) + 2)), cell)
    row.names(cell) <- NULL
    return (cell)
  }
  #####################END cell.data function#

  #####################BEGIN ljlevels.function function#
  ljlevels.function <- function (cell, internal.var) {
    ljlevels <- sapply(c(c(NA), c(cell[c(internal.var)])), nlevels)
    return (ljlevels)
  }
  #####################END ljlevels.function function#

  #####################BEGIN l.levels.j.function function#
  l.levels.j.function <- function (cell, internal.var) {
    l.levels.j <- (sapply(sapply(c(c(NA), c(cell[c(internal.var)])), levels), as.numeric))
    for (i in 1:length(internal.var)) {
      if (min(l.levels.j[[i+1]]) < 1) {
        stop(paste0("\n\nThe data.object can only encompass codes values from 1 (one):\n\n",
                    "\tTherefore, you must to avoid to use the 0 (zero) code ...\n\n\n"))
      }
    }
    return (l.levels.j)
  }
  #####################END l.levels.j.function function#

  #####################BEGIN parameter.FG function#
  parameter.FG <- function (initial.K, initial.gamma, cell, initial.gamma.object) {
    if (initial.gamma == c("random")) {
      FG <- as.data.frame(matrix(NA, nrow(cell), initial.K, byrow = T))
      for (i in 1:nrow(cell)) {
        ki <- sample(c(1:initial.K), initial.K, replace = FALSE, prob = NULL)
        random.gamma <- c(rep(as.double(0), initial.K))
        for (k in 1:initial.K) {
          if (k == 1) {
            random.gamma[k] <- as.double(stats::runif(1, min = as.double(0), max = as.double(1)))
          } else {
            random.gamma[k] <- as.double(stats::runif(1, min = as.double(0), max = (as.double(1) - sum(random.gamma))))
          }
        }
        for (k in 1:initial.K) {
          if (k == 1) {
            FG[i, ki[k]] <- random.gamma[k]
          } else if (k != 1) {
            if (sum(FG[i, 1:initial.K], random.gamma[k], na.rm = TRUE) <= as.double(1)) {
              if (k != length(ki)) {
                FG[i, ki[k]] <- random.gamma[k]
              } else if (k == length(ki)) {
                FG[i, ki[k]] <- (as.double(1) - (sum(FG[i, 1:initial.K], na.rm = TRUE)))
              }
            } else if (sum(FG[i, 1:initial.K], random.gamma[k], na.rm = TRUE) > as.double(1)) {
              FG[i, ki[k]] <- (as.double(1) - (sum(FG[i, 1:initial.K], na.rm = TRUE)))
            }
          }
        }
      }
    } else if (initial.gamma == c("pure1")) {
      FG <- as.data.frame(matrix(c(as.double(1), rep(as.double(0), (initial.K - 1))), nrow(cell), initial.K, byrow = T))
    } else if (initial.gamma == c("equal.values")) {
      FG <- as.data.frame(matrix(c(as.double(1) / initial.K), nrow(cell), initial.K, byrow = T))
    } else if (initial.gamma == c("gamma.object")) {
      FG <- initial.gamma.object
    }
    if (initial.gamma != c("gamma.object")) {
      FG <- cbind(rep(NA, nrow(cell)), FG)
      FG[1, ] <- NA
    } else {
      FG <- cbind(rep(NA, (nrow(cell) - 1)), FG)
      FG <- rbind(rep(NA, (initial.K + 1)), FG)
    }
    names(FG) <- c("Patterns", paste("k", 1:initial.K, sep = ""))
    FG[, "Patterns"] <- NA
    return (FG)
  }
  #####################END parameter.FG function#

  #####################BEGIN fit.P function#
  fit.P <- function (initial.K, initial.lambda, ljlevels, initial.lambda.object) {
    FITP <- matrix(as.integer(1), (initial.K + 1), length(ljlevels), byrow = T)
    FITP[1, ] <- NA
    FITP[, 1] <- NA
    if (initial.lambda == c("lambda.object")) {
      for (k in 2:(initial.K + 1)) {
        for (j in 2:length(ljlevels)) {
          for (l in 2:(ljlevels[[j]] + 1)) {
            if (initial.lambda.object[(k - 1), (j - 1), (l - 1)] < as.double(0)) {
              FITP[k, j] <- as.integer(0)
            }
          }
        }
      }
    }
    return (FITP)
  }
  #####################BEGIN fit.P function#

  #####################BEGIN parameter.FP function#
  parameter.FP <- function (initial.K, initial.lambda, ljlevels, initial.lambda.object) {
    FP <- array(c((as.double(1) / initial.K)), c((initial.K + 1), length(ljlevels), (max(ljlevels) + 1)), dimnames = list(c(paste("k", 0:initial.K, sep = "")), c(paste("j", 0:(length(ljlevels) - 1), sep = "")), c(paste("l", 0:max(ljlevels), sep = ""))))
    FP[1, , ] <- NA
    FP[, 1, ] <- NA
    FP[, , 1] <- NA
    if (initial.lambda == c("pure1")) {
      for (j in 2:length(ljlevels)) {
        FP[2, j, 2] <- as.double(1)
        for (l in 3:(ljlevels[[j]] + 1)) {
          FP[2, j, l] <- as.double(0)
        }
      }
    }
    if (initial.lambda == c("random")) {
      for (k in 2:(initial.K + 1)) {
        for (j in 2:length(ljlevels)) {
          li <- sample(c(2:(ljlevels[[j]] + 1)), ljlevels[[j]], replace = FALSE, prob = NULL)
          random.lambda <- c(rep(as.double(0), (ljlevels[[j]] + 1)))
          for (l in 2:(ljlevels[[j]] + 1)) {
            if (l == 2) {
              random.lambda[l] <- as.double(stats::runif(1, min = as.double(0), max = as.double(1)))
            } else {
              random.lambda[l] <- as.double(stats::runif(1, min = as.double(0), max = (as.double(1) - sum(random.lambda))))
            }
          }
          for (l in 2:(ljlevels[[j]] + 1)) {
            if (l == 2) {
              FP[k, j, 2:(ljlevels[[j]] + 1)] <- NA
              FP[k, j, li[l - 1]] <- random.lambda[l]
            } else if (l != 2) {
              if (sum(FP[k, j, 2:(ljlevels[[j]] + 1)], random.lambda[l], na.rm = TRUE) <= as.double(1)) {
                if ((l - 1) != length(li)) {
                  FP[k, j, li[l - 1]] <- random.lambda[l]
                } else if ((l - 1) == length(li)) {
                  FP[k, j, li[l - 1]] <- as.double(1) - sum(FP[k, j, 2:(ljlevels[[j]] + 1)], na.rm = TRUE)
                }
              } else if (sum(FP[k, j, 2:(ljlevels[[j]] + 1)], random.lambda[l], na.rm = TRUE) > as.double(1)) {
                FP[k, j, li[l - 1]] <- as.double(1) - sum(FP[k, j, 2:(ljlevels[[j]] + 1)], na.rm = TRUE)
              }
            }
          }
        }
      }
    }
    if (initial.lambda == c("lambda.object")) {
      for (k in 2:(initial.K + 1)) {
        for (j in 2:length(ljlevels)) {
          for (l in 2:(ljlevels[[j]] + 1)) {
            if (initial.lambda.object[(k - 1), (j - 1), (l - 1)] < as.double(0)) {
              FP[k, j, l] <- -(initial.lambda.object[(k - 1), (j - 1), (l - 1)])
            } else {
              FP[k, j, l] <- initial.lambda.object[(k - 1), (j - 1), (l - 1)]
            }
          }
        }
      }
    }
    for (k in 2:(initial.K + 1)) {
      for (j in 2:length(ljlevels)) {
        for (l in 2:(max(ljlevels) + 1)) {
          if (l > (max(ljlevels[[j]]) + 1)) {
            FP[k, j, l] <- NA
          }
        }
      }
    }
    return (FP)
  }
  #####################END parameter.FP function#

  #####################BEGIN v.order.K function#
  v.order.K <- function (initial.K, cell, ljlevels, FP) {
    N <- sum(cell$Freq, na.rm = TRUE)
    Zc <- 2.58
    for (l in 2:((max(ljlevels) + 1) - 1)) {
      v.order <- c(rep(as.double(0), initial.K))
      for (k in 2:(initial.K + 1)) {
        Pc1 <- 0
        Pc2 <- 0
        for (j in 2:length(ljlevels)) {
          if (l < (max(ljlevels[[j]]) + 1)) {
            if ((as.double(sum(FP[k, j, 2:l])) != as.double(0)) == TRUE) {
              if (v.order[k - 1] == as.double(0)) {
                v.order[k - 1] <- as.double(sum(FP[k, j, 2:l]) ^ (1 / N))
              } else {
                v.order[k - 1] <- prod(v.order[k - 1], as.double(sum(FP[k, j, 2:l]) ^ (1 / N)))
              }
            } else {
              Pc1 <- Pc1 + 1
              Pc2 <- Pc2 + ((length(ljlevels) - 1) + 2) - j
            }
          } else {
            l_l <- ((max(ljlevels[[j]]) + 1) - 1)
            if ((as.double(sum(FP[k, j, 2:l_l])) != as.double(0)) == TRUE) {
              if (v.order[k - 1] == as.double(0)) {
                v.order[k - 1] <- as.double(sum(FP[k, j, 2:l_l]) ^ (1 / N))
              } else {
                v.order[k - 1] <- prod(v.order[k - 1], as.double(sum(FP[k, j, 2:l_l]) ^ (1 / N)))
              }
            } else {
              Pc1 <- Pc1 + 1
              Pc2 <- Pc2 + ((length(ljlevels) - 1) + 2) - j
            }
          }
        }
        v.order[k - 1] <- (v.order[k - 1] / (1 + ((Pc2 / sum((length(ljlevels) - 1):1)) * Pc1)))
      }
      p.v.order <- matrix(NA, initial.K, initial.K)
      for (k in 1:initial.K) {
        for (k_k in 1:initial.K) {
          if (k != k_k) {
            p.v.order[k, k_k] <- ((v.order[k] - v.order[k_k]) / sqrt((((N * v.order[k]) + (N * v.order[k_k])) / (N + N)) * (1 - (((N * v.order[k]) + (N * v.order[k_k])) / (N + N))) * ((N + N) / (N * N))))
          }
        }
      }
      if (min(abs(p.v.order), na.rm = TRUE) > Zc) {
        break
      }
    }
    if (min(abs(p.v.order), na.rm = TRUE) > Zc) {
      v.order <- order(v.order, decreasing = TRUE)
    } else {
      v.order <- NULL
    }
    return (v.order)
  }
  #####################END v.order.K function#

  #####################BEGIN v.order.P function#
  v.order.P <- function (initial.K, ljlevels, v.order, beforeP) {
    afterP <- array(NA, c((initial.K + 1), length(ljlevels), (max(ljlevels) + 1)), dimnames = list(c(paste("k", 0:initial.K, sep = "")), c(paste("j", 0:(length(ljlevels) - 1), sep = "")), c(paste("l", 0:max(ljlevels), sep = ""))))
    k <- 1
    for (k_k in v.order) {
      k <- k + 1
      for (j in 2:length(ljlevels)) {
        for (l in 2:(max(ljlevels[[j]]) + 1)) {
          afterP[k, j, l] <- beforeP[(k_k + 1), j, l]
        }
      }
    }
    return(afterP)
  }
  #####################END v.order.P function#

  #####################BEGIN data.gamma function#
  data.gamma <- function (case.id, data.object, internal.var, case.weight,
                          initial.K, cell, ljlevels, l.levels.j,
                          order.K, v.order, IG, FG,
                          omega.fit, dec.char) {
    IG[, "Patterns"] <- cell$Patterns
    IG <- IG[-1, ]
    row.names(IG) <- NULL
    FG[, "Patterns"] <- cell$Patterns
    FG$FreqPatterns <- cell$Freq
    FG <- FG[-1,]
    row.names(FG) <- NULL
    if ((order.K == TRUE) && (!is.null(v.order))){
      v.order <- (v.order + 1)
      IG <- IG[, c(1, v.order)]
      FG <- FG[, c(1, v.order, (max(v.order) + 1))]
    }
    names(IG) <- c("Patterns", paste("initial_gik", 1:initial.K, sep = ""))
    names(FG) <- c("Patterns", paste("final_gik", 1:initial.K, sep = ""), "FreqPatterns")
    names.object <- names(data.object)
    data.object <- data.object[do.call(order, data.object[, internal.var]), ]
    data.object$Patterns <- do.call(paste, c(as.list(data.object[, internal.var]), sep=""))
    if (omega.fit == TRUE) {
      data.object <- merge(data.object, IG, by = c("Patterns"), all.y = TRUE)
      data.object <- merge(data.object, FG, by = c("Patterns"), all.y = TRUE)
    } else {
      data.object <- merge(data.object, IG, by = c("Patterns"))
      data.object <- merge(data.object, FG, by = c("Patterns"))
    }
    data.object <- data.object[c(c(names.object), c("Patterns"), c("FreqPatterns"), c(paste("initial_gik", 1:initial.K, sep = "")), c(paste("final_gik", 1:initial.K, sep = "")))]
    data.object <- data.object[do.call(order, list(data.object[, case.id])), ]
    if (omega.fit == TRUE) {
      if (!is.na(case.weight)) {
        data.object[case.weight][is.na(data.object[case.weight])] <- 1
      }
      st <- 0
      for (j in 1:(length(ljlevels) - 1)) {
        st <- (st + as.numeric((nchar(max(l.levels.j[[j + 1]])))))
        data.object[, c(internal.var[j])] <- as.data.frame(substr(data.object[, "Patterns"], st, ((st + as.numeric(nchar(max(l.levels.j[[j + 1]])))) - 1)))
      }
    }
    row.names(data.object) <- NULL
    nf <- 0
    repeat {
      nf <- (nf + 1)
      dataoutput <- paste (getwd(), "/GoMK", initial.K, "(", nf, ")", ".TXT", sep = "", collapse = NULL)
      logname <- paste (getwd(), "/LogGoMK", initial.K, "(", nf, ")", ".TXT", sep = "", collapse = NULL)
      if ((file.exists(dataoutput) == FALSE) && (file.exists(logname) == FALSE)) {
        break
      }
    }
    utils::write.table(data.object, file = dataoutput, quote = FALSE, sep = " ", eol = "\n", na = ".", dec = dec.char, row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")
    return (nf)
  }
  #####################END data.gamma function#

  #####################BEGIN loggom function#
  loggom <- function (case.id,
                      case.weight,
                      data.object,
                      internal.var,
                      initial.K, final.K,
                      initial.gamma, initial.lambda,
                      gamma.algorithm, lambda.algorithm,
                      gamma.fit, lambda.fit, order.K, omega.fit,
                      cell, ljlevels, l.levels.j, IP, FP, loglik, nf, dec.char, v.order) {
    output <- paste (getwd(), "/LogGoMK", initial.K, "(", nf, ")", ".TXT", sep = "", collapse = NULL)
    file.create(output)
    sink(output)
    cat(paste(date(), "\n\n", sep = "", collapse = NULL))
    summary.parameters(case.id,
                       case.weight,
                       data.object,
                       internal.var,
                       initial.K, final.K,
                       initial.gamma, initial.lambda,
                       gamma.algorithm, lambda.algorithm,
                       gamma.fit, lambda.fit, order.K, omega.fit)
    if (max(ljlevels, na.rm = TRUE) == 2) {
      cat(paste("\n\n*Note: All internal variables are dichotomous.\n", sep = "", collapse = NULL))
    }
    if (omega.fit == TRUE) {
      cat(paste("\n\n*Note ", (nrow(cell) - 1), " unique data patterns (I) has found for combinations of the all patterns.\n", sep = "", collapse = NULL))
    } else if (omega.fit == FALSE) {
      cat(paste("\n\n*Note ", (nrow(cell) - 1), " unique data patterns (I) has found in data object.\n", sep = "", collapse = NULL))
    }
    charnamevar <- max(sapply(internal.var, nchar))
    cat(paste("\n\nFrequency Table Original Data:\n", sep = "", collapse = NULL))
    for (j in 1:(length(ljlevels) - 1)) {
      if (!is.na(case.weight)) {
        n <- stats::xtabs(data.object[, case.weight] ~ data.object[, internal.var[j]], data.object)
      } else {
        aux <- internal.var[j]
        n <- table(data.object[, aux])
      }
      p <- prop.table(n)
      for (l in l.levels.j[[(j + 1)]]) {
        if (l == (min(l.levels.j[[(j + 1)]]))) {
          if (j == 1) {
            t <- do.call(paste, as.list(rep("", charnamevar)))
            cat(paste(t, " \t", "   ", "\tn\t%\n", sep = "", collapse = NULL))
          }
          t <- do.call(paste, as.list(rep("", (charnamevar - (nchar(internal.var[j]))))))
          cat(paste(internal.var[j], t, sep = "", collapse = NULL))
          cat(paste("\t", "l", l, "\t", sep = "", collapse = NULL))
        } else {
          t <- do.call(paste, as.list(rep("", charnamevar)))
          cat(paste(t, " \t", "l", l, "\t", sep = "", collapse = NULL))
        }
        if ((p[[l]] * 100) < 10) {
          cat(paste(format(n[[l]], nsmall = 0, decimal.mark = dec.char), "\t", "0", format(round((p[[l]] * 100), 3), nsmall = 3, decimal.mark = dec.char), "\n", sep = "", collapse = NULL))
        } else {
          cat(paste(format(n[[l]], nsmall = 0, decimal.mark = dec.char), "\t", format(round((p[[l]] * 100), 3), nsmall = 3, decimal.mark = dec.char), "\n", sep = "", collapse = NULL))
        }
      }
      cat(paste("\n", sep = "", collapse = NULL))
    }
    LJ <- 0
    for (i in 1:2) {
      if (i == 1) {
        cat(paste("\nPrimal Pure Type Probabilities:\n", sep = "", collapse = NULL))
      } else {
        cat(paste("\nLatter Pure Type Probabilities:\n", sep = "", collapse = NULL))
      }
      for (j in 2:length(ljlevels)) {
        if (i == 1) {
          LJ <- sum(LJ, max(l.levels.j[[j]]))
        }
        for (l in 2:(ljlevels[[j]] + 1)) {
          if ((j == 2) && (l == (min(l.levels.j[[j]]) + 1))) {
            t <- do.call(paste, as.list(rep("", charnamevar)))
            cat(paste(t, "\t", "   ", sep = "", collapse = NULL))
            for (k in 2:(initial.K + 1)) {
              cat(paste("\tk", (k - 1), "    ", sep = "", collapse = NULL))
            }
            cat(paste("\n", sep = "", collapse = NULL))
            t <- do.call(paste, as.list(rep("", (charnamevar - (nchar(internal.var[j - 1]))))))
            cat(paste(internal.var[j - 1], t, sep = "", collapse = NULL))
            cat(paste("\t", "l", l - 1, sep = "", collapse = NULL))
          } else if (l == (min(l.levels.j[[j]]) + 1)) {
            t <- do.call(paste, as.list(rep("", (charnamevar - (nchar(internal.var[j - 1]))))))
            cat(paste(internal.var[j - 1], t, sep = "", collapse = NULL))
            cat(paste("\t", "l", l - 1, sep = "", collapse = NULL))
          }
          for (k in 2:(initial.K + 1)) {
            if((k == 2) && (l != (min(l.levels.j[[j]]) + 1))) {
              t <- do.call(paste, as.list(rep("", charnamevar)))
              cat(paste(t, " \t", "l", l - 1, sep = "", collapse = NULL))
            }
            if (i == 1) {
              cat(paste("\t", format(round(IP[k, j, l], 4), nsmall = 4, decimal.mark = dec.char), sep = "", collapse = NULL))
            } else {
              cat(paste("\t", format(round(FP[k, j, l], 4), nsmall = 4, decimal.mark = dec.char), sep = "", collapse = NULL))
            }
          }
          cat(paste("\n", sep = "", collapse = NULL))
        }
        cat(paste("\n", sep = "", collapse = NULL))
      }
    }
    if (is.null(v.order)) {
      cat(paste("\n*Note: Could not organize pure type probabilities with a confidence interval of 99.0%.\n\n", sep = "", collapse = NULL))
    }
    AIC <- ((2 * ((initial.K * LJ) + (initial.K * (sum(cell$Freq, na.rm = TRUE))))) - (2 * loglik[2]))
    cat(paste("\nPrimal Log-Likelihood is:       \t", format(round(loglik[1], 4), nsmall = 4, decimal.mark = dec.char) , "\n", sep = "", collapse = NULL))
    cat(paste("\n\nLatter Log-Likelihood is:       \t", format(round(loglik[2], 4), nsmall = 4, decimal.mark = dec.char) , "\n", sep = "", collapse = NULL))
    cat(paste("\n\nAkaike Information Criterion:   \t", format(round(AIC, 4), nsmall = 4, decimal.mark = dec.char) , "\n", sep = "", collapse = NULL))
    cat(paste("\n\nLambda-Marginal Frequency Ratio (LMFR):\n", sep = "", collapse = NULL))
    table <- utils::capture.output({
      for (j in 2:length(ljlevels)) {
        if (omega.fit == FALSE) {
          if (!is.na(case.weight)) {
            n <- stats::xtabs(data.object[, case.weight] ~ data.object[, internal.var[(j - 1)]], data.object)
          } else {
            aux <- internal.var[(j-1)]
            n <- table(data.object[, aux])
          }
        } else {
          n <- stats::xtabs(cell[-1, "Freq"] ~ cell[-1, internal.var[(j - 1)]], cell)
        }
        p <- prop.table(n)
        for (l in l.levels.j[[j]]) {
          if (l == (min(l.levels.j[[j]]))) {
            if (j == 2) {
              t <- do.call(paste, as.list(rep("", charnamevar)))
              cat(paste(t, "\t", "   ", "\tn\t%", sep = "", collapse = NULL))
              for (k in 2:(initial.K + 1)) {
                cat(paste("\tk", (k - 1), "    ", sep = "", collapse = NULL))
              }
              for (k in 2:(initial.K + 1)) {
                cat(paste("\tk", (k - 1), "/%lj", sep = "", collapse = NULL))
              }
              cat(paste("\n", sep = "", collapse = NULL))
            }
            t <- do.call(paste, as.list(rep("", (charnamevar - (nchar(internal.var[(j - 1)]))))))
            cat(paste(internal.var[(j - 1)], t, sep = "", collapse = NULL))
            cat(paste("\t", "l", l, "\t", sep = "", collapse = NULL))
          } else {
            t <- do.call(paste, as.list(rep("", charnamevar)))
            cat(paste(t, " \t", "l", l, "\t", sep = "", collapse = NULL))
          }
          if ((p[[l]] * 100) < 10) {
            cat(paste(format(n[[l]], nsmall = 0, decimal.mark = dec.char), "\t", "0", format(round((p[[l]] * 100), 3), nsmall = 3, decimal.mark = dec.char), sep = "", collapse = NULL))
          } else {
            cat(paste(format(n[[l]], nsmall = 0, decimal.mark = dec.char), "\t", format(round((p[[l]] * 100), 3), nsmall = 3, decimal.mark = dec.char), sep = "", collapse = NULL))
          }
          for (k in 2:(initial.K + 1)) {
            cat(paste("\t", format(round(FP[k, j, (l + 1)], 4), nsmall = 4, decimal.mark = dec.char), sep = "", collapse = NULL))
          }
          for (k in 2:(initial.K + 1)) {
            cat(paste("\t", format(round((FP[k, j, (l + 1)] / p[[l]]), 4), nsmall = 4, decimal.mark = dec.char), sep = "", collapse = NULL))
          }
          cat(paste("\n", sep = "", collapse = NULL))
        }
        if (j != length(ljlevels)) {
          cat(paste("\n", sep = "", collapse = NULL))
        }
      }
    })
    writeLines(table)
    sink()
    if (is.null(v.order)) {
      cat(paste("\n*Note: Could not organize pure type probabilities with a confidence interval of 99.0%.\n", sep = "", collapse = NULL))
    }
    cat(paste("\n*Note: Saved frequency table, pure type probabilities and loglikelihood values to ", output, ".\n", sep = "", collapse = NULL))
    cat(paste("\n*Note: Saved original data with GoM scores to " , getwd(), "/GoMK", initial.K, "(", nf, ")", ".TXT", ".\n", sep = "", collapse = NULL))
    return(table)
  }
  #####################END loggom function#

  if (!(is.data.frame (data.object))) {
    stop("The data.object is not a data frame ...")
  }
  if (initial.K < 2) {
    stop("The initial.K information is wrong ...")
  }
  if (final.K < 2) {
    stop("The final.K information is wrong ...")
  }
  gamma.algorithm <- gamma.algorithm[1]
  initial.gamma <- initial.gamma[1]
  if ((gamma.fit != TRUE) & (gamma.fit != FALSE)) {
    stop("The gamma.fit information is wrong ...")
  }
  lambda.algorithm <- lambda.algorithm[1]
  initial.lambda <- initial.lambda[1]
  if ((lambda.fit != TRUE) & (lambda.fit != FALSE)) {
    stop("The lambda.fit information is wrong ...")
  }
  if (is.null(internal.var)) {
    internal.var <- names(data.object)[!(names(data.object) %in% c(c(case.id), c(case.weight)))]
  }
  if ((order.K != TRUE) & (order.K != FALSE)) {
    stop("The order.K information is wrong ...")
  }
  if ((omega.fit != TRUE) & (omega.fit != FALSE)) {
    stop("The omega.fit information is wrong ...")
  }
  dec.char <- substr(dec.char, 1, 1)
  if ((dec.char != ",") && (dec.char != ".")) {
    dec.char <- "."
  }
  pathfolder <- getwd()
  default.scipen <- options("scipen")
  default.digits <- options("digits")
  options(digits = 15)
  FINAL.PARAMETERS <- vector("list")
  for (initial.K in initial.K:final.K) {
    verify.parameters(case.id,
                      case.weight,
                      data.object,
                      internal.var,
                      initial.gamma, initial.lambda,
                      gamma.algorithm, lambda.algorithm,
                      initial.gamma.object, initial.lambda.object)
    summary.parameters(case.id,
                       case.weight,
                       data.object,
                       internal.var,
                       initial.K, final.K,
                       initial.gamma, initial.lambda,
                       gamma.algorithm, lambda.algorithm,
                       gamma.fit, lambda.fit, order.K, omega.fit)
    cell <- cell.data(data.object, case.weight, internal.var, omega.fit)
    ljlevels <- ljlevels.function(cell, internal.var)
    l.levels.j <- l.levels.j.function(cell, internal.var)
    baselevel <- c(rep(0, length(ljlevels)))
    for (j in 2:length(ljlevels)) {
      baselevel[j] <- l.levels.j[[j]][1]
    }
    names(baselevel) <- names(ljlevels)
    if (final.K > (nrow(cell) - 1)) {
      stop("The final.K information can not be greater than value of (I) ...")
    }
    cell <- as.matrix(cell)
    row.names(cell) <- NULL
    cell[, 1] <- NA
    cell <- apply(cell, 2, as.numeric)
    if (initial.gamma == c("gamma.object")) {
      if (!is.data.frame(initial.gamma.object)) {
        stop("The initial.gamma.object is not a data frame ...")
      }
    }

    ###############################C++###############################
    loglik <- vector(mode = "list", length = MC_iter)
    FG_list <- vector(mode = "list", length = MC_iter)
    FP_list <- vector(mode = "list", length = MC_iter)
    FITP_list <- vector(mode = "list", length = MC_iter)
    IG_list <- vector(mode = "list", length = MC_iter)
    cat("\n Monte Carlo Iteration: \n")
    pb = utils::txtProgressBar(min = 1, max = MC_iter, initial = 1, style = 3, width = 60)
    for(i in 1:MC_iter){
      FG <- parameter.FG(initial.K, initial.gamma, cell, initial.gamma.object)
      FG <- as.matrix(FG)
      IG <- as.data.frame(round(FG, 4))
      if (nrow(IG) != nrow(cell)) {
        stop("The number of lines of initial.gamma.object is wrong ...")
      }
      if (initial.lambda == c("lambda.object")) {
        if (!is.array(initial.lambda.object)) {
          stop("The initial.lambda.object is not a array ...")
        }
      }
      FITP <- fit.P(initial.K, initial.lambda, ljlevels, initial.lambda.object)
      FP <- parameter.FP(initial.K, initial.lambda, ljlevels, initial.lambda.object)
      FP <- as.array(FP)

      FG_list[[i]] <- FG
      FP_list[[i]] <- FP
      FITP_list[[i]] <- FITP
      IG_list[[i]] <- IG

      invisible(utils::capture.output(loglik[[i]] <- GoM_Model(baselevel, ljlevels, cell,
                                                        gamma.algorithm, FG, gamma.fit,
                                                        lambda.algorithm, FP, lambda.fit, FITP)))

      utils::setTxtProgressBar(pb,i)



    }
    cat("\n \n")
    maxloglik <- vector(length = MC_iter)
    for(i in 1:MC_iter){
      maxloglik[i] <- loglik[[i]][2]
    }
    FG <- FG_list[[which.max(maxloglik)]]
    FP <- FP_list[[which.max(maxloglik)]]
    FITP <- FITP_list[[which.max(maxloglik)]]
    IG <- IG_list[[which.max(maxloglik)]]
    loglik <- loglik[[which.max(maxloglik)]]
    IP <- array(NA, c((initial.K + 1), length(ljlevels), (max(ljlevels) + 1)), dimnames = list(c(paste("k", 0:initial.K, sep = "")), c(paste("j", 0:(length(ljlevels) - 1), sep = "")), c(paste("l", 0:max(ljlevels), sep = ""))))
    for (k in 2:(initial.K + 1)) {
      for (j in 2:length(ljlevels)) {
        for (l in 2:(ljlevels[[j]] + 1)) {
          IP[k, j, l] <- round(FP[k, j, l], 4)
        }
      }
    }
    #################################################################
    newfolder <- paste("K", initial.K, sep = "", collapse = NULL)
    if (file.exists(newfolder) == FALSE) {
      dir.create(newfolder, showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    setwd(paste(pathfolder, "/K", initial.K, sep = "", collapse = NULL))
    cell <- as.data.frame(cell)
    cell$Patterns <- do.call(paste, c(as.list(cell[c(internal.var)]), sep=""))
    FG <- as.data.frame(round(FG, 4))
    FP <- round(FP, 4)
    if (order.K == TRUE) {
      v.order <- v.order.K(initial.K, cell, ljlevels, FP)
      if (!is.null(v.order)) {
        IP <- v.order.P(initial.K, ljlevels, v.order, IP)
        FP <- v.order.P(initial.K, ljlevels, v.order, FP)
      }
    } else {
      v.order <- as.integer(0)
    }
    options(scipen = 9999999)
    options(digits = 4)
    nf <- data.gamma(case.id, data.object, internal.var, case.weight,
                     initial.K, cell, ljlevels, l.levels.j,
                     order.K, v.order, IG, FG,
                     omega.fit, dec.char)
    table <- loggom(case.id,
           case.weight,
           data.object,
           internal.var,
           initial.K, final.K,
           initial.gamma, initial.lambda,
           gamma.algorithm, lambda.algorithm,
           gamma.fit, lambda.fit, order.K, omega.fit,
           cell, ljlevels, l.levels.j, IP, FP, loglik, nf, dec.char, v.order)
    setwd(pathfolder)
    FG <- FG[-1, -1]
    row.names(FG) <- NULL
    if ((order.K == TRUE) && (!is.null(v.order))){
      FG <- FG[, v.order]
    }
    names(FG) <- c(paste("gik", 1:initial.K, sep = ""))
    FINAL.PARAMETERS[[paste0("K", initial.K)]]$Gik <- FG
    FP <- FP[-1, -1, -1]
    FINAL.PARAMETERS[[paste0("K", initial.K)]]$Pkjl <- FP
    FINAL.PARAMETERS[[paste0("K", initial.K)]]$Likelihood <- loglik[2] #adicionei o AIC a partir daqui
    LJ <- 0
    for(i in 2:length(ljlevels)){
      LJ <- sum(LJ, max(l.levels.j[[i]]))
    }
    AIC <- ((2 * ((initial.K * LJ) + (initial.K * (sum(cell$Freq, na.rm = TRUE))))) - (2 * loglik[2]))
    FINAL.PARAMETERS[[paste0("K", initial.K)]]$AIC <- AIC
    FINAL.PARAMETERS[[paste0("K", initial.K)]]$Table <- table
  }
  options(scipen = default.scipen[[1]])
  options(digits = default.digits[[1]])
  return(FINAL.PARAMETERS)
}
#'
#'
#' @examples
#' \dontrun{
#' data <- data.frame(x1 = round(stats::runif(n = 500, 1, 2), 0),
#'                    x2 = round(stats::runif(n = 500, 1, 3), 0),
#'                    x3 = round(stats::runif(n = 500, 1, 4), 0),
#'                    x4 = round(stats::runif(n = 500, 1, 5), 0),
#'                    Id = 1:500)
#' gom_ml(data.object = data, case.id = "Id", initial.lambda = "random")
#' }


