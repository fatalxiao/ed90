# bootBC.ci.R

# last update 28-DEC-2011

# The R software implementation (functions) is cited in Pace NL, Stylianou, MP: Advances in and Limitations of Up-and-down Methodology. Anesthesiology 2007; 107:144-152.
# R:A Language and Environment for Statistical Computing (R Foundation for Statistical Computing, Vienna, Austria) is available as Free Software under the terms of the Free Software Foundation's GNU General Public License.

# the Canty Bootstrap package does not provide bca (bias corrected accelerated)
# confidence intervals for parametric bootstrap objects
# this function implements the bias corrected (bc) confidence intervals used by Stylianou et al
# this function is tested for boot objects estimated by the package 'boot' version 1.3-3 under R version 2.14.0
# this function is anticipated to be upwardly compatible as the version of boot and R update
# the arguments tObserved and tBoot are values from the boot object
# tObserved = t0 is the vector of observed statistics
# tBoot = t is the matrix with R rows each of which is a bootstrap replicate of the statistics
# if the 95% bc confidence intervals of the ith statistic are desired from the boot object test.boot
# estimated from the databases test.df and testPava.df
# the call would be bootBC.ci(test.boot$t0[i], test.boot$t[, i]
# the output of bootBC.ci has named elements which should be self explanatory
# the bootstrap parametric resample function produces further doses strictly within the dose sample space, thus the CIs are similarly bounded

#' @title Estimate Confidence Interval of ED50 Using Isotonic Regression
#' @description Estimate confidence interval of ED50 using isotonic regression based on bootstrap method.
#' @param tObserved the vector of observed statistics.
#' @param tBoot The matrix with R rows each of which is a bootstrap replicate of the statistics.
#' @param conf Confidence level.
#' @import stats
#' @export
#' @examples
#' library(ed50)
#' library(boot)
#' pavaData <- preparePava(groupS)
#' bootResult <- boot(data = groupS,
#'               statistic = bootIsotonicRegression,
#'                       R = 10,
#'                     sim = 'parametric',
#'                 ran.gen = bootIsotonicResample,
#'                     mle = list(baselinePava = pavaData,
#'                                   firstDose = 2.5,
#'                           PROBABILITY.GAMMA = 0.5),
#'            baselinePava = pavaData,
#'       PROBABILITY.GAMMA = 0.5)
#' bootBC.ci(tObserved = bootResult$t0[3],
#'               tBoot = bootResult$t[, 3],
#'                conf = 0.95)


bootBC.ci <- function(tObserved, tBoot, conf = 0.95) {

  .call <- match.call()

  .tObserved = tObserved;
  .tBoot     = tBoot;
  .conf      = conf;
  .R         = length(tBoot);

  .mean           <- mean(.tBoot, na.rm = T);
  .median         <- median(.tBoot, na.rm = T);
  .bias           <- .mean - .tObserved;
  .se             <- sqrt(var(.tBoot, na.rm = T));
  .biasCorrection <- sum(.tBoot <= .tObserved)/(length(.tBoot) + 1);
  .z0             <- qnorm(.biasCorrection);
  .alphaL          <- (1 - .conf)/2;
  .alphaU          <- (1 + .conf)/2;
  .zalphaL         <- qnorm(.alphaL);
  .zalphaU         <- qnorm(.alphaU);
  .adjAlphaL       <- pnorm(2 * .z0 + .zalphaL);
  .adjAlphaU       <- pnorm(2 * .z0 + .zalphaU);
  .lowerSequence   <- trunc(.adjAlphaL * (.R + 1));
  .upperSequence   <- trunc(.adjAlphaU * (.R + 1));

  .lowerCI         <- ifelse(.lowerSequence == 0, min(.tBoot), sort(.tBoot)[.lowerSequence]);
  .upperCI         <- ifelse(.upperSequence == .R + 1, max(.tBoot), sort(.tBoot)[.upperSequence]);

  names(.R)                    <- c('Number of Boot Replications');
  names(.tObserved)            <- c('Original Statistic');
  names(.mean)                 <- c('Mean of Boot Replications');
  names(.median)               <- c('Median of Boot Replications');
  names(.bias)                 <- c('Bias of Original Statistic');
  names(.biasCorrection)        <- c('Bias Correction: (Boot Replicates <= Original Statistic)/(Total Replicates)');
  names(.z0)					<- c('z0');
  names(.alphaL)					<- c('alphaL');
  names(.alphaU)					<- c('alphaU');
  names(.zalphaL)				<- c('zalphaL');
  names(.zalphaU)				<- c('zalphaU');
  names(.adjAlphaL)				<- c('Adjusted Alpha Lower');
  names(.adjAlphaU)				<- c('Adjusted Alpha Upper');
  names(.lowerSequence)         <- c('Lower Sequence Number');
  names(.upperSequence)         <- c('Upper Sequence Number');

  names(.se)                   <- c('Standard Error of Boot Statistic');
  names(.lowerCI)                <- paste((100 * .alphaL), '% Bias Corrected Lower Bound', sep = '')
  names(.upperCI)                <- paste((100 * .alphaU), '% Bias Corrected Upper Bound', sep = '')

  return(c(.call, .tObserved, .mean, .median, .bias, .se, .biasCorrection,
    .z0, .alphaL, .alphaU, .zalphaL, .zalphaU, .adjAlphaL, .adjAlphaU, .lowerSequence, .upperSequence, .lowerCI, .upperCI, .R))
}
