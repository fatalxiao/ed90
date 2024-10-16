#' @title Compare ED50 Estimation of Independent Two-sample Case
#' @description Test the statistical  difference of two independent estimation results of ED50.
#' @param group1 A list object of ED50 estimation.
#' @param group2  Another list object of ED50 estimation to be compared with.
#' @param alpha The significant level of test. 0.05 is the defaut value.
#' @import stats
#' @export
#' @return The difference between two groups of ED50 estimation in terms of statistical
#' significance.
#' @references Noguchi, K., & Marmolejo-Ramos, F. (2016). Assessing equality of means using the
#' overlap of range-preserving confidence intervals. American Statistician, 70(4), 325-334.
#' @examples
#' library(ed50)
#' ans1 <- estimate(groupS$doseSequence, groupS$responseSequence, method = 'ModTurPoint')
#' ans2 <- estimate(groupSN$doseSequence, groupSN$responseSequence, method = 'Dixon-Mood')
#' compare(ans1, ans2)

compare <- function(group1, group2, alpha = .05)
{
  # Take out the result of ED50 estimates of the two methods
  mean1 <- group1$`Estimate of ED50`
  mean2 <- group2$`Estimate of ED50`
  # Take out the estimates of standard error of ED50 estimates
  std1 <- group1$`Standard Error of Estimate`
  std2 <- group2$`Standard Error of Estimate`

  # Caculate the test statistic and p value
  statistic <- (mean1 - mean2) * ((std1^2 + std2^2)^(-0.5))
  pValue <- 2 * (1 - pnorm(abs(statistic)))
  # Print the results
  if(pValue < alpha)
    cat('\nTest Result \n----------- ',
        '\nThere is significant difference between ED50 estiation of the two groups',
        '\nwith p-value =', pValue, 'and significance level =', alpha, '\n\t')
  if(pValue >= alpha)
    cat('\nTest Result \n----------- ',
        '\nThere is no significant difference between ED50 estiation of the two groups',
        '\nwith p-value =', pValue, 'and significance level =', alpha, '\n\t')
}
