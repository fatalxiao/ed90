#' @title Generate Simulation Data of Up-and-Down Experiment
#' @description The function is used to generate simulation data of up-and-down
#' experiment, and provide three cases that tolerance distribution obeys normal,
#' triangle or chi-square distribution.
#' @param number The number of experiments in a trail.
#' @param useTurPoint A logical value indicating whether the parameter \code{number} refers to the
#' amount of turning points. The default value is \code{FALSE}.
#' @param start The first dose level given in this trail.
#' @param doseStep A fix value  that represents the difference between two adjacent dose levels.
#' @param distribution The tolerance distribution, including normal, triangle and chi-square
#' distribution, and the default distribution is N(0, 1).
#' @param normalMean Parameter mean of normal distribution, the default value is 0.
#' @param normalStd Parameter std of normal distribution, the default value is 1.
#' @param triMean Parameter mean of triangle distribution, the default value is 0.
#' @param triWidth Parameter width of triangle distribution, the default value is 2.
#' @param chiDegree Parameter degree of freedom of chi-square distribution, the default value is 1.
#' @import stats
#' @import utils
#' @export
#' @return A data frame.
#' @examples
#' library(ed50)
#' generateData(number = 20, start = 2, doseStep = 0.2, distribution = 'Normal')
#' generateData(number = 40, start = 2, doseStep = 0.2, distribution = 'Chi-square')

generateData <- function(number,
                         useTurPoint = FALSE,
                         start,
                         doseStep = 1,
                         distribution = c('Normal', 'Triangle', 'Chi-square'),
                         normalMean = 0,
                         normalStd = 1,
                         triMean = 0,
                         triWidth = 2,
                         chiDegree = 1)
{
  # Get the selected distribution
  distribution <- tryCatch(match.arg(distribution), error = function(e) 'error')
  if(distribution == 'error')
  {
    return(warning('The distribution should be one of "Normal", "Triangle", "Chi-square"!'))
  }

  # Init dose sequence value and dose response label
  doseSequence <- NULL
  doseResponse <- NULL
  dose         <- start
  numTurPoint  <- 0
  numDose      <- 0

  repeat({
    # Generate tolerance data
    if(distribution == "Normal")
    {
      tolerance <- rnorm(1, normalMean, normalStd)
    }
    if(distribution == "Triangle")
    {
      u <- runif(1)
      if(u <= 0.5)
      {
        tolerance <- triMean - (0.5 - sqrt(0.5 * u)) * triWidth
      } else {
        tolerance <- triMean + (0.5 - sqrt(0.5 * (1-u))) * triWidth
      }
    }
    if(distribution == "Chi-square")
    {
      tolerance <- rchisq(1, chiDegree)
    }

    # Genarate dose and response data
    doseSequence <- c(doseSequence, dose)
    if(dose < tolerance)
    {
      dose <- dose + doseStep
      doseResponse <- c(doseResponse, 0)
    } else {
      dose <- dose - doseStep
      doseResponse <- c(doseResponse, 1)
    }

    # Set the condition to stop the function
    if(useTurPoint)
    {
      if(length(unique(tail(doseResponse, 2))) == 2)
      {
        numTurPoint <- numTurPoint + 1
      }
      if(numTurPoint == number) break
    } else {
      numDose <- numDose + 1
      if(numDose == number) break
    }
  })

  return(data.frame(doseSequence = round(doseSequence, 2),
                    doseResponse = doseResponse))
}
