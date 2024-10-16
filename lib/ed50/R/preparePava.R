#' @title Covert Data Using PAVA Algorithm
#' @description  Covert data using PAVA algorithm, the result is uesd for isotonic regression estimation.
#' @param data A data frame of dose experiments.
#' @import stats
#' @export
#' @examples
#' library(ed50)
#' preparePava(groupS)
#' preparePava(groupSN)

preparePava <- function(data) {

# tabulate sequential response and dose vector

  .responseSequence <- data$responseSequence;
  .doseSequence     <- data$doseSequence;

	.temp1.df		<- aggregate(.responseSequence, by =  list(.doseSequence), FUN = length);
  .temp1.df   <- data.frame(.nDoses = .temp1.df$Group.1, .nTrials = .temp1.df$x);

  .temp2.df		<- aggregate(.responseSequence, by =  list(.doseSequence), FUN = sum);
  .temp2.df   <- data.frame(.nDoses = .temp2.df$Group.1, .nEvents = .temp2.df$x);

  .temp3.df		<- merge(.temp1.df, .temp2.df, by = '.nDoses');

# apply code found on s-news listserver (Max Kuhn, 10-JUN-1999) to implement PAVA

  .nEvents <- .temp3.df$.nEvents;
  .nTrials <- .temp3.df$.nTrials;
  .nDoses  <- .temp3.df$.nDoses;

  k <- length(.nEvents);
  p <- matrix(0, k, 1);
  pp <- matrix(0, k, k);

  for (i in seq(along = .nEvents)) {
    .sum1 <- 0;
    .sum2 <- 0;
    for(j in i:k) {
      .sum1 <- .sum1 + .nEvents[j];
      .sum2 <- .sum2 + .nTrials[j];
      pp[i, j] <- .sum1/.sum2;
    }
    .tempMatrix <- as.matrix(pp[(1:i), (i:k)]);
    p[i] <- ifelse(i > 1, max(apply(.tempMatrix, 1, min)), min(.tempMatrix));
  }

  .tempDataFrame <- data.frame(
    naiveProbability = .nEvents/.nTrials,
    pavaProbability  = p,
    nEvents          = .nEvents,
    nTrials          = .nTrials,
    nDoses           = .nDoses
  )

  return(pavaDataFrame = .tempDataFrame);
}





