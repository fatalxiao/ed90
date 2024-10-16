#' @title Isotonic Regression Function
#' @description Function of isotonic regression.
#' @param data the same dataframe called by the boot function.
#' @param PROBABILITY.GAMMA the target effect probability in the BCD experiment; default = 0.5 and need not be specified.
#' @param baselinePava the dataframe prepared by the function preparePava.
#' @import stats
#' @export
#' @examples
#' library(ed50)
#' pavaData <- preparePava(groupS)
#' bootIsotonicRegression(data = groupS, PROBABILITY.GAMMA = 0.5, baselinePava = pavaData)

bootIsotonicRegression <- function(data, PROBABILITY.GAMMA = 0.5, baselinePava) {

# estimates muHat1, muHat2, muHat3, pmuHat1, pmuHat2, pmuHat3 of original data

  .pavaProbabilityBaseline  <- baselinePava$pavaProbability;
  .nDosesBaseline			<- baselinePava$nDoses;

if (min(.pavaProbabilityBaseline, na.rm = T) > PROBABILITY.GAMMA) {
  .muHat1Baseline      <- min(.nDosesBaseline, na.rm = T);
  .muHat1NextBaseline  <- min(.nDosesBaseline, na.rm = T);
  .muHat2Baseline      <- min(.nDosesBaseline, na.rm = T);
  .muHat3Baseline      <- min(.nDosesBaseline, na.rm = T);
  .pmuHat1Baseline     <- .pavaProbabilityBaseline[.nDosesBaseline = .muHat1Baseline];
  .pmuHat2Baseline     <- .pavaProbabilityBaseline[.nDosesBaseline = .muHat2Baseline];
  .pmuHat3Baseline     <- .pavaProbabilityBaseline[.nDosesBaseline = .muHat3Baseline];

} else if (max(.pavaProbabilityBaseline, na.rm = T) <= PROBABILITY.GAMMA) {
  .muHat1Baseline      <- max(.nDosesBaseline, na.rm = T);
  .muHat1NextBaseline  <- max(.nDosesBaseline, na.rm = T);
  .muHat2Baseline      <- max(.nDosesBaseline, na.rm = T);
  .muHat3Baseline      <- max(.nDosesBaseline, na.rm = T);
  .pmuHat1Baseline     <- .pavaProbabilityBaseline[.nDosesBaseline = .muHat1Baseline];
  .pmuHat2Baseline     <- .pavaProbabilityBaseline[.nDosesBaseline = .muHat2Baseline];
  .pmuHat3Baseline     <- .pavaProbabilityBaseline[.nDosesBaseline = .muHat3Baseline];

} else {
  .muHat1Baseline      <-
    max(.nDosesBaseline[.pavaProbabilityBaseline <= PROBABILITY.GAMMA], na.rm = T);
  .muHat1NextBaseline  <-
    min(.nDosesBaseline[.nDosesBaseline > .muHat1Baseline], na.rm = T);
  .pmuHat1Baseline     <-
	.pavaProbabilityBaseline[.nDosesBaseline == .muHat1Baseline];
  .pmuHat1NextBaseline <-
    .pavaProbabilityBaseline[.nDosesBaseline == .muHat1NextBaseline];

# code to handle the floating point problem of PAVA estimators equidistant from PROBABILITY.GAMMA

  isPmuHat2Equal       <-
     identical(all.equal(abs(.pmuHat1Baseline - PROBABILITY.GAMMA), abs(.pmuHat1NextBaseline - PROBABILITY.GAMMA)), T);
  isPmuHat2Unequal     <-
    !identical(all.equal(abs(.pmuHat1Baseline - PROBABILITY.GAMMA), abs(.pmuHat1NextBaseline - PROBABILITY.GAMMA)), T);

  if  (isPmuHat2Equal) {
    .muHat2Baseline    <-  .muHat1Baseline;
  } else if (isPmuHat2Unequal && abs(.pmuHat1Baseline - PROBABILITY.GAMMA) < abs(.pmuHat1NextBaseline - PROBABILITY.GAMMA)) {
    .muHat2Baseline    <-  .muHat1Baseline;
  } else {
    .muHat2Baseline    <-  .muHat1NextBaseline;
  }

  .pmuHat2Baseline     <-
    .pavaProbabilityBaseline[.nDosesBaseline == .muHat2Baseline];

  .muHat3Baseline      <-
    ((PROBABILITY.GAMMA - .pmuHat1Baseline) / (.pmuHat1NextBaseline - .pmuHat1Baseline)) * (.muHat1NextBaseline - .muHat1Baseline) + .muHat1Baseline;
  .pmuHat3Baseline     <- PROBABILITY.GAMMA

 }

# start of code to get PAVA estimators for each resampled data set

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

  for (i in 1:k) {
    .sum1 <- 0;
    .sum2 <- 0;
    for(j in i:k) {
      .sum1 <- .sum1 + .nEvents[j];
      .sum2 <- .sum2 + .nTrials[j];
      pp[i, j] <- .sum1/.sum2;
    }
    .tempMatrix2 <- as.matrix(pp[(1:i), (i:k)]);
    p[i] <- ifelse(i > 1, max(apply(.tempMatrix2, 1, min)), min(.tempMatrix2));
  }

  .tempDataFrame3 <- data.frame(
    .naiveProbability = .nEvents/.nTrials,
    .pavaProbability  = p,
    .nEvents          = .nEvents,
    .nTrials          = .nTrials,
    .nDoses           = .nDoses
  )

  .pavaProbability  <- .tempDataFrame3$.pavaProbability;
  .nDoses			<- .tempDataFrame3$.nDoses;

if (min(.pavaProbability, na.rm = T) > PROBABILITY.GAMMA) {
  .muHat1      <- min(.nDoses, na.rm = T);
  .muHat2      <- min(.nDoses, na.rm = T);
  .muHat3      <- min(.nDoses, na.rm = T);

  .pmuHat1     <- .pavaProbability[.nDoses == .muHat1Baseline];
  .pmuHat2     <- .pavaProbability[.nDoses == .muHat2Baseline];
  .pmuHat1Next <- .pavaProbability[.nDoses == .muHat1NextBaseline];
  .pmuHat3     <-
    ((.muHat3Baseline - .muHat1Baseline) / (.muHat1NextBaseline - .muHat1Baseline)) * (.pmuHat1Next - .pmuHat1) + .pmuHat1;

} else if (max(.pavaProbability, na.rm = T) <= PROBABILITY.GAMMA) {
  .muHat1      <- max(.nDoses, na.rm = T);
  .muHat2      <- max(.nDoses, na.rm = T);
  .muHat3      <- max(.nDoses, na.rm = T);

  .pmuHat1     <- .pavaProbability[.nDoses == .muHat1Baseline];
  .pmuHat2     <- .pavaProbability[.nDoses == .muHat2Baseline];
  .pmuHat1Next <- .pavaProbability[.nDoses == .muHat1NextBaseline];
  .pmuHat3     <-
    ((.muHat3Baseline - .muHat1Baseline) / (.muHat1NextBaseline - .muHat1Baseline)) * (.pmuHat1Next - .pmuHat1) + .pmuHat1;

} else {
  .muHat1			   <- max(.nDoses[.pavaProbability <= PROBABILITY.GAMMA], na.rm = T);
  .muHat1Next          <- min(.nDoses[.nDoses > .muHat1], na.rm = T);

  .pmuHat1     		   <- .pavaProbability[.nDoses == .muHat1];
  .pmuHat1Next  	   <- .pavaProbability[.nDoses == .muHat1Next];

  isPmuHat2Equal       <-  identical(all.equal(abs(.pmuHat1 - PROBABILITY.GAMMA), abs(.pmuHat1Next - PROBABILITY.GAMMA)), T);
  isPmuHat2Unequal     <- !identical(all.equal(abs(.pmuHat1 - PROBABILITY.GAMMA), abs(.pmuHat1Next - PROBABILITY.GAMMA)), T);

  if  (isPmuHat2Equal) {
    .muHat2            <-  .muHat1;
  } else if (isPmuHat2Unequal && abs(.pmuHat1 - PROBABILITY.GAMMA) < abs(.pmuHat1Next - PROBABILITY.GAMMA)) {
    .muHat2            <-  .muHat1;
  } else {
    .muHat2            <-  .muHat1Next;
  }

  .muHat3			<- ((PROBABILITY.GAMMA - .pmuHat1) / (.pmuHat1Next - .pmuHat1)) * (.muHat1Next - .muHat1) + .muHat1;

  .pmuHat1			<- .pavaProbability[.nDoses == .muHat1Baseline];
  .pmuHat2			<- .pavaProbability[.nDoses == .muHat2Baseline];
  .pmuHat1Next		<- .pavaProbability[.nDoses == .muHat1NextBaseline];
  .pmuHat3          <-
    ((.muHat3Baseline - .muHat1Baseline) / (.muHat1NextBaseline - .muHat1Baseline)) * (.pmuHat1Next - .pmuHat1) + .pmuHat1;
}

# if missing, the value of the estimates of muHati of the resampled data are imputed to the value of the original data
# if missing, the value of the estimates of pmuHati of the resampled data are imputed to 0 or 1

  estimates <- c(
    muHat1  = ifelse(is.na(ifelse(exists('.muHat1'),  .muHat1,  NA)), .muHat1Baseline,  .muHat1),
    muHat2  = ifelse(is.na(ifelse(exists('.muHat2'),  .muHat2,  NA)), .muHat2Baseline,  .muHat2),
    muHat3  = ifelse(is.na(ifelse(exists('.muHat3'),  .muHat3,  NA)), .muHat3Baseline,  .muHat3),
    pmuHat1 = ifelse(is.na(ifelse(exists('.pmuHat1'), .pmuHat1, NA)), ifelse(min(.nDoses, na.rm = T) > .muHat1Baseline, 0, 1), .pmuHat1),
    pmuHat2 = ifelse(is.na(ifelse(exists('.pmuHat2'), .pmuHat2, NA)), ifelse(min(.nDoses, na.rm = T) > .muHat2Baseline, 0, 1), .pmuHat2),
    pmuHat3 = ifelse(is.na(ifelse(exists('.pmuHat3'), .pmuHat3, NA)), ifelse(min(.nDoses, na.rm = T) > .muHat1Baseline, 0, 1), .pmuHat3)
  	);

}
