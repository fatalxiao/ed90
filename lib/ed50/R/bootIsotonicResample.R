#' @title The resample function of isotonic regression
#' @description The function is designed as an argument for the boot function of the Canty Bootstrap package.
#' @param data Original experiment data.
#' @param mle A list of additional arguments to be used by bootIsotonicResample.
#' @import stats
#' @export
#' @examples
#' library(ed50)
#' pavaData <- preparePava(groupS)
#' bootIsotonicResample(data = groupS,
#'                       mle = list(baselinePava = pavaData,
#'                                     firstDose = 2.5,
#'                             PROBABILITY.GAMMA = 0.5))

bootIsotonicResample <- function (data, mle) {

# labels the arguments to bootIsotonicResample

  .responseSequence		<- data$responseSequence;
  .doseSequence			<- data$doseSequence;
  .sequenceLength		<- length(.responseSequence)
  .shortSequenceLength  <- .sequenceLength - 1;
  .pavaProbability		<- mle$baselinePava$pavaProbability;
  .nDoses				<- mle$baselinePava$nDoses;
  .firstDose			<- mle$firstDose;
  PROBABILITY.GAMMA		<- mle$PROBABILITY.GAMMA;
  isGAMMALow			<- PROBABILITY.GAMMA < 0.5;
  isGAMMAMid			<- PROBABILITY.GAMMA == 0.5;
  isGAMMAHigh			<- PROBABILITY.GAMMA > 0.5;
  PROBABILITY.BETA		<-
  {
    if (PROBABILITY.GAMMA < 0.5)	{
	  PROBABILITY.BETA <- PROBABILITY.GAMMA/(1 -PROBABILITY.GAMMA);
	  } else if (PROBABILITY.GAMMA == 0.5)	{
	  PROBABILITY.BETA <- 1;
	  } else
	  PROBABILITY.BETA <- (1 - PROBABILITY.GAMMA)/PROBABILITY.GAMMA;
  }

# samples from uniform distribution to set probabilistically the response to the starting dose

  {
    if (runif(1) <= .pavaProbability[.nDoses == .firstDose])	{
	  .firstResponse <- 1;
	  } else													{
	  .firstResponse <- 0;
	  }
  }

# initializes resample vectors

  .resampleDoseSequence                <-
    c(.firstDose, rep(0, times = .shortSequenceLength));
  .resampleResponseSequence            <-
    c(.firstResponse, rep(0, times = .shortSequenceLength));

# the loop to generate the next dose including boundary conditions and probabilistic coin flip

  for (i in 1:.shortSequenceLength) {
    isResponsePositive  <- .resampleResponseSequence[i] == 1;
    isResponseNegative  <- !isResponsePositive;
    isDoseMinimum		<- .resampleDoseSequence[i] == min(.nDoses);
    isDoseMaximum		<- .resampleDoseSequence[i] == max(.nDoses);
    isDoseBetween		<- .resampleDoseSequence[i] != min(.nDoses) &&
      .resampleDoseSequence[i] != max(.nDoses);
    isBeta				<- runif(1) <= PROBABILITY.BETA;
    isNotBeta			<- !isBeta;

  {
    if (isGAMMALow && isDoseMinimum && isResponsePositive)							{
      .resampleDoseSequence[i + 1] <- .resampleDoseSequence[i];
	  } else if (isGAMMALow && isDoseMinimum && isResponseNegative && isNotBeta)	{
      .resampleDoseSequence[i + 1] <- .resampleDoseSequence[i];
	  } else if (isGAMMALow && isDoseMinimum && isResponseNegative && isBeta)		{
      .resampleDoseSequence[i + 1] <- .nDoses[match(.resampleDoseSequence[i], .nDoses) + 1];
	  } else if (isGAMMALow && isDoseBetween && isResponsePositive)					{
      .resampleDoseSequence[i + 1] <- .nDoses[match(.resampleDoseSequence[i], .nDoses) - 1];
	  } else if (isGAMMALow && isDoseBetween && isResponseNegative && isNotBeta)	{
	  .resampleDoseSequence[i + 1] <- .resampleDoseSequence[i];
	  } else if (isGAMMALow && isDoseBetween && isResponseNegative && isBeta)		{
      .resampleDoseSequence[i + 1] <- .nDoses[match(.resampleDoseSequence[i], .nDoses) + 1];
      } else if (isGAMMALow && isDoseMaximum && isResponseNegative)					{
	  .resampleDoseSequence[i + 1] <- .resampleDoseSequence[i];
      } else if (isGAMMALow && isDoseMaximum && isResponsePositive)					{
      .resampleDoseSequence[i + 1] <- .nDoses[match(.resampleDoseSequence[i], .nDoses) - 1];
      } else if (isGAMMAMid && isDoseMinimum && isResponsePositive)					{
	  .resampleDoseSequence[i + 1] <- .resampleDoseSequence[i];
      } else if (isGAMMAMid && isDoseMinimum && isResponseNegative)					{
      .resampleDoseSequence[i + 1] <- .nDoses[match(.resampleDoseSequence[i], .nDoses) + 1];
      } else if (isGAMMAMid && isDoseBetween && isResponsePositive)					{
      .resampleDoseSequence[i + 1] <- .nDoses[match(.resampleDoseSequence[i], .nDoses) - 1];
      } else if (isGAMMAMid && isDoseBetween && isResponseNegative)					{
      .resampleDoseSequence[i + 1] <- .nDoses[match(.resampleDoseSequence[i], .nDoses) + 1];
      } else if (isGAMMAMid && isDoseMaximum && isResponsePositive)					{
      .resampleDoseSequence[i + 1] <- .nDoses[match(.resampleDoseSequence[i], .nDoses) - 1];
      } else if (isGAMMAMid && isDoseMaximum && isResponseNegative)					{
	  .resampleDoseSequence[i + 1] <- .resampleDoseSequence[i];
      } else if (isGAMMAHigh && isDoseMinimum && isResponsePositive)				{
	  .resampleDoseSequence[i + 1] <- .resampleDoseSequence[i];
      } else if (isGAMMAHigh && isDoseMinimum && isResponseNegative)				{
      .resampleDoseSequence[i + 1] <- .nDoses[match(.resampleDoseSequence[i], .nDoses) + 1];
      } else if (isGAMMAHigh && isDoseBetween && isResponsePositive && isBeta)		{
      .resampleDoseSequence[i + 1] <- .nDoses[match(.resampleDoseSequence[i], .nDoses) - 1];
      } else if (isGAMMAHigh && isDoseBetween && isResponsePositive && isNotBeta)	{
	  .resampleDoseSequence[i + 1] <- .resampleDoseSequence[i];
      } else if (isGAMMAHigh && isDoseBetween && isResponseNegative)				{
      .resampleDoseSequence[i + 1] <- .nDoses[match(.resampleDoseSequence[i], .nDoses) + 1];
      } else if (isGAMMAHigh && isDoseMaximum && isResponsePositive && isBeta)		{
      .resampleDoseSequence[i + 1] <- .nDoses[match(.resampleDoseSequence[i], .nDoses) - 1];
      } else if (isGAMMAHigh && isDoseMaximum && isResponsePositive && isNotBeta)	{
	  .resampleDoseSequence[i + 1] <- .resampleDoseSequence[i];
      } else if (isGAMMAHigh && isDoseMaximum && isResponseNegative)						{
	  .resampleDoseSequence[i + 1] <- .resampleDoseSequence[i];
	  }
  }

# samples from uniform distribution to set probabilistically the response to the next dose

  {
    if (runif(1) <= .pavaProbability[.nDoses == .resampleDoseSequence[i + 1]])  {
	  .resampleResponseSequence[i + 1] <- 1;
	  } else															{
	  .resampleResponseSequence[i + 1] <- 0;
	  }
  }

  }

  estimates <- data.frame(
    responseSequence = .resampleResponseSequence,
    doseSequence     = .resampleDoseSequence);


}
