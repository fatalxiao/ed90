#' @title  Estimate ED50
#' @description Estimate 50 percent effective dose using different methods.
#' @param doseSequence A sequence of doses given in order
#' @param doseResponse A sequence of response results shown in order
#' @param confidence The confidence level of interval estimate
#' @param method The method used to estimate ED50, there are five methods here, respectively
#' Dixon-Mood, Choi (Choi's Original Turning Point), ModTurPoint (Modified Turning Point),
#' Logistic (Logistic Regression) and Isotonic (Isotonic Regression). The defaut is Dixon-Mood.
#' @param tpCiScale The scale level to enlarge the confidence interval estimated by Modified
#' Turning Point Method. The default value is \code{2.4/qnorm(0.975)}.
#' @param boot.n The number of boot process if Logistic method is chosen to estimate ED50.
#' @import stats
#' @import boot
#' @export
#' @return A list of estimation result consisting of method of estimation, ED50 estimate,
#' standard error of ED50 estimate, confidence level and estimate of confidence interval.
#' @references
#' Dixon, W. J., & Mood, A. M. (1948). A method for obtaining and analyzing
#' sensitivity data. Publications of the American Statistical Association, 43(241), 109-126.
#' Choi, S. C. (1990). Interval estimation of the ld50based on an up-and-down experiment.
#' Biometrics, 46(2), 485-492.
#' Pace, N. L., & Stylianou, M. P. (2007). Advances in and limitations of up-and-down
#' methodology: a precis of clinical use, study design, and dose estimation in anesthesia
#' research. Anesthesiology, 107(1), 144-52.
#' @examples
#' library(ed50)
#' estimate(groupS$doseSequence, groupS$responseSequence, method = 'Dixon-Mood')
#' estimate(groupS$doseSequence, groupS$responseSequence, method = 'Logistic', boot.n = 1000)

estimate <- function(doseSequence,
                      doseResponse,
                      confidence = .95,
                      method = c('Dixon-Mood', 'Choi', 'ModTurPoint', 'Logistic', 'Isotonic'),
                      tpCiScale = 2.4/qnorm(0.975),
                      boot.n = 10000)
{
   # First check out the type of data input
   # Whether the dose data or response data is a numeric vector...
   if(!is.vector(doseSequence))
     return(warning('The dose data input is not a vector!'))
   if(!is.numeric(doseSequence))
     return(warning('The type of dose data input is not numeric!'))
   if(!is.vector(doseResponse))
     return(warning('The response data input is not a vector!'))
   if(!is.numeric(doseResponse))
     return(warning('The type of response data input is not numeric!'))

   # Check if the lengths of dose data and response data input are the same
   doseLength <- length(doseSequence)
   responseLength <- length(doseResponse)
   if(doseLength != responseLength)
     return(warning('Dose data and response data input are of different length!'))

   # Get the method
  method <- tryCatch(match.arg(method), error = function(e) 'error')
  if(method == 'error')
    return(warning('The parameter method should be one of "Dixon-Mood", "ModTurPoint", "Logistic" and "Isotonic"!'))

  # Prepare for estimation
  # Compute the fixed dose step and check if there is a mistake.
  doseDiff <- round(diff(doseSequence), 10)
  doseStep <- unique(abs(doseDiff))
  if(length(doseStep) > 1)
    return(warning('The dose step is not fixed!'))
  if(length(doseStep) == 0)
    return(warning('There is only one dose sample input!'))
  if(length(doseStep) == 1 & doseStep == 0)
    return(warning('The dose step is set to be zero!'))

  # Estimate ED50 using Dixon-Mood method
  if(method == 'Dixon-Mood')
  {
    # Get the class of the target dose, failure one or success
    tmp1 <- table(doseResponse)
    tmp2 <- as.numeric(names(tmp1)[which.min(tmp1)])

    # Get the specific target dose data
    doseTarget <- doseSequence[doseResponse == tmp2]

    # Get the sample size of all levels of target doses and rank
    tmp3 <- table(doseTarget)
    tmp4 <- as.numeric(names(tmp3))
    tmp5 <- seq(min(tmp4), max(tmp4), doseStep)
    tmp6 <- setdiff(tmp5, tmp4)
    tmp7 <- length(tmp6)
    if(tmp7 != 0)
    {
      tmp8 <- rep(0, tmp7)
      names(tmp8) <- tmp6
      tmp3 <- c(tmp3, tmp8)
    }
    tmp3 <- tmp3[as.character(tmp5)]

    # Calculate the ED50 estimation using Dixon-Mood method
    n <- tmp3
    N <- sum(n)
    i <- seq_along(n) - 1
    A <- sum(i * n)
    B <- sum(i^2 * n)
    m <- min(tmp4) + doseStep * (A/N + 0.5 * (-1)^tmp2)

    # Calculate the standard error of ED50 estimate
    # But parameter G should be first determined
    s <- 1.62 * doseStep * ((N * B - A^2) / (N^2) + 0.029)
    ratio <- doseStep / s
    gTableOrigin <- data.frame(Ratio = seq(0.2, 5.0, 0.1),
                               G1 = c(0.92, 0.925, 0.94, 0.95, 0.96, 0.975, 0.98, 0.99,
                                      1, 1.01, 1.025, 1.04, 1.05, 1.07, 1.08, 1.1, 1.11,
                                      1.12, 1.13, 1.155, 1.17, 1.18, 1.195, 1.2, 1.2, 1.205,
                                      1.21, 1.215, 1.22, 1.222, 1.224, 1.226, 1.228, 1.23,
                                      1.231, 1.232, 1.233, 1.234, 1.235, 1.236, 1.237, 1.238,
                                      1.239, 1.24, 1.241, 1.242, 1.243, 1.244, 1.245),
                               G2 = c(0.92, 0.925, 0.94, 0.95, 0.96, 0.975, 0.98, 0.99, 1,
                                      1.01, 1.025, 1.04, 1.05, 1.07, 1.08, 1.09, 1.1, 1.11,
                                      1.12, 1.15, 1.175, 1.2, 1.22, 1.245, 1.285, 1.305, 1.33,
                                      1.37, 1.4, 1.45, 1.495, 1.53, 1.58, 1.63, 1.7, 1.75, 1.81,
                                      1.895, 1.95, 2.015, 2.1, 2.2, 2.29, 2.39, 2.49, 2.6, 2.75,
                                      2.91, 3.15))
    if(ratio < 0.2)
      return(warning('The dose step might be set too narrow!'))
    if(ratio > 5)
      return(warning('The dose step might be set too wide!'))
    if((ratio >= 0.2 & ratio <= 1.6 & !(m %in% tmp4)) | (m %in% tmp4))
    {
      mode <- loess(formula = G1 ~ Ratio, data = gTableOrigin)
      G <- as.vector(predict(mode, newdata = data.frame(Ratio = ratio)))
    }
    if(ratio > 1.6 & ratio <= 5 & !(m %in% tmp4))
    {
      mode <- loess(formula = G2 ~ Ratio, data = gTableOrigin[gTableOrigin$Ratio >= 1.6, ])
      G <- as.vector(predict(mode, newdata = data.frame(Ratio = ratio)))
    }
    sm <- G * s / sqrt(N)

    # Calculate the boundary of confidence interval
    lb <- m - qnorm(0.5 + 0.5 * confidence) * sm
    ub <- m + qnorm(0.5 + 0.5 * confidence) * sm

    # Summarise the whole result
    ans <- list('Method of Estimation'= 'Dixon-Mood',
                'Estimate of ED50' = m,
                'Standard Error of Estimate' = sm,
                'Value of Parameter G' = G,
                'Confidence Level' = paste0(100 * confidence, '%'),
                'Lower Bound' = lb,
                'Upper Bound' = ub)
  }

  # Estimate ED50 using modified turning point method
  if(method == 'ModTurPoint' | method == 'Choi')
  {
    # Tell all the turning points and compute the Ws
    tmp1 <- diff(doseDiff)
    tmp2 <- which(abs(tmp1) == (2 * doseStep))
    tmp3 <- doseResponse[responseLength] == doseResponse[responseLength-1]
    if((length(tmp1) == 0 | length(tmp2) == 0) & tmp3)
      return(warning('The dose data input has no turning points!'))
    doseW <- (doseSequence[tmp2] + doseSequence[tmp2 + 1]) / 2
    if(!tmp3)
    {
      addW <- (doseSequence[doseLength] + doseSequence[doseLength-1]) / 2
      doseW <- c(doseW, addW)
    }

    # Compute the number of Ws
    nW <- length(doseW)

    # Calculate ED50 estimate
    if(nW < 3)
    {
      return(warning(paste('There are only a few turning points:', doseW)))
    }
    meanW <- mean(doseW[-1])

    # define a, b, c
    if(method == 'ModTurPoint')
    {
      cenW <- doseW - meanW
    } else {
      cenW <- doseW
    }
    a <- cenW[2]^2 + cenW[nW]^2
    b <- sum(cenW[-c(1, nW)] * cenW[-c(1, 2)])
    c <- sum(cenW[-c(1, 2, nW)]^2)

    # Solve the systerm of equations
    rhoHatFunc <- function(rho)
    {
      (b - rho * c) * (1 - rho^2) - rho * (a - 2 * rho * b + (1 + rho^2) * c) / (nW - 1)
    }
    rhoHat <- tryCatch(uniroot(rhoHatFunc,
                               interval = c(-1+10^(-5), 1-10^(-5)),
                               tol = 1e-9)$root, error = function(e) 'error')
    if(rhoHat == 'error') return(warning('Problem occured in solving the equations!\nMaybe the number of turning points is not enough.'))
    rhoHat <- max(abs(rhoHat))
    sigHat2 <- (a - 2 * rhoHat * b + (1 + rhoHat^ 2) * c) / (nW - 1)

    # Compute standard error of ED50 estimate
    i <- 1:(nW - 2)
    sd <- ((sigHat2 / ((nW - 1) * (1 - rhoHat^2))) *
             (1 + sum((2 * (nW - i - 1) * rhoHat^i) / (nW - 1))))^(0.5)

    # Summarise ED50 estimate and its confidence interval
    zAlpha <- qnorm(0.5 + 0.5 * confidence) * tpCiScale
    lb <- meanW - zAlpha * sd
    ub <- meanW + zAlpha * sd
    ans <- list('Method of Estimation'= ifelse(method == 'ModTurPoint',
                                               'Modified Turning Point',
                                               "Choi's Original Turning Point"),
                'Number of turning points' = nW,
                'Estimate of ED50' = meanW,
                'Standard Error of Estimate' = sd,
                'Confidence Level' = paste0(100 * confidence, '%'),
                'Lower Bound' = lb,
                'Upper Bound' = ub)
  }

  # Estimate ED50 using logistic regression method
  if(method == 'Logistic')
  {
    # Write data into a data frame
    dataFrame <- data.frame(doseSequence = doseSequence,
                            doseResponse = doseResponse)

    # Calculate the ED50 estimate
    mode <- glm(formula = doseResponse ~ doseSequence,
                data = dataFrame,
                family = binomial)
    betaHat <- as.vector(coef(mode))
    muHat <- - betaHat[1] / betaHat[2]

    # Use oversampling to estimate ci
    # The boot function
    boot.fn <- function(data, index)
    {
      mode <- glm(formula = doseResponse ~ doseSequence,
                  data = data,
                  subset = index,
                  family = binomial)
      betaHat <- as.vector(coef(mode))
      muHat <- - betaHat[1] / betaHat[2]
      return(muHat)
    }
    set.seed(1)
    bootRes <- suppressWarnings(boot(data = dataFrame,
                                     statistic = boot.fn,
                                     R = boot.n))
    bootRes <- as.vector(bootRes$t)
    sd <- sqrt(var(bootRes))
    lb <- as.vector(quantile(bootRes, probs = (1 - confidence) / 2, na.rm = TRUE))
    ub <- as.vector(quantile(bootRes, probs = .5 + confidence / 2, na.rm = TRUE))

    # Summarise the whole result
    ans <- list('Method of Estimation'= 'Logistic Regression',
                'Estimate of ED50' = muHat,
                'Standard Error of Estimate' = sd,
                'Confidence Level' = paste0(100 * confidence, '%'),
                'Lower Bound' = lb,
                'Upper Bound' = ub)
  }

  # Estimate ED50 using isotonic regression method
  # Just use the programme written by Professor Pace
  if(method == 'Isotonic')
  {
    # Create a data frame
    dataFrame <- data.frame(doseSequence = doseSequence,
                        responseSequence = doseResponse)

    # Change the data into the special form using PAVA algorithm
    pavaData <- preparePava(dataFrame)

    # This the boot function
    bootResult <- boot(data = dataFrame,
                  statistic = bootIsotonicRegression,
                          R = boot.n,
                        sim = 'parametric',
                    ran.gen = bootIsotonicResample,
                        mle = list(baselinePava = pavaData,
                                      firstDose = doseSequence[1],
                              PROBABILITY.GAMMA = 0.5),
               baselinePava = pavaData,
          PROBABILITY.GAMMA = 0.5)

    # Get the prediction result of the confidence interval
    prediction <- bootBC.ci(tObserved = bootResult$t0[3],
                                tBoot = bootResult$t[, 3],
                                 conf = confidence)

    # Clean the prediction result
    predictionLength <- length(prediction)
    ans <- list('Method of Estimation' = 'Isotonic',
                'Estimate of ED50' = prediction$`Mean of Boot Replications`,
                'Standard Error of Estimate' = prediction$`Standard Error of Boot Statistic`,
                'Confidence Level' = paste0(100 * confidence, '%'),
                'Lower Bound' = prediction[[predictionLength - 2]],
                'Upper Bound' = prediction[[predictionLength - 1]])
  }

  return(ans)
 }
