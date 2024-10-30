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
