# 安装和加载所需的包
if (!require("ed50")) {
    install.packages("ed50")
}
library(ed50)
if (!require("boot")) {
    install.packages("boot")
}
library(boot)

groupS <- read.csv("./groupS.csv", 1, encoding = 'UTF-8')
groupS

doseSequence <- groupS$doseSequence
doseResponse <- groupS$responseSequence
confidence <- .95
tpCiScale <- 2.4 / qnorm(0.975)
boot.n <- 2000

# Create a data frame
dataFrame <- data.frame(doseSequence = doseSequence,
                        responseSequence = doseResponse)
dataFrame

# Change the data into the special form using PAVA algorithm
pavaData <- preparePava(dataFrame)
pavaData

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
bootResult

# Get the prediction result of the confidence interval
prediction <- bootBC.ci(tObserved = bootResult$t0[3],
                        tBoot = bootResult$t[, 3],
                        conf = confidence)
prediction

# Clean the prediction result
predictionLength <- length(prediction)
ans <- list('Method of Estimation' = 'Isotonic',
            'Estimate of ED50' = prediction$`Mean of Boot Replications`,
            'Standard Error of Estimate' = prediction$`Standard Error of Boot Statistic`,
            'Confidence Level' = paste0(100 * confidence, '%'),
            'Lower Bound' = prediction[[predictionLength - 2]],
            'Upper Bound' = prediction[[predictionLength - 1]])
ans
