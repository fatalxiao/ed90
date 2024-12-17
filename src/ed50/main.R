# 安装和加载所需的包
if (!require("ed50")) {
  install.packages("ed50")
}
library(ed50)
if (!require("boot")) {
  install.packages("boot")
}
library(boot)

# 剂量
# doseSequence <- c(2.5,2.7,2.9,3.1,3.3,3.1,3.3,3.5,3.3,3.5,3.7,3.5,3.3,3.5,3.7,3.9,3.7,3.5,3.7,3.9,3.7,3.9,3.7,3.5,3.7,3.9,3.7,3.9,3.7,3.5,3.7,3.9,3.7,3.5,3.3,3.5)
# 对应剂量的反应
# responseSequence <- c(0,0,0,0,1,0,0,1,0,0,1,1,0,0,0,1,1,0,0,1,0,1,1,0,0,1,0,1,1,0,0,1,1,1,0,0)
# 创建一个数据框
# groupS <- data.frame(doseSequence = doseSequence, responseSequence = responseSequence)

# 读 CSV 文件
groupS <- read.csv("./loading volume data.csv", 1, encoding='UTF-8')
# print(groupS)

# Change the data into the special form using PAVA algorithm
pavaData <- preparePava(groupS)
pavaData

# confidence <- .95
#
# # This the boot function
# bootResult <- boot(data = groupS,
#               statistic = bootIsotonicRegression,
#                       R = 2000,
#                     sim = "parametric",
#                 ran.gen = bootIsotonicResample,
#                     mle = list(baselinePava = pavaData,
#                                   firstDose = 7,
#                           PROBABILITY.GAMMA = 0.5),
#            baselinePava = pavaData,
#       PROBABILITY.GAMMA = 0.5)
# bootResult
#
# # Get the prediction result of the confidence interval
# prediction <- bootBC.ci(tObserved = bootResult$t0[3],
#                             tBoot = bootResult$t[, 3],
#                              conf = confidence)
# prediction
#
# # Clean the prediction result
# predictionLength <- length(prediction)
# ans <- list('Method of Estimation' = 'Isotonic',
#             'Estimate of ED50' = prediction$`Mean of Boot Replications`,
#             'Standard Error of Estimate' = prediction$`Standard Error of Boot Statistic`,
#             'Confidence Level' = paste0(100 * confidence, '%'),
#             'Lower Bound' = prediction[[predictionLength - 2]],
#             'Upper Bound' = prediction[[predictionLength - 1]])
# ans

