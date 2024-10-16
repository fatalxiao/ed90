# 创建数据框示例
groupS <- read.csv("./groupS.csv", 1, encoding='UTF-8')
groupS

doses <- groupS$doseSequence
responses <- groupS$responseSequence

data <- data.frame(doses, responses)

# 拟合Logistic模型
logistic_model <- glm(responses ~ doses, data = data, family = binomial)

# 定义预测函数
ed90_predict_logistic <- function(model, p = 0.9) {
  q <- qlogis(p)
  b <- coef(model)
  ed90 <- (q - b[1]) / b[2]
  return(ed90)
}

# 计算ED90
ed90_value_logistic <- ed90_predict_logistic(logistic_model)
print(ed90_value_logistic)
