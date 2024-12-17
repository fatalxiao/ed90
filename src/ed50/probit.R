# 创建数据框示例
groupS <- read.csv("./loading volume data.csv", 1, encoding = 'UTF-8')
groupS

doses <- groupS$doseSequence
responses <- groupS$responseSequence

data <- data.frame(doses, responses)

# 拟合Probit模型
probit_model <- glm(responses ~ doses, data = data, family = binomial(link = "probit"))

# 定义预测函数
ed90_predict <- function(model, p = 0.9) {
    q <- qnorm(p)
    b <- coef(model)
    ed90 <- (q - b[1]) / b[2]
    return(ed90)
}

# 计算ED90
ed90_value <- ed90_predict(probit_model)
print(ed90_value)
