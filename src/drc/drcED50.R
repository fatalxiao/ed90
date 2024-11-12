# 安装和加载所需的包
if (!require("drc")) {
    install.packages("drc")
}
library(drc)

groupS <- read.csv("./PIEB.csv", 1, encoding = 'UTF-8')
groupS

doses <- groupS$doseSequence
responses <- groupS$responseSequence

# 创建一个数据框
# data <- data.frame(response = responses, dose = doses)
data <- data.frame(response = c(rep(1, 8), 0, 1, 0, rep(1, 5), rep(0, 2), rep(1, 2), 0, 1, 0, rep(1, 3), 0, rep(1, 7), 0, rep(1, 3), rep(0, 2), rep(1, 8), 0, 1),
                   dose = c(rep(25, 7), rep(22.5, 2), rep(25, 2), rep(27.5, 6), 30, rep(32.5, 3), rep(35, 2), rep(37.5, 4), rep(40, 5), rep(37.5, 3), rep(40, 4), 42.5, rep(45, 9), 47.5))
data

# 使用drc包的drm()函数拟合剂量-反应模型
# LL.4()表示4参数的对数-逻辑模型，其中参数为: b=c,b=d，e=c(下限),f=d(上限)
model <- drm(data, fct = LL.4())
model

# 打印模型概要获取更多信息
summary(model)

# 使用ED()函数计算ED50
# 第一个参数是模型，第二个参数是你想要计算的药效比例（在这里是50）
ed90 <- ED(model, 90, interval = "delta")
ed90
