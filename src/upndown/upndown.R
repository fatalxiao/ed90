# 安装和加载所需的包
if (!require("upndown")) {
    install.packages("upndown")
}
library(upndown)

# 创建数据框示例
groupS <- read.csv("./loading volume data.csv", 1, encoding = 'UTF-8')
groupS

# 准备数据
dose <- groupS$doseSequence
response <- groupS$responseSequence

ed90_estimate <- dixonmood(x = dose, y = response, full = TRUE, flip = FALSE)
print(ed90_estimate)
