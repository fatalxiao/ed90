# 安装和加载所需的包
if (!require("drc")) {
  install.packages("drc")
}
library(drc)

# 假设你的数据如下:
doses <- c(0.1, 0.3, 1, 3, 10, 30) # 剂量
responses <- c(15, 20, 50, 80, 95, 98) # 对应剂量的反应

# 创建一个数据框
data <- data.frame(dose = doses, response = responses)

# 使用drc包的drm()函数拟合剂量-反应模型
# LL.4()表示4参数的对数-逻辑模型，其中参数为: b=c,b=d，e=c(下限),f=d(上限)
model <- drm(response ~ dose, data = data, fct = LL.4())

# 打印模型概要获取更多信息
summary(model)

# 使用ED()函数计算ED90
# 第一个参数是模型，第二个参数是你想要计算的药效比例（在这里是90）
ed90 <- ED(model, 90, interval = "delta")

# 显示ED90的值
print(ed90)
